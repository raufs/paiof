import os
import sys
import logging
import traceback
import subprocess
import multiprocessing
from collections import defaultdict
from operator import itemgetter
import statistics

def createLoggerObject(log_file):
	"""
	Description:
	This function creates a logging object.
	********************************************************************************************************************
	Parameters:
	- log_file: Path to file to which to write logging.
	********************************************************************************************************************
	Returns:
	- logger: A logging object.
	********************************************************************************************************************
	"""

	logger = logging.getLogger('task_logger')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(log_file)
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	return logger

def closeLoggerObject(logObject):
	"""
	Description:
	This function closes a logging object.
	********************************************************************************************************************
	Parameters:
	- logObject: A logging object.
	********************************************************************************************************************
	"""

	handlers = logObject.handlers[:]
	for handler in handlers:
		handler.close()
		logObject.removeHandler(handler)

def setupDirectories(directories):
	"""
	Description:
	This is a generalizable function to create directories.
	********************************************************************************************************************
	Parameters:
	- dictionaries: A list of paths to directories to create or recreate (after removing).
	********************************************************************************************************************
	"""
	try:
		assert (type(directories) is list)
		for d in directories:
			if os.path.isdir(d):
				os.system('rm -rf %s' % d)
			os.system('mkdir %s' % d)
	except Exception as e:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)


def multiProcess(inputs):
	"""
	Genralizable function to be used with multiprocessing to parallelize list of commands. Inputs should correspond
	to space separated command (as list), with last item in list corresponding to a logging object handle for logging
	progress.
	"""
	input_cmd = inputs
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
	except Exception as e:
		sys.stderr.write(traceback.format_exc())

def setupDIAMONDdbs(og_seq_dir, diamond_db_dir, logObject, cpus):
	"""
	Function to setup individual DIAMOND dbs for self-blast later on.
	"""
	try:
		build_db_cmds = []
		for f in os.listdir(og_seq_dir):
			if f.endswith('.fa'): 
				og_seq = og_seq_dir + f
				diamond_db = diamond_db_dir + f + '.dmnd'
				diamond_db_cmd = ['diamond', 'makdedb', '--ignore-warnings', '--in', og_seq, '-d', diamond_db]
				build_db_cmds.append(diamond_db_cmd)
		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, build_db_cmds)
		p.close()
	except Exception as e:	
		sys.stderr.write('Problem with setting up DIAMOND dbs for orthogroups.\n')
		logObject.error('Problem with setting up DIAMOND dbs for orthogroups.')
		sys.stderr.write(traceback.format_exc())

def blastOGs(og_seq_dir, diamond_db_dir, blast_results_dir, logObject, cpus, evalue_cutoff=1e-3):
	"""
	Function to run self-blast.
	"""
	try:
		blastp_cmds = []
		for f in os.listdir(og_seq_dir):
			if f.endswith('.fa'): 
				og_seq = og_seq_dir + f
				diamond_db = diamond_db_dir + f + '.dmnd'
				results_file = blast_results_dir + f.split('.txt')
				assert (os.path.isfile(diamond_db) and os.path.isfile(og_seq))
				blastp_cmd = ['diamond', 'blastp', '--ignore-warnings', '--threads', '1', '--very-sensitive',
							  '--query', og_seq, '--db', diamond_db, '--outfmt', '6', 'qseqid', 'sseqid',
					   		  'pident', 'bitscore', '-k0', '--out', results_file, '--evalue', str(evalue_cutoff)]
				blastp_cmds.append(blastp_cmd)
		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, blastp_cmds)
		p.close()
	except Exception as e:	
		sys.stderr.write('Problem with performing self-BLASTp.\n')
		logObject.error('Problem with performing self-BLASTp.')
		sys.stderr.write(traceback.format_exc())

def calculatePairwiseAAI(blast_results_dir, species_id_file, sequence_id_file, final_results_tsv, logObject):
	"""
	Function to run self-blast.
	"""
	try:
		species_id_to_name = {}
		with open(species_id_file) as osif:
			for line in osif:
				line = line.strip()
				ls = line.split()
				sid = ls[0][:-1]
				sname = ls[1]
				species_id_to_name[sid] = sname

		sequence_name_to_species_id = {}
		genome_gene_counts = defaultdict(int)
		with open(sequence_id_file) as osif:
			for line in osif:
				line = line.strip()
				ls = line.split()
				sid = ls[0].split('_')[0]
				gene_name = ls[1]
				sequence_name_to_species_id[gene_name] = sid
				genome_gene_counts[sid] += 1

		genome_to_genome_aais = defaultdict(lambda: defaultdict(list))
		for f in os.listdir(blast_results_dir):
			blast_results_file = blast_results_dir + f

			results = []
			with open(blast_results_file) as obrf:
				for line in obrf:
					line = line.strip()
					query, subject, pident, bitscore = line.split('\t')
					query_sample = sequence_name_to_species_id[query]
					subject_sample = sequence_name_to_species_id[subject]
					results.append([query, subject, query_sample, subject_sample, float(bitscore), float(pident)])
			
			# The AAI is currently computed by taking the percent identity to the best match of a gene 
			# from sample/genome A to a gene in genome B belonging to the same ortholog group. This 
			# doesn't necessary indicate the two genes are RBHs. Best hit is determined based on bitscore.
			query_samples_already_matched = defaultdict(set)
			query_sample_genes_already_matched = defaultdict(set)
			for hit in sorted(results, key=itemgetter(4), reverse=True):
				qgene, sgene, q, s, _, pid = hit
				if s in query_samples_already_matched[qgene]: continue
				if sgene in query_sample_genes_already_matched[q]: continue
				genome_to_genome_aais[q][s].append(pid)
				query_samples_already_matched[qgene].add(s)
				query_sample_genes_already_matched[q].add(sgene)

		final_results_tsv_handle = open(final_results_tsv, 'w')
		final_results_tsv_handle.write('\t'.join(['Genome_1', 'Genome_2', 'AAI', 'Stdev_AAI', 'Genes_Shared', 'Genes_Shared_Normalized_by_Genome_1_Genes']) + '\n')
		for q in genome_to_genome_aais:
			query_name = species_id_to_name[q]
			query_genome_gene_count = genome_gene_counts[q]
			for s in genome_to_genome_aais:
				if q == s: continue
				subject_name = species_id_to_name[s]
				matches = genome_to_genome_aais[q][s]
				match_count = len(matches)
				match_prop = match_count/float(query_genome_gene_count)
				aai = statistics.mean(matches)
				stdev_aai = statistics.stdev(matches)
				printlist = [str(x) for x in [query_name, subject_name, aai, stdev_aai, match_count, match_prop]]
				final_results_tsv_handle.write('\t'.join(printlist) + '\n')
		final_results_tsv_handle.close()
		
	except Exception as e:
		sys.stderr.write('Problem with processing self-BLASTp results to determine pairwise AAIs.\n')
		logObject.error('Problem with processing self-BLASTp results to determine pairwise AAIs.')
		sys.stderr.write(traceback.format_exc())