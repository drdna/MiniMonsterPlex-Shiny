# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 13:49:59 2023

@author: treyd
"""
#importing packages
import os
import glob
import re
import subprocess 
import multiprocessing
from pathlib import Path

def auto_bowtie2(outPut_Folder, input_file, fileNum,threads):
	try:
		print(fileNum, " is entering the pipeline")
		# bowtie2 + samtools sort call
		command = ['bowtie2', '--no-unal',
			'-p', str(threads),
			'-x', 'index/70-15_small_index',
			'-U', input_file,
			'--local', '--very-sensitive-local',
			'2>', f'{outPut_Folder}/bowtie_out/{fileNum}_alignment_summary.txt',
			'|',
			'samtools', 'sort',
			'-',
			'-@', '2',
			'-O', 'bam',
			'-o', f'{outPut_Folder}/bowtie_out/{fileNum}hits.bam']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True,
					text=True)
	except subprocess.CalledProcessError as e:
		# Clean up failed files
		if os.path.isfile(f'{outPut_Folder}/bowtie_out/{fileNum}hits.bam'):
			os.remove(f'{outPut_Folder}/bowtie_out/{fileNum}hits.bam')
		if os.path.isfile(f'{outPut_Folder}/bowtie_out/{fileNum}_alignment_summary.txt'):
			os.remove(f'{outPut_Folder}/bowtie_out/{fileNum}_alignment_summary.txt')
		
		# Raise an exception with detailed error info instead of quitting
		error_details = (f"Something went wrong with bowtie2 for file: {fileNum}.\n"
						 f"Command: {e}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)


def parse_alignment_summary(fileNum, outPut_Folder):
	with open(fileNum, 'r') as f:
		content = f.read()

	# Extract total reads
	total_match = re.search(r'(\d[\d,]*) reads; of these:', content)
	if not total_match:
		raise ValueError(f"Total reads not found in {fileNum}")
	total_reads = int(total_match.group(1).replace(',', ''))

	# Extract aligned 1 time
	aligned_once_match = re.search(r'(\d[\d,]*) \([\d.]+%\) aligned exactly 1 time', content)
	if not aligned_once_match:
		raise ValueError(f"Aligned exactly 1 time not found in {fileNum}")
	aligned_once = int(aligned_once_match.group(1).replace(',', ''))

	# Extract aligned >1 times
	aligned_multiple_match = re.search(r'(\d[\d,]*) \([\d.]+%\) aligned >1 times', content)
	if not aligned_multiple_match:
		raise ValueError(f"Aligned >1 times not found in {fileNum}")
	aligned_multiple = int(aligned_multiple_match.group(1).replace(',', ''))

	# Compute total aligned
	aligned_total = aligned_once + aligned_multiple
	aligned_fraction = aligned_total / total_reads

	# Return CSV string
	return f"{Path(fileNum).stem.replace('_alignment_summary','')},{total_reads},{aligned_total},{aligned_fraction:.4f}"


def auto_mpileup(outPut, fileNum, threads):
	try:
		command = ['bcftools',
				'mpileup',
				'--threads',
				str(threads),
				'-d',
				'100000',
				#'-R',
				#'MonsterPlexRegionsFileSuperCont_small_index.txt',
				'--annotate',
				'FORMAT/AD',
				'-f',
				'index/70-15_small.fasta',
				os.path.join(outPut,'bowtie_out',f'{fileNum}hits.bam'),
				'>>',
				f'{outPut}/mpileup_out/{fileNum}.vcf']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True,
					text=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(f'{outPut}/mpileup_out/{fileNum}.vcf'):
			os.remove(f'{outPut}/mpileup_out/{fileNum}.vcf')
		error_details = (f"Something went wrong with mpileup for file: {fileNum}.\n"
						 f"Command: {e}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)

	
def auto_call(outPut, fileNum):
	try:
		command = ['bcftools',
				'call',
				'-c',
				'--ploidy',
				'1',
				os.path.join(outPut, 'mpileup_out', f'{fileNum}.vcf'),
				'-o',
				f'{outPut}/call_out/{fileNum}call.vcf']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True,
					text=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(f'{outPut}/call_out/{fileNum}call.vcf'):
			os.remove(f'{outPut}/call_out/{fileNum}call.vcf')
		error_details = (f"Something went wrong with bcftools call for file: {fileNum}.\n"
						 f"Command: {e}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)

	
def auto_bedtools(outPut, fileNum):
	try:
		command = ['bedtools',
				'genomecov',
				'-ibam',
				os.path.join(outPut, 'bowtie_out', f'{fileNum}hits.bam'),
				'-bg',
				'>',
				f'{outPut}/coverage_out/{fileNum}cover.bed']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True,
					text=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(f'{outPut}/coverage_out/{fileNum}cover.bed'):
			os.remove(f'{outPut}/coverage_out/{fileNum}cover.bed')
		error_details = (f"Something went wrong with bedtools genomecov for file: {fileNum}.\n"
						 f"Command: {e}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)

def auto_bgzip(outPut, fileNum,file_ex):
	try:
		command =['bgzip',
				os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}')]
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True,
					text=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}.gz')):
			os.remove(os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}.gz'))
		error_details = (f"Something went wrong with bgzip for file: {fileNum}.\n"
						 f"Command: {e}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)

	try:
		command = ['tabix',
				os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}.gz')]
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True,
					text=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}.gz.tbi')):
			os.remove(os.path.join(outPut,'mpileup_out',f'{fileNum}{file_ex}.gz.tbi'))
		error_details = (f"Something went wrong with tabix for file: {fileNum}.\n"
						 f"Command: {e}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)

def autoVCFZip(outPut, file, fileNum):
	#bg zip the bcftools call result file
	try:
		command = ['bgzip',
				os.path.join(outPut,'call_out',f'{fileNum}call.vcf')]
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True,
					text=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(os.path.join(outPut,'call_out',f'{fileNum}call.vcf.gz')):
			os.remove(os.path.join(outPut,'call_out',f'{fileNum}call.vcf.gz'))
		error_details = (f"Something went wrong with bgzip for file: {fileNum}.\n"
						 f"Command: {e}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)

	#tabix the call results
	try:
		command = ['tabix',
				os.path.join(outPut,'call_out',f'{fileNum}call.vcf.gz')]
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True,
					text=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(os.path.join(outPut,'call_out',f'{fileNum}call.vcf.gz.tbi')):
			os.remove(os.path.join(outPut,'call_out',f'{fileNum}call.vcf.gz.tbi'))
		error_details = (f"Something went wrong with tabix for file: {fileNum}.\n"
						 f"Command: {' '.join(command)}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)
	with open(f'{outPut}/fastqListCall.txt', 'a') as append:
		file_name = Path(file)
		file_name = file_name.name.split('.')[0]
		append.write(f'{outPut}/call_out/{file_name}call.vcf.gz\n')

def autoMerge(outPut,project_name):
	try:
		command = ['bcftools',
				'merge',
				'-l',
				f'{outPut}/fastqListCall.txt',
				'-o',
				f'{outPut}/merge_out/{project_name}MergedCallAll.vcf']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True,
					text=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile(os.path.join(outPut,'merge_out',f'{project_name}MergedCallAll.vcf')):
			os.remove(os.path.join(outPut,'merge_out',f'{project_name}MergedCallAll.vcf'))
		error_details = (f"Something went wrong with bcftools merge.\n"
						 f"Command: {e}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)
		
def sampleBuilder(outPut,metadata_file,project_name):
	sites =[]
	sitesUsed =[]
	#reads a list of sites you want and only looks at data from there
	with open('MonsterPlexSitesList_small_index.txt', 'r') as read:
		for line in read:
			sites.append(line.strip('\n'))

	with open(f'{outPut}/merge_out/{project_name}MergedCallAll.vcf', 'r') as read:
		seqs = list()
		check = False
		for line in read:
			if line.split('\t')[0] == '#CHROM':
				print("header past seqs made")
				fqList = line.split('\t')
				for n in range(9,len(fqList)):
					seqs.append([fqList[n], ''])
					check = True
			#elif check:
			elif check and line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1] in sites:
				#this creates a horizontal split of the line
				lineList = line.strip('\n').split('\t')
				for n in range(9,len(lineList)):
					fields = lineList[n].split(':')
					if fields[0] == '.':
						seqs[n - 9][1] += "N"
					elif len(fields[2].split(',')) == 1:
						if fields[0] == '0':
							if int(fields[2]) > 5:
								seqs[n - 9][1] += lineList[3]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							#this checks alt
						elif fields[0] == '1':
							if int(fields[2]) > 5:
								seqs[n - 9][1] += lineList[4]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						else:
							raise RuntimeError(f"Something went wrong with {seqs[n][0]}\n{lineList[n]} at site {line.strip().split()[0]} {line.strip().split()[1]}")
					#this checks cases were both ref and alt are registered
					elif len(fields[2].split(',')) >= 2:
						#this creates a list out of the AD field
						AD = fields[2].split(',')
						#this checks if ref is blank
						if AD[0] == '.':
							if int(AD[1]) > 5:
								seqs[n - 9][1] += lineList[4].split(',')[0]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						#this checks if alt is blank
						elif AD[1] == '.':
							if int(AD[0]) > 5:
								seqs[n - 9][1] += lineList[3]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						#checks if ref is greater then alt
						elif int(AD[0]) > int(AD[1]):
							if int(AD[0]) > (int(AD[1]) * 20):
								seqs[n - 9][1] += lineList[3]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						#checks if alt is greater than ref
						elif int(AD[1]) > int(AD[0]):
							if int(AD[1]) > (int(AD[0]) * 20):
								seqs[n - 9][1] += lineList[4].split(',')[0]
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
							else:
								seqs[n - 9][1] += "N"
								sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						elif int(AD[1]) == int(AD[0]):
							seqs[n - 9][1] += lineList[3]
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						else:
							raise RuntimeError(f"Something went wrong with {seqs[n][0]}\n{lineList[n]} at site {line.strip().split()[0]} {line.strip().split()[1]}")
					else:
						raise RuntimeError(f"Something went wrong with {seqs[n][0]}\n{lineList[n]} at site {line.strip().split()[0]} {line.strip().split()[1]}")

		sample_metadata = metaDataBuilder(metadata_file)
		
		print(sample_metadata)
		
		os.makedirs(f'{outPut}/built_fasta', exist_ok=True)
		
		with open(f'{outPut}/built_fasta/{project_name}builtSeqMeta.fasta', 'a') as writeSeq:
			for read in seqs:
				file_path = Path(read[0])
				seqID = file_path.name.split('.')[0].split('hits')[0]
				if len(seqID.split("_")) > 1:
					seqID = f'{"-".join(seqID.split("_"))}'
				if (seqID) in sample_metadata:
					seqSpecies = sample_metadata[seqID][0]
					seqHost = sample_metadata[seqID][1]
					seqLineage = sample_metadata[seqID][2]
					seqCountry = sample_metadata[seqID][3]
					writeSeq.write(f'>{seqID}_{seqSpecies}_{seqHost}_{seqLineage}_{seqCountry}\n{read[1]}\n')
				else:
					writeSeq.write('>' + seqID
								+ '_._._._.' + '\n' + read[1] + '\n')
						
def metaDataBuilder(metadata_file):
	metaData = {}
	with open(metadata_file, 'r') as read:
		for line in read:
			ID = line.split(',')[0].strip('\n')
			species = line.split(',')[1].strip('\n')
			host = line.split(',')[2].strip('\n')
			lineage = line.split(',')[3].strip('\n')
			country = line.split(',')[4].strip('\n')
			metaData[ID] = [species, host, lineage, country]
		return metaData


def main(project, metadata_file, complete=False):
	outPut_Folder = os.path.join('Projects',project,"output")
	metadata_file_name = os.path.join('Projects',project,"metadata",metadata_file)
	input_folder = os.path.join('Projects',project,"newFastq")
	threads = multiprocessing.cpu_count()
	if threads > 8:
		threads = 8

	filtered = False
	fileList = glob.glob(f'{input_folder}/*.gz')

	print(f"{'Variable':<20} {'Value'}")
	print(f"{'-'*20} {'-'*40}")
	print(f"{'Project':<20} {project}")
	print(f"{'Metadata File':<20} {metadata_file_name}")
	print(f"{'Output Folder':<20} {outPut_Folder}")
	print(f"{'Input Folder':<20} {input_folder}")
	print(f"{'Complete mode':<20} {complete}")

	os.makedirs(outPut_Folder,exist_ok=True)
	os.makedirs(os.path.join('Projects',project,'completed_fastq'),exist_ok=True)
	os.makedirs(os.path.join('Projects',project,'processed_vcf/'),exist_ok=True)
	os.makedirs(os.path.join(outPut_Folder,'bowtie_out'),exist_ok=True)
	os.makedirs(os.path.join(outPut_Folder,'mpileup_out'),exist_ok=True)
	os.makedirs(os.path.join(outPut_Folder, "processedAlignSumm"),exist_ok=True)
	os.makedirs(os.path.join(outPut_Folder,'call_out'),exist_ok=True)
	os.makedirs(os.path.join(outPut_Folder,'coverage_out'),exist_ok=True)
	os.makedirs(os.path.join(outPut_Folder, 'merge_out'),exist_ok=True)

	alignment_summary_file = f'{outPut_Folder}/alignment_summary.csv'
	if os.path.exists(alignment_summary_file):
		os.remove(alignment_summary_file)

	file_name_list = []
	#everything now runs on a file by file basis skipping steps if an output file already exists
	for file in fileList:
		file_path = Path(file)
		fileNum = file_path.name.split('.')[0]
		file_name_list.append(fileNum)
		if not os.path.isfile(os.path.join(outPut_Folder,'bowtie_out',f'{fileNum}hits.bam')):
			auto_bowtie2(outPut_Folder, file,fileNum, threads)

		summary_string = parse_alignment_summary(f'{outPut_Folder}/bowtie_out/{fileNum}_alignment_summary.txt', outPut_Folder)
		with open(alignment_summary_file, "a") as out_file:
			out_file.write(summary_string + "\n")

		if not os.path.isfile(os.path.join(outPut_Folder,'bowtie_out',f'{fileNum}hits.bam.csi')):
			try:
				command = ['samtools', 'index', os.path.join(outPut_Folder,'bowtie_out',f'{fileNum}hits.bam')]
				subprocess.run(command, check=True, capture_output=True, text=True)
			except subprocess.CalledProcessError as e:
				error_details = (f"Something went wrong with samtools index for file: {fileNum}.\n"
								 f"Command: {e}\n"
								 f"Return Code: {e.returncode}\n"
								 f"Stderr: {e.stderr}")
				raise RuntimeError(error_details)

		if not os.path.isfile(os.path.join(outPut_Folder,'mpileup_out',f'{fileNum}.vcf')):
			auto_mpileup(outPut_Folder, fileNum, threads)
		
		if not os.path.isfile(os.path.join(outPut_Folder,'call_out',f'{fileNum}call.vcf')):
			auto_call(outPut_Folder, fileNum)

		if not os.path.isfile(f'{outPut_Folder}/coverage_out/{fileNum}cover.bed'):
			auto_bedtools(outPut_Folder, fileNum)

		file_extension = ".vcf"
		if not os.path.isfile(os.path.join(outPut_Folder,'mpileup_out',f'{fileNum}{file_extension}.gz.tbi')):
			auto_bgzip(outPut_Folder, fileNum, file_extension)

		if not os.path.isfile(os.path.join(outPut_Folder,'call_out',f'{fileNum}call.vcf.gz.tbi')):	
			autoVCFZip(outPut_Folder, file, fileNum)

	if not os.path.isfile(os.path.join(outPut_Folder,'merge_out',f'{project}MergedCallAll.vcf')):
		autoMerge(outPut_Folder,project)

	if not os.path.isdir(os.path.join(outPut_Folder,'built_fasta')):
		sampleBuilder(outPut_Folder,metadata_file_name,project)
	
	return file_name_list

