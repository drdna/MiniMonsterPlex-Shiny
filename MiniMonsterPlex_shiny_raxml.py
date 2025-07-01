import os
import shutil
import subprocess
import glob
import gzip
from pathlib import Path

def fasta_filter(outPut,included_isolates,project_name):
	to_write = []
	with open(f'{outPut}/built_fasta/{project_name}builtSeqMeta.fasta','r') as read:
		lines = read.readlines()
		for i in range(0,len(lines)):
			if lines[i][0] == '>':
				print(lines[i].split('_')[0].split(">")[1].strip())
				if lines[i].split('_')[0].split(">")[1].strip() in included_isolates:
					to_write.append([lines[i],lines[i+1]])
	with open(f'{outPut}/built_fasta/{project_name}builtSeqFiltered.fasta','a') as write:
		for isolate in to_write:
			write.write(f'{isolate[0]}{isolate[1]}')


def autoRAxML(outPut,filtered,project_name):
	wd = os.getcwd()
	#command for running RAXML
	command = ['raxmlHPC',
		'-w',
		os.path.join(wd,outPut),
		'-p', '1234',
		'-f', 'a',
		'-x', '1234',
		'-s', f'{outPut}/built_fasta/{project_name}builtSeqFiltered.fasta' if filtered else f'{outPut}/built_fasta/{project_name}builtSeqMeta.fasta',
		'-n', 'miniMonsterPlex.raxml',
		'-m', 'GTRGAMMA',
		'-#', '1000']
	try:
		subprocess.run(' '.join(command),
				shell=True,
				check=True,
				capture_output=True,
				text=True)
	except subprocess.CalledProcessError as e:
		raxmls = glob.glob(os.path.join(outPut,'*.raxml'))
		for file in raxmls:
			os.remove(file)
		error_details = (f"Something went wrong with RAxML.\n"
						 f"Command: {' '.join(command)}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)
	
	os.makedirs(os.path.join(outPut,'raxml_out'), exist_ok=True)
	raxmls = glob.glob(os.path.join(outPut,'*miniMonsterPlex.raxml*'))
	for file in raxmls:
		shutil.move(file,os.path.join(outPut,'raxml_out', Path(file).name))
		
def raxmlGate(outPut_Folder,filtered, project_name):
	if filtered:
		fasta_path = f'{outPut_Folder}/built_fasta/{project_name}builtSeqFiltered.fasta'
		with open(fasta_path, 'r') as read:
			cnt = sum(1 for line in read if line.startswith(">"))
		if cnt < 4:
			print("Tree building requires a minimum of 4 isolates, but you have fewer than 4 after filtering.")
			print("RAxML will not be run on the filtered set.")
			return None # Special value to indicate skipping RAxML
		else:
			return True # Proceed with filtered
	else:
		return False # Proceed with unfiltered

def mlTree(outPut_Folder,project_name):
	try:
		command = ['Rscript',
			'--vanilla',
			'MLtree.R',
			f'{outPut_Folder}/raxml_out/RAxML_bestTree.miniMonsterPlex.raxml']
		subprocess.run(' '.join(command),
					shell=True,
					check=True,
					capture_output=True,
					text=True)
	except subprocess.CalledProcessError as e:
		if os.path.isfile('NA.pdf'):
			os.remove('NA.pdf')
		error_details = (f"Something went wrong with the MLtree.R script.\n"
						 f"Command: {e}\n"
						 f"Return Code: {e.returncode}\n"
						 f"Stderr: {e.stderr}")
		raise RuntimeError(error_details)
	os.makedirs(os.path.join(outPut_Folder,'tree_out'), exist_ok=True)
	shutil.move('NA.pdf',f'{outPut_Folder}/tree_out/{project_name}_tree.pdf')

#series of lines for cleaing up left over temp data
def cleanup(outPut,input_folder,complete,project_name):
	print(input_folder)
	print(os.path.join('Projects',project_name,'completed_fastq/'))
	files_to_move = glob.glob(os.path.join(input_folder,'*.gz'))
	for file in files_to_move:
		shutil.move(file,os.path.join('Projects',project_name,'completed_fastq/'))

	with open('totalMergedCall.vcf', 'a') as f:
		with open(f'{outPut}/merge_out/{project_name}MergedCallAll.vcf','r') as read:
			for line in read:
				f.write(line)

	with open(f'{outPut}/merge_out/{project_name}MergedCallAll.vcf', 'rb') as f_in:
		with gzip.open(f'{outPut}/merge_out/{project_name}MergedCallAll.vcf.gz', 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
	shutil.move(f'{outPut}/merge_out/{project_name}MergedCallAll.vcf.gz',os.path.join('Projects',project_name,'processed_vcf/'))

	shutil.rmtree(f'{outPut}/bowtie_out/')
	shutil.rmtree(f'{outPut}/coverage_out/')
	shutil.rmtree(f'{outPut}/call_out/')
	shutil.rmtree(f'{outPut}/mpileup_out/')
	shutil.rmtree(f'{outPut}/merge_out/')
	os.remove(f'{outPut}/fastqListCall.txt')

	
	# with open('totalFasta.mfa','a') as f:
	# 	with open(f'{outPut}/built_fasta/{project_name}builtSeqMeta.fasta','r') as read:
	# 		for line in read:
	# 			f.write(line)
	# #does the combined anyalysis of all the data
	# if complete:
	# 	command = ['raxmlHPC',
	# 		 '-p',
	# 		 '1234',
	# 		 '-f',
	# 		 'a',
	# 		 '-x',
	# 		 '1234',
	# 		 '-s',
	# 		 'totalFasta.mfa',
	# 		 '-n',
	# 		 'miniMonsterPlex_full.raxml',
	# 		 '-m',
	# 		 'GTRGAMMA',
	# 		 '-#',
	# 		 '1000']
	# 	subprocess.run(' '.join(command),
	# 			 shell=True,
	# 			 check=True)
	# 	command = ['Rscript',
	# 	   '--vanilla',
	# 	   'MLtree.R',
	# 	   'RAxML_bestTree.miniMonsterPlex_full.raxml']
	# 	subprocess.run(' '.join(command),
	# 			shell=True,
	# 			check=True)
	# 	files_to_delete = glob.glob(os.path.join(os.getcwd(),'*.raxml'))
	# 	for file in files_to_delete:
	# 		os.remove(file)
	# 	files_to_delete = glob.glob(os.path.join(os.getcwd(),'*.raxml.support'))
	# 	for file in files_to_delete:
	# 		os.remove(file)
def main(project, included_isolates,complete=False):
	outPut_Folder = os.path.join('Projects',project,"output")
	input_folder = os.path.join('Projects',project,"newFastq")

	filtered = False
	fileList = glob.glob(f'{input_folder}/*.gz')


	if len(included_isolates) >= 1:
		fasta_filter(outPut_Folder, included_isolates,project)
		filtered = True

	raxml_decision = raxmlGate(outPut_Folder,filtered, project)

	if raxml_decision is not None and not os.path.isdir(os.path.join(outPut_Folder,'raxml_out')):
		autoRAxML(outPut_Folder, raxml_decision, project)

	if os.path.exists(f'{outPut_Folder}/raxml_out/RAxML_bestTree.miniMonsterPlex.raxml'):
		mlTree(outPut_Folder,project)

	cleanup(outPut_Folder,input_folder,complete,project)