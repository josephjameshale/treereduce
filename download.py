import argparse
import subprocess
import os
import re

def master(input_file,temp_dir,output_dir,finished_list):
	# ensure temp dir is empty
	for f in os.listdir(temp_dir):
		subprocess.run(['rm','-r',temp_dir + f])
		#os.remove(temp_dir + f)
	# generate list of runs to download
	run_list = make_run_list(input_file)
	run_list = [(run_name,data_type) for run_name,data_type in run_list if run_name not in finished_list]
	# this removes any runs that were already processed
	for run_name,data_type in run_list:
		if data_type == 'run':
			download_run(run_name,temp_dir,output_dir)
		if data_type == 'assembly':
			download_assembly(run_name,temp_dir,output_dir + 'assemblies/')

def download_assembly(run_name,temp_dir,output_dir):
	temp_assembly_dir = temp_dir + run_name + '/'
	temp_assembly_zip = temp_assembly_dir + run_name + '.zip'
	final_assembly_dir = temp_assembly_dir + 'ncbi_dataset/data/' + run_name + '/'
	commandlist = []
	commandlist.append(['mkdir',temp_assembly_dir])
	# this makes a separate directory in the temp folder for this assembly
	commandlist.append(['datasets','download','genome','accession',run_name,'--include', 
		'gff3,rna,cds,protein,genome,seq-report','--filename',temp_assembly_zip])
	# this downloads the assembly as a single zip file in this directory
	commandlist.append(['unzip',temp_assembly_zip,'-d',temp_assembly_dir])
	# this unzips all contents into the assembly dir
	commandlist.append(['mv',final_assembly_dir,output_dir])
	# this moves the entire data folder to the output directory, keeping its original name
	commandlist.append(['mv',temp_assembly_dir + 'ncbi_dataset/data/assembly_data_report.jsonl',output_dir + f'{run_name}/'])
	commandlist.append(['mv',temp_assembly_dir + 'ncbi_dataset/data/dataset_catalog.json',output_dir + f'{run_name}/'])
	# this moves the two metadata files into the assembly directory
	for c in commandlist:
		subprocess.run(c)

def download_run(run_name,temp_dir,output_dir):
	# download files
	# names should look like [run_name]_1.fastq
	command = ['fasterq-dump','--outdir',temp_dir,'--mem','5G','--split-3',
		 '--threads','1','--skip-technical',run_name]
	subprocess.run(command)
	forward_file_missing,reverse_file_missing = True,True
	for f in os.listdir(temp_dir):
		if '_1.fastq' in f:
			move_and_compress(f,temp_dir,output_dir,read='1')
			forward_file_missing = False
		elif '_2.fastq' in f:
			move_and_compress(f,temp_dir,output_dir,read='2')
			reverse_file_missing = False
	if any([forward_file_missing,reverse_file_missing]):
		print(f'could not locate all files for {run_name}')

def move_and_compress(filename,input_dir,output_dir,read):
	# converts [run_name]_1.fastq to [run_name]_R1.fastq
	new_filename = filename.replace('_' + read + '.fastq','_R' + read + '.fastq')
	subprocess.run(['mv',input_dir + filename,output_dir + new_filename])
	subprocess.run(['gzip',output_dir + new_filename])

def make_run_list(input_file):
	run_list = []
	with open(input_file,'r') as fh:
		for line in fh:
			data_type,id,tree = line.strip().split('\t')
			run_list.append((id,data_type))
			# data_type should be either 'run' or 'assembly'
			# if data_type == 'run':
			# 	run_list.append((id,'run'))
			# elif data_type == 'assembly':
			# 	assembly_list.append(id)
	return(run_list)

def make_finished_list(dirlist):
	# this looks at multiple directories and returns a list of the finished accession numbers present in them
	finished_list = set()
	if dirlist is not None:
		for dirpath in dirlist:
			finished_list.update(check_run_list(dirpath))
	return(finished_list)

def check_run_list(raw_data_dir):
	finished_list = set()
	for filename in os.listdir(raw_data_dir):
		if '.fastq.gz' in filename:
			filename2 = re.split('_R[1,2].fastq.gz',filename)[0]
			finished_list.add(filename2)
		elif os.path.isdir(raw_data_dir + filename):
			finished_list.add(filename)
	return(finished_list)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'--input','-i',type=str,
		help='''Provide a one-column tab-separated file consisting of run types, accession numbers, and phylogeny trees. See the provided example for formatting.''',
		required=True
		)
	parser.add_argument(
		'--temp','-t',type=str,
		help='''Provide a path to a temporary directory.''',
		default='./temp/'
		)
	parser.add_argument(
		'--output_dir','-od',type=str,
		help='''Provide a path to a directory where the downloaded files will be stored.''',
		required=True
		)
	parser.add_argument(
		'--finished_directory_list','-fdl',type=str,nargs='+',
		help='''Provide a list of directories containing raw reads from samples that you do NOT want to download again. Any named subdirectory in this folder will
		also be skipped, such as finished assembly downloads.''',
		default=None
		)
	args = parser.parse_args()
	output_dir_assembly = args.output_dir + 'assemblies/'
	for p in [args.output_dir,args.temp,output_dir_assembly]:
		if p is not None and not os.path.isdir(p):
			subprocess.run(['mkdir','-p',p])
	finished_name_list = make_finished_list(args.finished_directory_list)
	master(input_file=args.input,temp_dir=args.temp,output_dir=args.output_dir,finished_list = finished_name_list)

if __name__ == "__main__":
	main()
