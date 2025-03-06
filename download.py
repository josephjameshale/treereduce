
import subprocess
import os

def main(input_file,temp_dir,output_dir):
	# ensure temp dir is empty
	for f in os.listdir(temp_dir):
		os.remove(f)
	# generate list of runs to download
	run_list = make_run_list(input_file)
	for r in run_list:
		# download files
		# names should look like [run_name]_1.fastq
		command = ['fasterq-dump','--outdir',temp_dir,'--mem','5G','--split-3',
			 '--threads','1','--skip-technical',r]
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
			print(f'could not locate all files for {r}')

def move_and_compress(filename,input_dir,output_dir,read):
	# converts [run_name]_1.fastq to [run_name]_R1.fastq
	new_filename = filename.replace('_' + read + '.fastq','_R' + read + '.fastq')
	subprocess.run(['mv',input_dir + filename,output_dir + new_filename])
	subprocess.run(['gzip',output_dir + new_filename])

def make_run_list(input_file):
	l = []
	with open(input_file,'r') as fh:
		for line in fh:
			data_type,id,tree = line.strip().split('\t')
			if data_type == 'run':
				l.append(id)
			elif data_type == 'assembly':
				# this is where a function to convert assemblies to runs needs to go
				continue
	return(l)

if __name__ == "__main__":
	inp = '/scratch/esnitkin_root/esnitkin0/jjhale/treereduce/sra_run_list_t60.tsv'
	outp = '/scratch/esnitkin_root/esnitkin0/jjhale/public/c_auris_ncbi_pathogen/sra_t60_NoAssemblies_03-06-25/'
	temp = '/scratch/esnitkin_root/esnitkin0/jjhale/public/c_auris_ncbi_pathogen/temp/'
	main(input_file=inp,temp_dir=temp,output_dir=outp)

