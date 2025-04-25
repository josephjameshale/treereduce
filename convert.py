
import random
import subprocess
import os


def main(input_file,output_file,conversion_table):
	cdict = make_conv(conversion_table)
	with open(input_file,'r') as fh_in,open(output_file,'w') as fh_out:
		for line in fh_in:
			isolate,tree = line.strip().split('\t')
			run_data = cdict[isolate] + [tree]
			_ = fh_out.write('\t'.join(run_data) + '\n')
			if run_data[1] == '':
				print(run_data)

def make_conv(input_file):
	od = {}
	with open(input_file,'r') as fh:
		header = fh.readline()
		isolate_pos,run_pos,assembly_pos = find_header_index(header,'Isolate'),find_header_index(header,'Run'),find_header_index(header,'Assembly')
		for line in fh:
			line2 = line[:-1].split('\t')
			isolate,run,assembly = line2[isolate_pos],line2[run_pos],line2[assembly_pos]
			if run != '':
				od[isolate] = ['run',run]
			else:
				od[isolate] = ['assembly',assembly]
	return(od)

def find_header_index(list_inp,search_pattern):
	hits = [ind for ind,st in enumerate(list_inp.strip().split('\t')) if st == search_pattern]
	if len(hits) == 1:
		return(hits[0])
	else:
		print(list_inp)
		quit(1)

if __name__ == "__main__":
	inp = 'results/results_t30/summary/isolate_list_t30.tsv'
	outp = 'sra_run_list_t30.tsv'
	conv = 'c_auris_isolates_browser_all_02-17-25.tsv'
	main(input_file=inp,output_file=outp,conversion_table=conv)

