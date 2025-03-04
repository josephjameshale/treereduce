
import random
import subprocess
import os

# main still needs a report or output about how many tips were removed from each tree

def main(input_dir,output_dir,thresh,rseed):
	flist = []
	# list all files with the expected format and give them each a prefix
	for f in os.listdir(input_dir):
		if '.newick_tree.newick' in f:
			flist.append((input_dir + f,f.split('.newick_tree.newick')[0]))
	with open(output_dir + 'summary/isolate_list_t' + str(thresh) + '.tsv','w') as fh_final, open(output_dir + 'summary/summary.tsv','w') as fh_sum:
		fh_sum.write('snp_cluster\toriginal_tree_size\tgotree_tips_removed\tnew_tree_size\tisolates_selected\n')
		for treepath,prefix in flist:
			prunedmatrix_path,original_tree_size,tips_removed,new_tree_size,FIX_FLAG,firstname = reduce(treepath,output_dir,thresh,prefix,rseed)
			isolate_list = make_isolate_list(prunedmatrix_path,FIX_FLAG,firstname,thresh)
			for i in isolate_list:
				_ = fh_final.write('\t'.join([i,prefix]) +'\n')
			data = [str(x) for x in [prefix,original_tree_size,tips_removed,new_tree_size,len(isolate_list)]]
			_ = fh_sum.write('\t'.join(data) + '\n')

def make_isolate_list(matrix,FIX_FLAG,firstname,thresh):
	# this makes a list of isolates from the provided distance matrix
	# if the distance matrix is empty, assume that gotree incorrectly pruned all names and return firstname
	# if FIX_FLAG is True, remove the first instance of a name with a small distance (there should only be one pair of these)
	isolate_list = []
	with open(matrix,'r') as fh:
		next(fh)
		for line in fh:
			d = line.strip().split('\t')
			# name of the current isolate
			tipname = d[0]
			distances = [float(x) for x in d[1:]]
			if FIX_FLAG and any([x < thresh for x in distances]):
				# if this file was tagged for fixing, skip the first name where a distance below the threshold is found
				FIX_FLAG = False
				continue
			if tipname != '' and float_check(tipname) is False:
				isolate_list.append(tipname)
	if isolate_list == []:
		# in some cases, the distance matrix will not have any names in it 
		isolate_list.append(firstname)
	return(isolate_list)

def float_check(inp_string):
	try:
		v = float(inp_string)
	except ValueError:
		return(False)
	return(True)


def reduce(treepath,output_dir,thresh,prefix,rseed):
	# remove quotes
	cleantree_path = output_dir + 'temp/' + prefix + '_cleaned_tree.newick'
	cleantree(treepath,cleantree_path)
	# generate distance matrix (this assumes your tree contains at least 2 nodes)
	matrix_path = output_dir + 'temp/' + prefix + '_distances.tsv'
	make_matrix(cleantree_path,matrix_path)
	# determine names to remove
	to_prune,matrix_size,firstname = find_prune(matrix_path,thresh,rseed)
	# remove these names from the matrix (with exceptions made for subsetting to exactly 2 tips)
	prunedtree_path = output_dir + prefix + '_pruned_cleaned_tree.newick'
	original_tree_size,tips_removed,new_tree_size,FIX_FLAG = gotree_prune(cleantree_path,matrix_size,to_prune,prunedtree_path)
	# recalculate distance matrix
	prunedmatrix_path = output_dir + prefix + '_pruned_distances.tsv'
	make_matrix(prunedtree_path,prunedmatrix_path,treesize=new_tree_size)
	return([prunedmatrix_path,original_tree_size,tips_removed,new_tree_size,FIX_FLAG,firstname])

def gotree_prune(treepath,matrix_size,to_prune,prunedtree_path):
	# given a tree and a list of names, use gotree to remove those tips from the tree
	# note that gotree has two limitations: 
	# 1) it cannot remove the first name in a small tree
	# 2) it cannot prune a tree down to 2 tips
	# limitation 1 should not be a problem, as the first name is never added to the list of names to prune
	# address limitation 2 first, randomly removing a name if it will result is a two-node tree (this may result in to_prune being empty)
	FIX_FLAG = False
	if len(to_prune) > 0 and matrix_size - len(to_prune) == 2:
		# remove a random element from the list and flag this to be fixed later
		_ = to_prune.pop(random.randint(0,len(to_prune)-1))
		FIX_FLAG = True
	# return the same tree if there is nothing to prune
	if len(to_prune) == 0:
		command = ['cp',treepath,prunedtree_path]
		print(' '.join(command))
		subprocess.run(command)
	else:
		command = ['singularity','exec','gotree_0-4-4.sif','gotree','prune','-i',treepath,'-o',prunedtree_path] + to_prune
		print(' '.join(command))
		subprocess.run(command)
	original_tree_size = matrix_size
	tips_removed = len(to_prune)
	new_tree_size = matrix_size - len(to_prune)
	return([original_tree_size,tips_removed,new_tree_size,FIX_FLAG])

def cleantree(filepath,outpath):
	# removes troublesome quotation marks from names
	with open(filepath,'r') as fhin, open(outpath,'w') as fhout:
		for line in fhin:
			fhout.write(line.replace("'",""))

def make_matrix(treepath,outpath,treesize = 2):
	if treesize > 1:
		command = ['singularity','exec','gotree_0-4-4.sif','gotree','matrix','-i',treepath,'-o',outpath,'--format','newick']
		print(' '.join(command))
		subprocess.run(command)
	else:
		# gotree cannot generate a matrix when the tree has only one node, so this makes the matrix manually
		with open(treepath,'r') as fhin, open(outpath,'w') as fhout:
			name1 = fhin.readline().strip().split(';')[0]
			fhout.write('1\n')
			fhout.write(name1 + '\t0.000000000000\n')

def find_prune(matrix_path,thresh,rseed):
	# this returns a list of the names to prune from the provided distance matrix, which can be empty
	# it also returns the size of the distance matrix
	# make a dictionary of index->isolate name, identify first name
	indict,firstname = make_indict(matrix_path)
	to_prune = set()
	rm_lines = random_matrix_lines(matrix_path,rseed)
	# this takes random lines from the distance matrix
	# this can only be done because indict is NOT randomized
	for line in rm_lines:
		d = line.strip().split('\t')
		# name of the current isolate
		tipname = d[0]
		if tipname in to_prune:
			# if this is an isolate that is already marked for removal, do not check its similarity with others
			continue
		# distances between the current isolate and each other (including itself)
		distances = [float(x) for x in d[1:]]
		# return the name of each isolate with a distance below the threshold, unless this is the comparison to itself
		# indict[dist_index] gives the isolate name that corresponds to that position in the vector
		close_names = [indict[dist_index] for dist_index,dist_val in enumerate(distances) if dist_val < thresh and indict[dist_index] != tipname]
		to_prune.update(close_names)
	return(list(to_prune),len(indict),firstname)

def random_matrix_lines(filepath,rseed):
	with open(filepath,'r') as fh:
		next(fh)
		rlines = fh.readlines()
	random.seed(rseed)
	random.shuffle(rlines)
	return(rlines)


def make_indict(filepath):
	# this returns a dictionary where the index is the key and the name of the isolate is the value
	# this also returns the name of the first isolate in the list, which should correspond to the leftmost isolate in the tree
	# this is only done to avoid the gotree issues for small tree sizes, where the leftmost name can't be removed
	d,firstname = {},''
	with open(filepath,'r') as fh:
		next(fh)
		c = 0
		for line in fh:
			isolate_name = line.strip().split('\t')[0]
			d[c] = isolate_name
			c+=1
	firstname = random.choice(list(d.values()))
	# this selects a random name to be used if the final distance matrix is empty for any reason
	return(d,firstname)

if __name__ == "__main__":
	# make sure to end paths with '/' here
	for thr in [5,10,20,30,40,50,60,70,80,100,200,500]:
		input_dir = 'trees/newick_trees/'
		output_dir = 'results/results_t' + str(thr) + '/' 
		subprocess.run(['mkdir','-p',output_dir + 'summary/'])
		subprocess.run(['mkdir','-p',output_dir + 'temp/'])
		main(input_dir,output_dir,thresh = thr,rseed = 777)


