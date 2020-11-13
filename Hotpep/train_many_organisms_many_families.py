#!/usr/bin/env python
#Runs parallel_group_many_proteins_many_patterns.rb for many organisms and types of proteins
#####################################
##Revised by Le Huang on 12/24/2018##
#####################################
from subprocess import call
import os
import os.path
import sys

###
#train_many_organisms_many_families.py [inputFolder] [threads] [hits] [freq]
###

organism_array = ["Chaetomium_globosum_cbs_148_51"]
#Chaetomium_globosum_cbs_148_51 fungus_fungus

cazyme_array = ["CE", "GH", "AA", "PL", "GT", "CBM"]
threads = 1
peptide_length = 6 #length of conserved peptides
hit_cut_off = 3 #number of conserved peptides necessary to classify a protein
freq_cut_off = 1.0 #minimum sum of frequencies necessary to classify a protein
list_multidomain_enzymes = "yes"

if len(sys.argv) > 1:
	organism_array = [sys.argv[1]]
if len(sys.argv) > 2:
	threads = int(sys.argv[2])
if len(sys.argv) > 4:
	hit_cut_off = int(sys.argv[3])
	freq_cut_off = float(sys.argv[4])
#Start Add by Le Huang 12/24/2018
if os.path.exists(organism_array[0]+'/Results/output.txt'):
	call(['rm', organism_array[0]+'/Results/output.txt'])
#End Add by Le Huang 12/24/2018

##Start Delete by Le Huang
#try:
# 	call(['rm', organism_array[0]+'/Results/output.txt'])
# except:
# 	pass
## End Delete by Le Huang
FILE_DIRNAME = os.path.dirname(os.path.realpath(__file__))
cazy_patterns_dir = os.path.join(FILE_DIRNAME, "CAZY_PPR_patterns")
parallel_group_many_proteins_many_patterns_noDNA_src = os.path.join(FILE_DIRNAME, "parallel_group_many_proteins_many_patterns_noDNA.py")

for protein_dir_name in organism_array:
	print("Screening "+protein_dir_name+" for")
	for cazy_class in cazyme_array:
		print(cazy_class)
		peptide_dir_name = os.path.join(cazy_patterns_dir, cazy_class)
		variables =  [threads, protein_dir_name, peptide_dir_name, peptide_length, hit_cut_off, freq_cut_off]
		listed_variables = " ".join([str(x) for x in variables])
		call("{} {}".format(parallel_group_many_proteins_many_patterns_noDNA_src, listed_variables), shell=True)
		#call(["add_functions_orf.py", protein_dir_name, peptide_dir_name])
		var1 = 1
		while var1 <= threads:
			try:
				os.remove(protein_dir_name+"/thread"+str(var1)+".txt")
			except:
				pass
			var1 += 1

	if list_multidomain_enzymes == "yes":
		list_multidomain_proteins_src = os.path.join(FILE_DIRNAME, "list_multidomain_proteins.py")
		call("{} ".format(list_multidomain_proteins_src)+protein_dir_name+" "+"_".join(cazyme_array), shell=True)
print("\nScreened\n"+"\n".join(organism_array))
print("for proteins of the types\n"+", ".join(cazyme_array))
