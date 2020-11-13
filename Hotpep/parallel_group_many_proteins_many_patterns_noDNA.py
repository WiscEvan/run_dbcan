#!/usr/bin/env python

from subprocess import call
import sys
import natsort
import os
import os.path
import re

FILE_DIRNAME = os.path.dirname(os.path.realpath(__file__))

fungus_dirpath = os.path.join(FILE_DIRNAME, "fungus_fungus")
if len(sys.argv) > 2:
	fungus_dirpath = sys.argv[2].replace("?", " ")
peptide_dir_name = os.path.join(FILE_DIRNAME, "CAZY_PPR_patterns", "GH")
if len(sys.argv) > 3:
	peptide_dir_name = sys.argv[3].replace("?", " ")

threads = 8
peptide_length = 6 #length of conserved peptides
hit_cut_off = 3 #number of conserved peptides necessary to classify a protein
freq_cut_off = 1.0 #minimum sum of frequencies necessary to classify a protein

if len(sys.argv) > 1:
	threads = int(sys.argv[1])
if len(sys.argv) > 4:
	peptide_length = int(sys.argv[4])
if len(sys.argv) > 5:
	hit_cut_off = int(sys.argv[5])
if len(sys.argv) > 6:
	freq_cut_off = float(sys.argv[6])

class Protein:
	def __init__(self, sequ):
		self.seq = sequ.upper()
		self.dna = None
		self.name = None
		self.peptides = None
		self.hits = None
		self.freq = None
		self.group = None
		self.subp = None
		self.accession = None
		self.neighbour_seqs = None
		self.ec = None #added by Le Feb 19, 2020
		
	
def callCustom(args):
	return call(args, shell=True)
	
print ("Assigning proteins to groups")
bact_group_many_proteins_many_patterns_src = os.path.join(FILE_DIRNAME, "bact_group_many_proteins_many_patterns.py")
input_args = [fungus_dirpath.replace(" ", "?"), peptide_dir_name.replace(" ", "?"), peptide_length, hit_cut_off, freq_cut_off]
command = "{src} {thread_num} {fungus} {peptide} {length} {hit_cutoff} {freq_cutoff}"
formatted_command = command.format(
	src=bact_group_many_proteins_many_patterns_src,
	thread_num=1,
	fungus=fungus_dirpath.replace(" ", "?"),
	peptide=peptide_dir_name.replace(" ", "?"),
	length=peptide_length,
	hit_cutoff=hit_cut_off,
	freq_cutoff=freq_cut_off,
)
call(formatted_command, shell=True)
print("Collecting Results")

pep_list_array = []
try:
	fams_filepath = os.path.join(peptide_dir_name, "large_fams.txt")
	f = open(fams_filepath, 'r')
except:
	fams_filepath = os.path.join(peptide_dir_name, "fam_list.txt")
	f = open(fams_filepath, 'r')
for line in f:
	pep_list_array.append(line.rstrip())
f.close()
pep_list_hash = {}
for fam in pep_list_array:
	pep_list_hash[fam]=[]
natsort.natsorted(pep_list_array)
thread_number = 1
fam = ""

filepath = os.path.join(fungus_dirpath, "thread{}.txt".format(thread_number))
f = open(filepath, 'r').readlines()
for x in range(len(f)):
	line = f[x].rstrip()
	if line.startswith("Family"):
		fam = line.split(" ")[-1]
	elif line.startswith(">"):
		p = Protein(f[x+1].rstrip())
		p.name = line
		p.peptides = f[x+2].rstrip()
		p.hits = int(f[x+3].rstrip())
		p.freq = float(f[x+4].rstrip())
		p.group = int(f[x+5].rstrip())
		#add by Le start Feb 19, 2020
		fam_main = " ".join(re.findall("[a-zA-Z]+", fam))
		fam_group_ec_fpath = "{peptide}/{family}/{family}_group_ec.txt".format(peptide=peptide_dir_name, family=fam)
		if os.path.exists(fam_group_ec_fpath):
			kk = open(fam_group_ec_fpath).readlines()[p.group-1].rstrip().split("\t")
			if len(kk) > 1:
				p.ec = kk[1]
			else:
				p.ec = "NA"
		else:
			p.ec = "NA"
		#add by Le end Feb 19, 2020
		pep_list_hash[fam].append(p)

output_dir_name = fungus_dirpath+"/Results"
if not os.path.exists(output_dir_name):
	call(["mkdir", output_dir_name])
for fam in pep_list_array:
	hit_array = pep_list_hash[fam]
	#if len(hit_array) > 0:
	fam_file = open(output_dir_name+"/output.txt", "a")
	hit_array.sort(key= lambda x: (x.group, -x.freq, -x.hits))
	for p in hit_array:
		#added ec number by Le start
		fam_file.write(fam+ '\t' +str(p.group)+"\t"+p.name.split(' ')[0][1:]+"\t"+str(p.freq)+"\t"+str(p.hits)+"\t"+p.peptides+"\t"+p.ec+"\n")
		#added ec number by Le end
	fam_file.close()

