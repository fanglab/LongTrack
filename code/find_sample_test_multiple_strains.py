#This takes each forward and reverse read of a metagenomics sample, and finds if it has kmers beloning to a strain. If yes, then we save these details
#The forward and reverse reads should match. ELSE ERROR. MAKE SURE AND VERIFY THAT READS ARE TOGETHER IN FILES
import sys
import subprocess
import os
import pdb
import HTSeq
import itertools
import gc
kmer_len = 31
sample_path = sys.argv[1] #Should be path where files are kept
kmer_strain_path_begin = sys.argv[2] #path where kmer of strains are kept. 
fileresults_begin = sys.argv[3]
sample_name = sys.argv[4]
total_strains = sys.argv[5]
strains = set() #Which strains are being processed
for elem in range(6,6+int(total_strains)):
  strains.add(sys.argv[elem])
dict_strain_kmer = {}
for strain in strains:	#Kmers of each strain and the numerical id. Numerical id will be saved as it is smaller in size.
  dict_strain_kmer [strain] = {}
  handle = file(kmer_strain_path_begin+strain+"_kmcdb_dump_withpos","r")
  content = handle.readlines()
  handle.close()
  for entry in content:
    entry = entry.rstrip('\n').split('\t')
    dict_strain_kmer[strain][entry[0]] = str(entry[1])
  del content

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def getkmer_count(seq):	#Kmers present in a read and their counts
  kmer_count = {}
  for x in range(len(seq)+1-kmer_len):
    kmer = seq[x:x+kmer_len]
    kmerr = "".join(complement.get(base, base) for base in reversed(kmer))
    lexico_small = sorted([kmer,kmerr])[0]
    kmer_count[lexico_small] = kmer_count.get(lexico_small, 0) + 1
  return kmer_count
dict_handle = {}
for strain in strains:
  dict_handle[strain] = file(fileresults_begin+strain,"w")

sample_f = sample_path+"_PE1.fasta"
sample_r = sample_path+"_PE2.fasta"
in1 = iter( HTSeq.FastaReader( sample_f ))
in2 = iter( HTSeq.FastaReader( sample_r ))
counter = 0
for read1, read2 in itertools.izip( in1, in2 ):
  counter += 1
  kmercount = {}
  kmercount["1"] = getkmer_count(read1.seq)
  kmercount["2"] = getkmer_count(read2.seq)
  str_write = {}
  for strain in strains:
    str_write = {}
    for direction in ['1','2']:
#GWNJ-0842:904:GW210914000:8:1101:16224:1432 1:N:0:NATCAG
      read_name = read1.name.split(" ")[0]+" "+direction
      str_write[direction] = read_name
      for kmer in kmercount[direction]:
        id = dict_strain_kmer[strain].get(kmer,"NA")
        if not id == "NA":
	  str_write[direction] += '|'+' '+id+' '+str(kmercount[direction][kmer])
    if "|" in str_write["1"] or "|" in str_write["2"]:
      dict_handle[strain].write(str_write["1"]+'\n')
      dict_handle[strain].write(str_write["2"]+'\n')
  if counter%10000 == 0:
    for strain in strains:
   #   print "Finished processing "+str(counter)+" reads for strain "+strain+' sample '+sample_name
      sys.stdout.flush() 
      dict_handle[strain].flush()
for strain in strains:
  print "Successfully_finished_processing_strain_"+strain+"_samplename_"+sample_name
  dict_handle[strain].close()
