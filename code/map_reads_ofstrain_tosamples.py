#This maps the reads assigned to each strain, back to the strain. After considering the high quality mapped reads, we only work with the unique mapped read (for a pileup at a strain) to normalize for different sequencing depth and also to check for relatively uniform pileup across the strain. 
import sys
import pdb
import subprocess
import random
import numpy
import os
import copy

filename_reads = sys.argv[1]
strainname = sys.argv[2]

MAG = sys.argv[3]
metagenome = sys.argv[4]
conflict_table = sys.argv[5]

sample_names = set()
allfiles_sample = os.listdir(metagenome)
for filename in allfiles_sample:
  if "_sample" in filename:
    sample_name = filename.split("_sample")[0]
    sample_names.add(sample_name)

allsample_name = sorted(sample_names)

filename_save_location = sys.argv[6]

conflict = {}
handle_conflict = open(conflict_table,"r")
content = handle_conflict.readlines()
for entry in content:
  entry = entry.rstrip('\n').split('\t')
  if not entry[0] in conflict:
    conflict[entry[0]] = set()
    for sample in entry[1].split(','):
      conflict[entry[0]].add(sample)

handle_results = open(filename_save_location+"/"+strainname+"_overallsummary","w")
reads = {}
reads_final = {}

handle = open(filename_reads,"r") #Save the fasta file for mapping the reads back to the strain
content = handle.readlines()
handle.close()
total_reads_strain = 0
for entry in content:
   entry = entry.split('\t')
   if entry[2]==strainname:
     if not entry[1] in reads:
       reads[entry[1]] = set()
       reads_final[entry[1]] = set()
     reads[entry[1]].add(entry[0])
     reads_final[entry[1]].add(entry[0])
     total_reads_strain += 1
if not os.path.isfile(filename_save_location+"/"+strainname+"_PE1.fasta"):
  handle1_w = open(filename_save_location+"/"+strainname+"_PE1.fasta","w")
  handle2_w = open(filename_save_location+"/"+strainname+"_PE2.fasta","w")
  for samplename in reads:

    handle1 = open(metagenome+"/"+samplename+"_sample_PE1.fasta","r")
    handle2 = open(metagenome+"/"+samplename+"_sample_PE2.fasta","r")
    content1 = handle1.readlines()
    content2 = handle2.readlines()
    total_reads = len(content1)
    index = 0
    while(index<total_reads):
      entry = content1[index].strip('>').split(' ')[0]
      part = entry.split(" ")[0]
      if samplename+"_part_"+part in reads[samplename]:   
        handle1_w.write(">"+samplename+"_part_"+part+" 1\n")
        handle1_w.write(content1[index+1])
        handle2_w.write(">"+samplename+"_part_"+part+" 2\n")
        handle2_w.write(content2[index+1])
        handle1_w.flush()
        handle2_w.flush()
        reads[samplename].remove(samplename+"_part_"+part)	
      index += 2 
  handle1_w.close()
  handle2_w.close()

if not os.path.isfile(filename_save_location+"/"+strainname+".1.bt2"):
  command_run = "bowtie2-build "+MAG+"/"+strainname+".fna "+filename_save_location+"/"+strainname
  subprocess.call(command_run,shell=True)
  print(command_run)
if not os.path.isfile(filename_save_location+"/"+strainname+".sam") or os.path.isfile(filename_save_location+"/"+strainname+".sam") == 0 :
  command_run = "bowtie2 -x "+filename_save_location+"/"+strainname+" -f --very-sensitive --no-mixed --no-discordant -1 "+filename_save_location+"/"+strainname+"_PE1.fasta -2 "+filename_save_location+"/"+strainname+"_PE2.fasta -S "+filename_save_location+"/"+strainname+".sam"
  subprocess.call(command_run,shell=True)
  print(command_run)

samfilename = filename_save_location+"/"+strainname+".sam"
handle = open(samfilename,"r")
dict_read_info = {}
content = handle.readlines()
index = 0
dict_read_info = {}
while(index<len(content)):
  if not content[index][0] == "@":
    entry = content[index].split('\t')
    if entry[2] == "*":
      index += 1
    else:
      pos1 = int(entry[3])+int(entry[8])
      index += 1
      entry = content[index].split('\t')
      pos2 = int(entry[3])+int(entry[8])
      dict_read_info[entry[0]] = [entry[2],min(pos1,pos2),max(pos1,pos2)] #This has the information for a squencing read, its contig for a strain and positions between which it is mapped
  index += 1

def find_unique_reads(dict_read_info,specific_reads): #This returns a dictionary of contigs->positions->read for a sample. We only consider a single read mapping to a position, if there are multiple other present too. We also return the ids of such unique reads for a sample 
  dict_regions = {}
  correct_reads = set()
  for read in dict_read_info:
    if read in specific_reads:
      pos1 = dict_read_info[read][1]
      pos2 = dict_read_info[read][2]
      if not dict_read_info[read][0] in dict_regions:
        dict_regions[dict_read_info[read][0]] = {}
      if not min(pos1,pos2) in dict_regions[dict_read_info[read][0]] and not max(pos1,pos2) in dict_regions[dict_read_info[read][0]]:
        for pos in range(min(pos1,pos2),max(pos1,pos2)):
          dict_regions[dict_read_info[read][0]][pos] = read
        correct_reads.add(read)
      else:
        for pos in range(min(pos1,pos2),max(pos1,pos2)):
          if not pos in dict_regions[dict_read_info[read][0]]:
            dict_regions[dict_read_info[read][0]][pos] = read
            correct_reads.add(read)
      specific_reads.remove(read)
      if len(specific_reads)==0:
        break
  return dict_regions,correct_reads

str1 = ""
str2 = ""
shared_reads = set()
dict_regions_sample = {}
unique_sample = {}

for samplename in allsample_name: #We further finds reads shared between unrelated samples. These reads likely map to the non-unique region of the strain, so we will remove them from all
  if samplename in reads_final:
    dict_regions_sample[samplename],unique_sample[samplename] = find_unique_reads(dict_read_info,reads_final[samplename])
  else:
    dict_regions_sample[samplename] = {}
    unique_sample[samplename] = set()

for samplename1 in allsample_name:
  for samplename2 in conflict[samplename1]:
    for region in dict_regions_sample[samplename1]:
      if region in dict_regions_sample[samplename2]:
        temp_pos = set(dict_regions_sample[samplename1][region].keys()).intersection(set(dict_regions_sample[samplename2][region].keys()))
        for pos in temp_pos:
          shared_reads.add(dict_regions_sample[samplename1][region][pos])
          shared_reads.add(dict_regions_sample[samplename2][region][pos])

for samplename1 in allsample_name: #Save the reads for a strain, and also the reads in the conflict samples
  unique_1 = unique_sample[samplename1]
  unique_c = set()
  for samplename2 in conflict[samplename1]:
    unique_c = unique_c.union(unique_sample[samplename2])
  str1 += str(len(unique_1.difference(shared_reads)))+','+str(len(unique_c.difference(shared_reads)))+ '\t'

handle_results.write(str1.rstrip('\t')+'\n')
handle_results.close()
print("Successfully_Finished\t"+strainname+'\t'+filename_save_location)
