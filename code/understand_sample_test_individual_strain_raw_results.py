#This is a first pass of results saved on each strain and metagenomic sample. Goal is to find k-mers for a strain that were even observed (to reduce the size, so when we process all strains together our memory requirement remains reasonable), further refine the informative k-mer set by removing those present on a read with <95% informative k-mers, and save such reads after the 1st pass.
#Next we will process this information for all strains together, further remove k-mers and reads that have k-mer ids from multiple strains, and then assign each read (in a metagenomic sample) to a strain uniquely.
import sys
import pdb
import os
import gc
import random

MAG = sys.argv[1]
unique_kmer = sys.argv[2]
metagenome = sys.argv[3]
raw_results_path = sys.argv[4]


sample_names = set()
allfiles_sample = os.listdir(metagenome)
for filename in allfiles_sample:
  if "_sample" in filename:
    sample_name = filename.split("_sample")[0]
    sample_names.add(sample_name)

allsample_name = sorted(sample_names)

strain_process = sys.argv[5]
file_results = sys.argv[6]
handle_read_seen = file(file_results+"_reads","w")

def remove_ids_current(read_info): #If sequencing read has uninformative k-mers then all the kmers from that read are assigned uninformative
  if ("|") in read_info:
    entry = read_info.split("|")
    for temp in entry[1:]:
     temp = temp.split()
     ids_strain[temp[0]] -= 0.1
     if ids_strain[temp[0]] <0:
       ids_strain[temp[0]] = 0.0

def check_ids(read_info):#If 95% of ids seen in a strain (these ids are which are informative after initial scrubbing) are informative (subset of ids seen after initial scrubbing, as now we further refine this set of informative kmers). 
  good = 0
  total = 0
  if ("|") in read_info:
    entry = read_info.split("|")
    for temp in entry[1:]:
      temp = temp.split()
      total += int(temp[1])*ids_strain_deduplicated[temp[0]]
      good += ids_strain[temp[0]]*int(temp[1])
    if total == 0:
	return("Unknown",good,total)
    elif good >= 0.95*total:
      return ("Good",good,total)
    else:
      return ("Bad",good,total)	
  else:
    return ("Unknown",good,total)

ids_strain = {} #Original set of informative kmers after scrubbing. Kmers are further assigned as noninformative.
ids_strain_deduplicated = {} #Original set of informative kmers after scrubbing. This remains static and is not updated
ids_seen_strain = set()
bad_read = set()
total_processed_read = 0
counter_read_added = 0

allfiles_rawresults = os.listdir(raw_results_path)
allfiles_rawresults_workwith = set()
for filename in allfiles_rawresults:
 if "_strain_" in filename:
    samplename_infilename = filename.split("_strain_")[0]
    strainname_infilename = filename.split("_strain_")[1]#.split("_raw_results")[0]
    if samplename_infilename in allsample_name and strainname_infilename == strain_process:
        allfiles_rawresults_workwith.add(filename)
filename = strain_process+"_kmcdb_dump"
handle = file(MAG+"/"+filename,"r") #Open the raw kmer file for each strain
content = handle.readline()
counter = 0
while(content):
  ids_strain[str(counter)] = 0
  ids_strain_deduplicated[str(counter)] = 0
  counter += 1
  content = handle.readline()

filename = strain_process+"_kmcdb_dump_withpos"
handle = file(unique_kmer+"/"+filename,"r") #Open the scrubbed kmer file for each strain
content = handle.readline()
while(content):
  content = content.rstrip('\n').split('\t')
  ids_strain[content[1]] = 1
  ids_strain_deduplicated[content[1]] = 1
  content = handle.readline()
handle.close()
del content
gc.collect()

for filename in allfiles_rawresults_workwith:
  samplename_infilename = filename.split("_strain_")[0]
  strainname_infilename = filename.split("_strain_")[1]#.split("_raw_results")[0]
  handle_raw = file(raw_results_path+"/"+filename,"r")
  counter_filename = 0
  counter = 0
  entry1 = handle_raw.readline()
  counter += 1
  entry2 = handle_raw.readline()
  while(entry1 and entry2):
    if ("|") in entry1:
      temp1 = entry1.split("|")[0]
      read1 = samplename_infilename+"_part_"+temp1.split(" ")[0]
      direction1 = temp1.split(" ")[1]
    else:
      read1 = samplename_infilename+"_part_"+entry1.rstrip('\n').split(" ")[0]
      direction1 = entry1.rstrip('\n').split(" ")[1]

    if ("|") in entry2:
      temp2 = entry2.split("|")[0]
      read2 = samplename_infilename+"_part_"+temp2.split(" ")[0]
      direction2 = temp2.split(" ")[1]
    else:
      read2 = samplename_infilename+"_part_"+entry2.rstrip('\n').split(" ")[0]
      direction2 = entry2.rstrip('\n').split(" ")[1]

    if read1 == read2 and direction1 != direction2:
      read = read1
      direction1 = "1"
      direction2 = "2"
    else:
      print "Something wrong here as reads do not match "+read1+"\t"+read2
    total_processed_read += 1
    flag1,good1,total1 = check_ids(entry1)
    flag2,good2,total2 = check_ids(entry2)

    if (flag1 == "Bad" and flag2 == "Bad") or (flag1 == "Bad" and good2<6) or (flag2 == "Bad" and good1<6):#If reads have very few k-mers, or <95% of kmers are informative, then we assign all kmers on this read as uninformative    
        remove_ids_current(entry1)
        remove_ids_current(entry2)
	bad_read.add(read)	 
    if (good1>=6 and flag1 == "Good") and (good2>=6 and flag2=="Good"): #6 is 5% of seq read len #If read has informative kmers, then we save this read for further processing.
	counter_filename += 1
        counter_read_added += 1 
	str_write = read+' '+"1"
	if flag1 == "Good" and ("|") in entry1:
          local_entry1 = entry1.split("|")
          for temp in local_entry1[1:]:
            temp = temp.split()
	    if ids_strain_deduplicated[temp[0]] == 1:
 	      str_write += "|"+temp[0]
	      ids_seen_strain.add(temp[0])
        handle_read_seen.write(str_write+'\n')

        str_write = read+' '+"2"
        if flag2 == "Good" and ("|") in entry2:
          local_entry2 = entry2.split("|")
          for temp in local_entry2[1:]:
            temp = temp.split()
            if ids_strain_deduplicated[temp[0]] == 1:
              str_write += "|"+temp[0]
              ids_seen_strain.add(temp[0])
        handle_read_seen.write(str_write+'\n')
	
    counter += 1
    entry1 = handle_raw.readline()
    counter += 1
    entry2 = handle_raw.readline()

  print filename+" total_reads until now "+str(counter_read_added)+" total_read_processed "+str(total_processed_read)+" this sample reads "+str(counter_filename)
      
orig_kmers = len(ids_strain)
seen_kmers = len(ids_seen_strain)
seen_goodkmers = 0
for identifier in ids_seen_strain:
   seen_goodkmers += ids_strain[identifier]
handle = file(file_results+"_kmers","w")

str_write = ""
for identifier in ids_seen_strain:
  str_write += identifier+' '+str(ids_strain[identifier])+','
handle.write(str_write.rstrip(',')+'\n')

str_write = ""
for read in bad_read:
   str_write += read+','
handle.write(str_write.rstrip(',')+'\n')

handle.close()
handle_read_seen.close()

print "Finished_processing_strain "+strain_process+" had original kmers "+str(orig_kmers)+' but only these firststep good kmers seen '+str(seen_kmers)+'\t out of which these are good '+str(seen_goodkmers)
