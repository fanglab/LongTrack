#In the previou step, we did  a first pass of results saved on each strain and metagenomic sample. Goal was to find k-mers for a strain that were even observed (to reduce the size, so when we process all strains together our memory requirement remains reasonable), further refine the informative k-mer set by removing those present on a read with <95% informative k-mers, and save such reads after the 1st pass.
#Here, we will process this information from the previous step for all strains together, further remove k-mers and reads that have k-mer ids from multiple strains, and then assign each read (in a metagenomic sample) to a strain uniquely.
import sys
import os
import pdb
import gc
import copy


MAG = sys.argv[1]
metagenome = sys.argv[2]
conflict_table = sys.argv[3]

scratch_path = sys.argv[4]
file_results = scratch_path+"/results_raw"

sample_names = set()
allfiles_sample = os.listdir(metagenome)
for filename in allfiles_sample:
  if "_sample" in filename:
    sample_name = filename.split("_sample")[0]
    sample_names.add(sample_name)
allsample_name = sorted(sample_names)

strain_names = set()
files = os.listdir(MAG)
for filename in files:
  if  ".fna" in filename:
    strain = filename.split('.fna')[0]
    strain_names.add(strain)
allstrain_name = sorted(strain_names)

reads_save = file_results+"_reads"
handle_reads = file(reads_save,"w")


conflict = {}
handle_conflict = file(conflict_table,"r")
content = handle_conflict.readlines()
for entry in content:
  entry = entry.rstrip('\n').split('\t')
  if not entry[0] in conflict:
    conflict[entry[0]] = set()
    for sample in entry[1].split(','):
      conflict[entry[0]].add(sample)

for sample in allsample_name:
  print sample+'\t'+str(sorted(conflict[sample]))

def delete_read(read):		#Delete a read and set all kmer ids on it to 0
  for direction in ['1','2']:
   if direction in dict_read[read]:
    strain_read = dict_read[read]['strain']
    for identifier in dict_read[read][direction]:
      ids_strain[strain_read][identifier] -= 0.1
      if ids_strain[strain_read][identifier] < 0:
        ids_strain[strain_read][identifier] = 0.0
  del dict_read[read]

def remove_ids_current(read_info,strain):	#If a sequencing read is not considered good (or uniquely map to a strain), then we assign all kmer ids on that to be 0, i.e. non informative.
  if ("|") in read_info:
    entry = read_info.rstrip('\n').split("|")
    for temp in entry[1:]:
     ids_strain[strain][temp] -= 0.1
     if ids_strain[strain][temp] < 0:
       ids_strain[strain][temp] = 0.0

def check_ids(read_info,strain):  #If total kmers on a read (after previous step) are less than 95% informative
  good = 0
  total = 0
  if ("|") in read_info:
    entry = read_info.rstrip('\n').split("|")
    total = len(entry)-1
    for temp in entry[1:]:
      good += ids_strain[strain][temp]
    if total == 0:
        return("Unknown",good,total)
    elif good >= 0.95*total:
      return ("Good",good,total)
    else:
      return ("Bad",good,total)	
  else:
    return ("Unknown",good,total)

def process_shared_read(read,strainname,entry1,entry2):		#if a read is shared, i.e. has kmers beloning to multiple strains then ignore the new read (but clean its kmers) if its kmers are <95% informative. In case they are, then likely this read is bad. In that case delete this read and set the kmers for that to be informative. Set the kmers seen for the other strain too as noninformative
  flag1,good1,total1 = check_ids(entry1,strainname)
  flag2,good2,total2 = check_ids(entry2,strainname)
  if (flag1 == "Bad" and flag2 == "Bad") or (flag1 == "Bad" and good2<6) or (flag2 == "Bad" and good1<6):
    delete_read(read)
    bad_read.add(read)
    remove_ids_current(entry1,strainname)
    remove_ids_current(entry2,strainname)
  if (good1>=6 and flag1 == "Good") and (good2>=6 and flag2=="Good"):
    c11 = len(dict_read[read]["1"])
    c12 = len(dict_read[read]["2"])
    c21 = total1
    c22 = total2
    if (c11 >0 and c12>0):
      if (c21>0 and c22>0 and flag1 == "Good" and flag2 == "Good"):
        delete_read(read)
        bad_read.add(read)
        remove_ids_current(entry1,strainname)
        remove_ids_current(entry2,strainname)

def process_new_read(read,strainname,entry1,entry2): #A read is considered if it is >95% informative kmers
  flag1 = "Bad"
  local_entry1 = entry1.rstrip('\n').split("|")
  total1 = len(local_entry1)-1
  good1 = 0
  for temp in local_entry1[1:]:
    good1 += ids_strain[strainname][temp]
  if total1 == 0:
    flag1 = "Unknown"
  elif good1>=0.95*total1 and total1>0:
    flag1 = "Good"
  else:
    flag1 = "Bad"

  flag2 = "Bad"
  local_entry2 = entry2.rstrip('\n').split("|")
  total2 = len(local_entry2)-1
  good2 = 0
  for temp in local_entry2[1:]:
    good2 += ids_strain[strainname][temp]
  if total2 == 0:
    flag2 = "Unknown"
  elif good2>=0.95*total2 and total2>0:
    flag2 = "Good"
  else:
    flag2 = "Bad"

  if (good1>=6 and flag1 == "Good") and (good2>=6 and flag2=="Good"): #6 is 5% of seq read len
    dict_read[read] = {}
    dict_read[read]['strain'] = strainname

    dict_read[read][direction1] = []
    if flag1 == "Good":
      local_entry1 = entry1.rstrip('\n').split("|")
      for temp in local_entry1[1:]:
        dict_read[read][direction1].append(temp)

    dict_read[read][direction2] = []
    if flag2 == "Good":
      local_entry2 = entry2.rstrip('\n').split("|")
      for temp in local_entry2[1:]:
        dict_read[read][direction2].append(temp)

  if (flag1 == "Bad" and flag2 == "Bad") or (flag1 == "Bad" and good2<6) or (flag2 == "Bad" and good1<6):
    bad_read.add(read)
    remove_ids_current(entry1,strainname)
    remove_ids_current(entry2,strainname)

def get_kmer_read_save(read,strainname,entry1,entry2): #Goal of this module is to save the kmers and their counts for a strain and for each metagenomic sample. We will use this to ignore kmers belonging to a strain, that are observed in unrelated samples (based on conflict)
  sample = read.split("_part")[0]
  flag1 = "Bad"
  local_entry1 = entry1.rstrip('\n').split("|")
  total1 = len(local_entry1)-1
  good1 = 0
  for temp in local_entry1[1:]:
    good1 += ids_strain[strainname][temp]
  if total1 == 0:
    flag1 = "Unknown"
  elif good1>=0.95*total1 and total1>0:
    flag1 = "Good"
  else:
    flag1 = "Bad"


  flag2 = "Bad"
  local_entry2 = entry2.rstrip('\n').split("|")
  total2 = len(local_entry2)-1
  good2 = 0
  for temp in local_entry2[1:]:
    good2 += ids_strain[strainname][temp]
  if total2 == 0:
    flag2 = "Unknown"
  elif good2>=0.95*total2 and total2>0:
    flag2 = "Good"
  else:
    flag2 = "Bad"

  if (good1>=6 and flag1 == "Good") and (good2>=6 and flag2=="Good"): #6 is 5% of seq read len
    if flag1 == "Good":
      local_entry1 = entry1.rstrip('\n').split("|")
      for temp in local_entry1[1:]:
	if not temp in dict_counts[strainname][sample]:
	  dict_counts[strainname][sample][temp] = 0
	dict_counts[strainname][sample][temp] += 1
    if flag2 == "Good":
      local_entry2 = entry2.rstrip('\n').split("|")
      for temp in local_entry2[1:]:
        if not temp in dict_counts[strainname][sample]:
          dict_counts[strainname][sample][temp] = 0
        dict_counts[strainname][sample][temp] += 1

dict_read = {}  #Info about reads etc. Read[direction] and id of kmers
bad_read = set()
dict_counts = {} #Info about strain and sample and kmers seen there. For each strain and sample, the ids that are present in that sample and their count
ids_strain = {} #Strain and ids for each strain and if they are 1 or 0
strain_info = {} #Strain final reads and stats. For each strain, the count of reads in a sample
dict_reads_seen = {} #Reads that are seen so can find common ones and process them separately

for strain in allstrain_name:
  strain_info[strain] = {}
  dict_counts[strain] = {}
  for sample in allsample_name:
    strain_info[strain][sample] = 0
    dict_counts[strain][sample] = {}
 
allfiles = os.listdir(scratch_path)		#This is done so we know all about which ids were seen and which reads were bad. Read kmer info from single run of strain so ids are known and bad reads too.
for filename in allfiles:
  if "_results_kmers" in filename:
    strainname = filename.split("_results_kmers")[0]
    if strainname in allstrain_name:
      if not strainname in ids_strain:
        ids_strain[strainname] = {}
      handle = file(scratch_path+"/"+filename,"r")
      content = handle.readline()
      content = content.rstrip('\n').split(',')
      if len(content)>1:
        for temp in content:
	  temp = temp.split()
	  ids_strain[strainname][temp[0]] = float(temp[1])
      print strainname+'\t'+str(len(ids_strain[strainname]))+'\t'+str(sum(ids_strain[strainname].values()))
      del content
      content = handle.readline()
      content = content.rstrip('\n').split(',')
      for temp in content:
	bad_read.add(temp)
      del content
      gc.collect()


######################## Find the reads that are present in other strains too. We will next process such reads to further prune ids
follow_upreads = set()
strains_shared = set()

for filename in allfiles:
  if "_results_reads" in filename:
    strainname = filename.split("_results_reads")[0]
    if strainname in allstrain_name:
      handle = file(scratch_path+"/"+filename,"r")
      entry1 = handle.readline()
      entry2 = handle.readline()
      while(entry1 and entry2):
        if ("|") in entry1:
          temp1 = entry1.split("|")[0]
          read1 = temp1.split(' ')[0]
        else:
          read1 = entry1.rstrip('\n').split(" ")[0]
        if ("|") in entry2:
          temp2 = entry2.split("|")[0]
          read2 = temp2.split(" ")[0]
        else:
          read2 = entry2.rstrip('\n').split(" ")[0]
        read = read1
        if not read in dict_reads_seen:
  	  dict_reads_seen[read] = set()
	dict_reads_seen[read].add(strainname)
	if len(dict_reads_seen[read])>1:
	  follow_upreads.add(read)
          strains_shared = strains_shared.union(dict_reads_seen[read])
        entry1 = handle.readline()
        entry2 = handle.readline()

print "Total follow up reads "+'\t'+str(len(follow_upreads))
follow_upreads = follow_upreads.difference(bad_read)
print "Total follow up reads after removing badreads "+'\t'+str(len(follow_upreads))
print "Total strains which share "+'\t'+str(len(strains_shared))
del dict_reads_seen
del strains_shared
gc.collect()
##################### Process these reads that are shared with other strains and choose only 1 for them. Also find the kmers present in each strain and sample
for filename in allfiles:
  if "_results_reads" in filename:
    strainname = filename.split("_results_reads")[0]
    if strainname in allstrain_name:
      handle = file(scratch_path+"/"+filename,"r")
      entry1 = handle.readline()
      entry2 = handle.readline()
      while(entry1 and entry2):
        if ("|") in entry1:
          temp1 = entry1.split("|")[0]
          read1 = temp1.split(' ')[0]
	  direction1 = temp1.split(' ')[1]
        else:
          read1 = entry1.rstrip('\n').split(" ")[0]
          direction1 = entry1.rstrip('\n').split(" ")[1]
        if ("|") in entry2:
          temp2 = entry2.split("|")[0]
          read2 = temp2.split(" ")[0]
          direction2 = temp2.split(" ")[1]
        else:
          read2 = entry2.rstrip('\n').split(" ")[0]
          direction2 = entry2.rstrip('\n').split(" ")[1]
        read = read1
        if read in follow_upreads:
          if read in dict_read:
	    process_shared_read(read,strainname,entry1,entry2)
          elif not read in dict_read:
	    process_new_read(read,strainname,entry1,entry2)
        get_kmer_read_save(read,strainname,entry1,entry2)	  

        entry1 = handle.readline()
        entry2 = handle.readline()

for strain in allstrain_name:  #This step will find kmers shared between unrelated samples and will set them to 0.
  for sample_name1 in allsample_name:
    share_counts = {}
    for sample_name2 in conflict[sample_name1]: 
      strain_share = set(dict_counts[strain][sample_name1].keys()).intersection(set(dict_counts[strain][sample_name2].keys()))
      for id in strain_share:
	if not id in share_counts:
	  share_counts[id] = 0
	share_counts[id] += dict_counts[strain][sample_name2][id]
    for id in share_counts:
      ids_strain[strain][id] -= 0.1*share_counts[id]
      if ids_strain[strain][id] < 0:
        ids_strain[strain][id] = 0.0
del dict_counts
gc.collect()

################################## Do this now for each strain. Read all the reads and repeat the process of throwing away reads unless we reach steady state

dict_read_copy = copy.deepcopy(dict_read) #Since dict_read at this point only has info for reads that corresponded to multiple strains, but were processed already to have a unique strain. So we copy that read to the dict_read for the strain. The other strains that havenot been processed yet (as only for a single strain) will be processed now too.
del dict_read

for strain in allstrain_name:
  dict_read = {}
  for filename in allfiles:
    if "_results_reads" in filename:
      strainname = filename.split("_results_reads")[0]
      if strainname == strain:
        handle = file(scratch_path+"/"+filename,"r")
        entry1 = handle.readline()
        entry2 = handle.readline()
        while(entry1 and entry2):
          if ("|") in entry1:
            temp1 = entry1.split("|")[0]
            read1 = temp1.split(' ')[0]
            direction1 = temp1.split(' ')[1]
          else:
            read1 = entry1.rstrip('\n').split(" ")[0]
            direction1 = entry1.rstrip('\n').split(" ")[1]
          if ("|") in entry2:
            temp2 = entry2.split("|")[0]
            read2 = temp2.split(" ")[0]
            direction2 = temp2.split(" ")[1]
          else:
            read2 = entry2.rstrip('\n').split(" ")[0]
            direction2 = entry2.rstrip('\n').split(" ")[1]
          read = read1
          if read in bad_read:
            remove_ids_current(entry1,strainname) #Remove ids from current strain
            remove_ids_current(entry2,strainname) #Remove ids from current strain
          elif read in follow_upreads:
	    if read in dict_read_copy:
 	      if dict_read_copy[read]['strain'] == strainname: #For future purposes we would have processed this shared read
	        dict_read[read] = copy.deepcopy(dict_read_copy[read])
	        del dict_read_copy[read]
	  else:
            process_new_read(read,strainname,entry1,entry2)
          entry1 = handle.readline()
          entry2 = handle.readline()

  superflag_repeat = "T"		#Now we redo and iteratively go through each read again, as kmers have been set to non-informative, and more reads will be pruned (and kmers again set to non-informative)
  counter = 0
  bad_half_read = set()
  while( superflag_repeat == "T"):
    if counter >10: #If already ran 10 times and difference is so small then can stop the process
      old_count = count_ids
      count_ids = sum(ids_strain[strain].values())
      if (old_count-count_ids)<=0.01*old_count:
        superflag_repeat = "F"
        break
    counter += 1 
    superflag_repeat = "F"
    print counter
    count_strain = len(dict_read)
    count_ids = sum(ids_strain[strain].values())
    print strain+'\t'+str(count_strain)+'\t'+str(len(ids_strain[strain]))+'\t'+str(sum(ids_strain[strain].values()))

    allreads = dict_read.keys()
    for read in allreads:
      total1 = len(dict_read[read]['1'])
      good1 = 0
      delete_read_flag1 = "F"
      for temp in dict_read[read]['1']:
        good1 += ids_strain[strain][temp] 

      total2 = len(dict_read[read]['2'])
      good2 = 0
      delete_read_flag2 = "F"
      for temp in dict_read[read]['2']:
        good2 += ids_strain[strain][temp]

      if good1<0.95*total1:
        delete_read_flag1 = "T"
      if good2<0.95*total2:
        delete_read_flag2 = "T"

      if (delete_read_flag1 == 'T' and delete_read_flag2 == 'T') or (delete_read_flag1 == 'T' and good2<6) or (delete_read_flag2 == 'T' and good1<6):      
        for direction in ['1','2']:
          if direction in dict_read[read]:
            for identifier in dict_read[read][direction]:
              ids_strain[strain][identifier] -= 0.1
              if ids_strain[strain][identifier] <0:
                ids_strain[strain][identifier] = 0.0
        del dict_read[read]
        superflag_repeat = "T"

      elif total1 == 0:
        dict_read[read]["1"] = []
        bad_half_read.add( (read,"1") )

      elif total2 == 0:
        dict_read[read]["2"] = []
        bad_half_read.add( (read,"2") )
  #Check to see if both half reads present in bad_half_reads
  for (read,direction) in bad_half_read:
    if direction == "1":
      otherdirection = "2"
    else:
      otherdirection = "1"
    if read in dict_read:
      if len(dict_read[read][otherdirection]) == 0:
        del dict_read[read]
      else:
        dict_read[read][direction] = []
  gc.collect()
  for read in dict_read:
    sample = read.split("_part")[0]
    if len(dict_read[read]["1"])>=6 or len(dict_read[read]["2"])>=6:
      strain_info[strain][sample] += 1
      str_write = read+'\t'+sample+'\t'+strain+'\t'+str(len(dict_read[read]["1"]))+'\t'+str(len(dict_read[read]["2"]))
      handle_reads.write(str_write+'\n')
  del dict_read
  gc.collect()

handle_reads.close()
handle = file(file_results,"w")
str_write = "strain"
for sample_name in allsample_name:
  str_write += '\t'+sample_name+"_reads"
handle.write(str_write+'\n')
for strain in allstrain_name:
  str_write = strain
  for sample_name in allsample_name:
    str_write += '\t'+str(strain_info[strain][sample_name])
  handle.write(str_write+'\n')
handle.close()
