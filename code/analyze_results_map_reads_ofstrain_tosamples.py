#Goal is to find the confidence of occurence of a strain in each metagenomic sample. We work with the mapped reads earlier, and then compare each sample's to its conflict (i.e. unrelated samples). Based on this we assign the confidence scores.
#Simply put if they are 9.5times less occuring then in their conflict samples, then the current sample does not have a strain, and the confidence is 0. For other cases, we find this confidence score accordingly.
import sys
import os
import pdb
import numpy
import copy

MAG = sys.argv[1]
metagenome = sys.argv[2]
conflict_table = sys.argv[3]
output_dir = sys.argv[4]

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
allstrain_name=sorted(strain_names)

conflict = {}
handle_conflict = open(conflict_table,"r")
content = handle_conflict.readlines()
for entry in content:
  entry = entry.rstrip('\n').split('\t')
  if not entry[0] in conflict:
    conflict[entry[0]] = set()
    for sample in entry[1].split(','):
      conflict[entry[0]].add(sample)

to_merge = set()
for sample1 in allsample_name:
  for sample2 in allsample_name:
    if not sample1 == sample2:
      if conflict[sample1].issubset(conflict[sample2]) and conflict[sample2].issubset(conflict[sample1]): #The idea behind this is to merge the negative controls/i.e. conflict samples if they share the same conflicts. The bigger picture is that since all samples do not have the same number of negative controls, so we try to only compare/check against unique negative controls. For ex: Donor in FMT has 3 negative controls. But each neg control has 2 neg controls and alll FMT samples in its conflict.
        to_merge.add( (sample1,sample2) )

handle1 = file(output_dir+"/results_readdistribution_actualreads_confidencescores","w")
handle2 = file(output_dir+"/results_readdistribution_actualreads","w")
str_write = "strain"+'\t'
for sample in allsample_name:
  str_write+=sample+'\t'
handle1.write(str_write.rstrip('\t')+'\n')
handle2.write(str_write.rstrip('\t')+'\n')
actual_reads = {}
for strain in allstrain_name:
  handle = file(output_dir+"/"+strain+"_overallsummary","r")
  content = handle.readlines()

  str_write = ""
  index = 0
  for entry in content[0].rstrip('\n').split('\t'):
    str_write += entry.split(',')[0]+'\t'
    actual_reads[allsample_name[index]] = float(entry.split(',')[0])
    index += 1

  flag_repeat = "T"
  while(flag_repeat=="T"): #Goal here is to identify those samples which have so few reads compared to negative control, that we can set them to 0. dict_del helps in this case. 
    # For example: Donor has v few reads compared to a neg control. But the neg control has very few reads compared to a timepoint 1 year ago. In this case, only the neg control will be set to 0. But the donor will not be.
    flag_repeat = "F"
    dict_del = {}
    for sample in allsample_name:
      for conf_sample in conflict[sample]:
        if 9.5*(actual_reads[sample]+.1)< (actual_reads[conf_sample]+.1) and (actual_reads[conf_sample]-actual_reads[sample])>=5:
          if not sample in dict_del:
            dict_del[sample] = set()
          dict_del[sample].add(conf_sample)
    for sample in allsample_name:
      if sample in dict_del:
        if len(dict_del[sample].difference(set(dict_del.keys())))>0:
          actual_reads[sample] = 0
        else:
          flag_repeat = "T"
  for sample in dict_del.keys():
    actual_reads[sample] = 0

  conflict_new = copy.deepcopy(conflict) #Basically we merge certain conflict samples together. As some samples have more/some less conflicts. This makes the comparison to find the confidence score of occurence more uniform
  for sample1 in allsample_name:
    for sample2 in allsample_name:
      if (sample1,sample2) in to_merge and sample1 in conflict_new and sample2 in conflict_new:
        if actual_reads[sample1]>=actual_reads[sample2]:
          if sample2 in conflict_new:
            del conflict_new[sample2]
        elif sample1 in conflict_new:
          del conflict_new[sample1]

  str_write = ""
  str_write2 = ""
  for sample in allsample_name:
    str_write2 += str(int(actual_reads[sample]))+'\t'
    if actual_reads[sample] >=1:  
      conf_stats = []
      counter = 0
      all_conf = 0
      for conf_sample in conflict[sample]:
        if conf_sample in conflict_new:
          all_conf += 1
          if (actual_reads[sample]+.1) >= 9.5*(actual_reads[conf_sample]+.1) and (actual_reads[sample]-actual_reads[conf_sample])>=1:
            counter += 1
            conf_stats.append((actual_reads[sample]+.1)/(actual_reads[sample]+actual_reads[conf_sample]+.2))
          else:
            conf_stats.append((actual_reads[sample]+.1)/(actual_reads[sample]+actual_reads[conf_sample]+.2))
      if counter== all_conf:
        str_write += '1'+'\t'
      else:
        str_write += str(round(numpy.mean(conf_stats),2))+'\t'
    else:
      str_write += '0'+'\t'
  handle1.write(strain+'\t'+str_write.rstrip('\t')+'\n')
  handle2.write(strain+'\t'+str_write2.rstrip('\t')+'\n')
handle1.close()
handle2.close()
