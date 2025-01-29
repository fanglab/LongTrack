#!/bin/bash

# ./LongTrack_demo.sh [MAG] [unique_kmer] [metagenome] [conflict_table] [output_dir]
# Usage example: sh LongTrack_demo.sh 

if [ "$#" -ne 0 ]; then
    echo "Usage: sh LongTrack_demo.sh"
    exit 1
fi

MAG=../Data/MAG
unique_kmer=../Data/unique_kmer
metagenome=../Data/metagenome
conflict_table=../Data/conflict_table
output_dir=../Tracking_results

# Load necessary modules
module load python/2.7.16 bowtie2/2.2.8

# set the bowtie path if needed
# PATH=$PATH:[bowtie_path]
# export PATH

date
echo "LongTrack demo Running"

# Create the output directory if it does not exist
mkdir -p "$output_dir"

# Initialize an array to hold all strains
strains=()

# Iterate over each .fna file and extract the strain name
for i in "$MAG"/*.fna
do
  strain=$(basename "$i" .fna)
  strains+=("$strain")
done

# Convert the array to a space-separated string
strains_set=$(printf "%s " "${strains[@]}")

# Count the number of strains
strain_count=${#strains[@]}

# Iterate over each metagenome sample and run the first Python script in parallel
for i in "$metagenome"/*_PE1.fasta
do
  sample=$(basename "$i" _sample_PE1.fasta)
  python find_sample_test_multiple_strains.py "$metagenome/${sample}_sample" \
        "$unique_kmer/" \
        "$output_dir/${sample}_strain_" \
        "$sample" \
        "$strain_count" \
        $strains_set&
done

# Wait for all background processes to finish
wait

# Run the second Python script for each strain
for strain in $strains_set
do
  python understand_sample_test_individual_strain_raw_results.py \
        "$MAG" "$unique_kmer" "$metagenome" "$output_dir" \
        "$strain" \
        "$output_dir/${strain}_results" \
        > "$output_dir/${strain}_results_kmers_reads"
done

# Run the third Python script
python understand_sample_test_allstrains_nextstep_raw_results.py "$MAG" "$metagenome" "$conflict_table" "$output_dir"

# Run the fourth Python script for each strain
for strain in $strains_set
do
  python map_reads_ofstrain_tosamples.py \
        "$output_dir/results_raw_reads" \
        "$strain" \
        "$MAG" "$metagenome" "$conflict_table" "$output_dir"
done

# Run the fifth Python script
python analyze_results_map_reads_ofstrain_tosamples.py "$MAG" "$metagenome" "$conflict_table" "$output_dir"

# Remove all files in the output directory except those matching the specified pattern
find "$output_dir" -type f ! -name 'results_readdistribution_actualreads*' -exec rm -f {} +

python LongTrack_result.py $output_dir

date
echo "Successfully_Finished_All"

