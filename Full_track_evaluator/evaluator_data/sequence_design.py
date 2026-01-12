#Dec 11th
import pandas as pd

#import pyfaidx
import pysam
import random
import numpy as np
random.seed(10)
# from Bio import SeqIO
# from Bio.Seq import Seq

##Load in .bed files for Enformer and keep test set regions only and extend to 196,608
# input_bed = "/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Full_track_evaluator/evaluator_data/enformer_human_hg38_sequences.bed"
# output_bed = "/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Full_track_evaluator/evaluator_data/enformer_human_hg38_sequences_test196kb.bed"
# EXTEND = 32768
# with open(input_bed) as infile, open(output_bed, "w") as outfile:
#     for line in infile:
#         if line.strip() == "" or line.startswith("#"):
#             continue

#         fields = line.strip().split("\t")
#         if len(fields) < 4:
#             continue

#         chrom, start, end, label = fields[:4]

#         # Keep only test entries
#         if label != "test":
#             continue

#         start = int(start)
#         end = int(end)

#         # Extend each side by 32,768 bp
#         new_start = start - EXTEND
#         new_end = end + EXTEND

#         # Write out updated BED line (preserve extra columns if they exist)
#         new_fields = [chrom, str(new_start), str(new_end), label] + fields[4:]
#         outfile.write("\t".join(new_fields) + "\n")


# ##Load in .bed files for Borzoi and keep test set regions only
# input_bed = "/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Full_track_evaluator/evaluator_data/sequences_human_borzoi.bed"
# output_bed = "//scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Full_track_evaluator/evaluator_data/sequences_human_borzoi_test.bed"

# with open(input_bed) as infile, open(output_bed, "w") as outfile:
#     for line in infile:
#         if line.strip() == "" or line.startswith("#"):
#             continue

#         fields = line.strip().split("\t")

#         # Check if there is a 4th column and it's equal to "fold3"
#         if len(fields) >= 4 and fields[3] == "fold3":
#             outfile.write(line)


fasta_path = "/arc/project/st-cdeboer-1/iluthra/hg38.fa"
genome = pysam.FastaFile(fasta_path)
regions= pd.read_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Full_track_evaluator/evaluator_data/enformer_borzoi_test_merged.bed", sep="\t", header=None, names=["chrom", "start", "end"])

sequences = []
for i in range(0,len(regions)):
    chrom_current = regions['chrom'].iloc[i]
    start_current = regions['start'].iloc[i]
    end_current =regions['end'].iloc[i]
    seq = genome.fetch(chrom_current, start_current,  end_current)
    seq = seq.upper()
    print(len(seq))
    sequences.append(seq)

regions['sequence'] = sequences
regions["region"] = regions["chrom"] + ":" + regions["start"].astype(str) + "-" + regions["end"].astype(str)

print(regions)

regions.to_parquet("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Full_track_evaluator/evaluator_data/enformer_borzoi_test_seqs.parquet", engine='pyarrow', compression='snappy')


# import pyBigWig

# bigwig_file = pyBigWig.open("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/Full_track_evaluator/evaluator_data/ENCFF972GVB.bigWig")

# chrom = "chr8"
# chrom_length = 620000 # GRCh38 length
# #split chromosome into 100kb bins
# segment_size = 150_000

# # Build sequence_coordinates as a flat dict
# sequence_coordinates = {}
# for i, start in enumerate(range(0, chrom_length, segment_size), 1):
#     end = min(start + segment_size, chrom_length)

#     key = f"seq{i}_{chrom}_{start}_{end}"
#     sequence_coordinates[key] = [chrom, start, end]

# #print(sequence_coordinates)



# sequence_measurement_dataframe = pd.DataFrame(columns = ['seq_id','sequence', 'measurements'])

# for seq_id, (chrom, start, end) in sequence_coordinates.items():
#     #print(f"ID: {seq_id}, Chromosome: {chrom}, Start: {start}, End: {end}")
#     seq = genome.fetch(chrom, start,  end)
#     seq = seq.upper()
#     #print(len(seq))
#     #print(seq)
#     measurements = bigwig_file.values(chrom, start, end)
#     #print(len(measurements))
#     data_2_add = pd.DataFrame({'seq_id': seq_id, 'sequence': [seq], 'measurements': [measurements]})
#     sequence_measurement_dataframe = pd.concat([sequence_measurement_dataframe, data_2_add])

# sequence_measurement_dataframe.to_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/Full_track_evaluator/evaluator_data/ENCFF972GVB.tsv", index = False, sep = '\t')
