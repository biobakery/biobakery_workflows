from collections import Counter
from Bio.SeqIO.QualityIO import FastqGeneralIterator #Biopython 1.51 or later
import sys

# Use: python deconcatenate.py input_filename prefix

input_filename = str(sys.argv[1])

output_paired_forward_filename = str(sys.argv[2]) + "_paired_1.fastq"
output_paired_reverse_filename = str(sys.argv[2]) + "_paired_2.fastq"
output_orphan_forward_filename = str(sys.argv[2]) + "_unmatched_1.fastq"
output_orphan_reverse_filename = str(sys.argv[2]) + "_unmatched_2.fastq"

f_suffix = "/1"
r_suffix = "/2"

f_suffix_crop = -len(f_suffix)
def get_name(title):
    name = title.split()[0]
    return name[:f_suffix_crop]

print("Scanning the file to build a list of names...")
ids = [get_name(title) for title, seq, qual in FastqGeneralIterator(open(input_filename))]
paired_lookup = Counter(ids)

print("Processing concatenated file...")
forward_handle = open(output_paired_forward_filename, "w")
orphan_forward_handle = open(output_orphan_forward_filename, "w")
reverse_handle = open(output_paired_reverse_filename, "w")
orphan_reverse_handle = open(output_orphan_reverse_filename, "w")
for title, seq, qual in FastqGeneralIterator(open(input_filename)):
    name = get_name(title)
    if paired_lookup[name] > 1:
        if title[-1] == "1":
            forward_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        else:
            reverse_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    else:
        if title[-1] == "1":
            orphan_forward_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        else:
            orphan_reverse_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

forward_handle.close()
reverse_handle.close()
orphan_forward_handle.close()
orphan_reverse_handle.close()
