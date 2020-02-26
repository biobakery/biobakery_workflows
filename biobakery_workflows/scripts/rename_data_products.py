# This script will replace the sequence ids with sample ids in all data products
# in the "output_folder", copied from the input_folder using the mappings in the "mapping_file"
# To run: $ python rename_data_products.py infolder outfolder mapping_file.tsv

import sys
import os
import subprocess
import time

input_folder = sys.argv[1]
output_folder = sys.argv[2]
mapping_file = sys.argv[3]

command1 = "sed -i 's/^old/new/g' "
command2 = "sed -i 's/\told/\tnew/g' "

# copy the existing files to the new folder
cmmd=["cp","-r",input_folder,output_folder]
print(cmmd)
subprocess.check_output(cmmd)

# read in the mappings
file_handle = open(mapping_file)
header = file_handle.readline()
mapping={}
for line in file_handle:
    data = line.rstrip().split("\t")
    mapping[data[0]]=data[1]

# rename all sample names in all files in product folder
for file in os.listdir(output_folder):
    file = os.path.join(output_folder,file)
    for seq_id, sample_id in mapping.items():
        for command in [command1, command2]:
            new_command = command.replace("old",seq_id).replace("new",sample_id)+file
            print(new_command)
            subprocess.check_output(new_command,shell=True)
            time.sleep(1)

