
# This script will create a new folder of symlinked files with the new id naming convention using the mapping file.
# $ python rename_files_to_sample_ids.py infolder outfolder mapping file

import os
import sys
import subprocess

infolder=os.path.abspath(sys.argv[1])
outfolder=os.path.abspath(sys.argv[2])

mapping_file=sys.argv[3]

# read in all of the sample ids
ids={}
with open(mapping_file) as file_handle:
    header = file_handle.readline()
    for line in file_handle:
        data = line.rstrip().split("\t")
        ids[data[1]]=data[0]

mapped_ids=set()
unmapped_ids=0
for filename in os.listdir(infolder):
    if not "contam" in filename and not filename.endswith(".log"):
        basename = filename.split(".")[0].split("_")[0]
        if not basename in ids:
            print "Unable to find basename in mapping file: "+basename
            unmapped_ids+=1
        else:
            mapped_ids.add(basename)
            newfilename=ids[basename]+".fastq"
            cmmd=["ln","-s",os.path.join(infolder,filename),os.path.join(outfolder,newfilename)]
            #print cmmd
            subprocess.check_call(cmmd)

print "Total mapped ids: "+ str(len(list(mapped_ids)))
print "Total unmapped ids: "+ str(unmapped_ids)
