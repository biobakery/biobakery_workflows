import cPickle as pickle
import bz2

from anadama2 import Workflow

# To run
# $ python pull_out_reads_by_species_metaphlan_results.py --input input_sam --output output_fastq
# This will look for the metaphlan sam output files named *_bowtie2.sam in the input folder and write
# files *_metaphlan_marker_aligned_subset.fasta to the output folder (one for each input sam file).
# The fasta reads will be any of the sample reads that map to a marker associated with one of the 
# species in the "--species-list" file. This file should have one species per line and be formatted 
# with the metaphlan species naming convention. More specifically, the species file should list
# one per line with metaphlan format (ie "s__Gemella_sanguinis") and for unknown species
# include the genus in this file (ie "s__Gemella_unclassified" should be included in the file as "g__Gemella").
# The metaphlan pkl database is also required for this script to run and can be provided 
# with the option "--pkl-database". 

SAM_READ_NAME_INDEX = 0
SAM_REFERENCE_NAME_INDEX = 2
SAM_SEQ_INDEX = 9

workflow = Workflow()

# input folder should have sam alignment files from metaphlan run
workflow.add_argument("pkl-database", desc="MetaPhlAn2 pkl database", default="metaphlan_db/mpa_v30_CHOCOPhlAn_201901.pkl")
workflow.add_argument("species-list", desc="the list of species to pull reads for", default="species_list.txt")
workflow.add_argument("input-tag-extension", desc="the file name tag and extension", default="_bowtie2.sam")
args = workflow.parse_args()

def find_reads(task):
    # read in the species
    with open(args.species_list) as file_handle:
        species_list = [taxon.rstrip() for taxon in file_handle.readlines()]

    db = pickle.load(bz2.BZ2File(args.pkl_database, 'r'))

    marker_to_species={}
    for marker,info in db['markers'].items():
        if info['clade'] in species_list:
            marker_to_species[marker]=info['clade']

    # read in the sam file and pull out the reads that align with the markers
    with open(task.targets[0].name, "w") as file_handle_write:
        with open(task.depends[0].name) as file_handle:
            for line in file_handle:
                if not line.startswith("@"):
                    data=line.rstrip().split("\t")
                    reference=data[SAM_REFERENCE_NAME_INDEX]
                    if reference in marker_to_species.keys():
                        seq_id = ";".join([data[SAM_READ_NAME_INDEX],marker_to_species[reference]])
                        seq = data[SAM_SEQ_INDEX]
                        file_handle_write.write("\n".join([">"+seq_id,seq])+"\n")

# for each of the input files write the fasta file of reads
for infile in workflow.get_input_files(extension=args.input_tag_extension):
    outfile = workflow.name_output_files(infile).replace(args.input_tag_extension,"_metaphlan_marker_aligned_subset.fasta")
    workflow.add_task(
        find_reads,
        depends=infile,
        targets=outfile)

workflow.go()

