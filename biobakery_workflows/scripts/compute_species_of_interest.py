# This script will compute a list of species that can be profiled by StrainPhlAn or PanPhlAn based on min coverage, abundance, and prevalence. 

# To run:
# $ python compute_species_of_interest.py --input input/ --taxonomic-profile input/metaphlan_taxonomic_profiles.tsv --output output/
# There will be a file per sample that includes the species, relative abundance, and coverage.
# A final file will include all of the species that pass the filters plus the number of samples they are found in.

import sys
import os

import gzip
import pickle
import bz2

from anadama2 import Workflow

def read_file_n_lines(file,n):
    line_set=[]

    if file.endswith(".gz"):
        open_function=gzip.open
    else:
        open_function=open

    with open_function(file) as file_handle:
        for line in file_handle:
            if len(line_set) == n:
                yield line_set
                line_set=[]
            line_set.append(line)
    
    # yield the last set
    if len(line_set) == n:
        yield line_set

def get_species_genome_size(pkl_database):
    genome_sizes={}
    db = pickle.load(bz2.BZ2File(pkl_database, 'r'))
    for info in db['taxonomy'].items():
        species="|".join(info[0].split("|")[0:-1])
        length=int(info[1][1])
        genome_sizes[species]=length

    return genome_sizes

def count_nt(infile):
    total_nt=0
    for seq_set in read_file_n_lines(infile,4):
        total_nt+=int(len(seq_set[1].rstrip()))

    return total_nt

def get_taxonomy_by_sample(taxonomic_profile):
    with open(taxonomic_profile) as file_handle:
        samples=[name.replace("_taxonomic_profile","") for name in file_handle.readline().rstrip().split("\t")]
        species={ key: [] for key in samples[1:] }        
        for line in file_handle:
            data=line.split("\t")
            if "s__" in data[0] and not "t__" in data[0]:
                for i, value in enumerate(data[1:]):
                    species[samples[i+1]].append([data[0],float(value)])
    return species

workflow = Workflow()

workflow.add_argument("taxonomic-profile", desc="the merged MetaPhlAn taxonomic profile", required=True)
workflow.add_argument("pkl-database", desc="MetaPhlAn pkl database", default="metaphlan_database/mpa_v30_CHOCOPhlAn_201901.pkl")
workflow.add_argument("min-coverage", desc="Min coverage to select a strain of interest", default=5.0)
workflow.add_argument("min-abundance",desc="the min abundance to select a strain of interest", default="0.0001")
workflow.add_argument("min-prevalence",desc="the min prevalence to select a strain of interest", default="0.1")
workflow.add_argument("input-extension",desc="the extension for the input files", default=".fastq.gz")
workflow.add_argument("unknown-estimate",desc="the estimate of unknown reads (default based on human stool samples from HMP)", default=0.45)

# get the arguments from the command line
args = workflow.parse_args()

def get_strains(task):
    # get all strain genome sizes
    species_genome_sizes=get_species_genome_size(task.depends[2].name)

    # get taxonomy by sample
    taxonomy_by_sample=get_taxonomy_by_sample(task.depends[1].name)

    # compute the coverage for each strain that meets the min abundance plus factor in the unknown read estimate
    total_sample_nt=count_nt(task.depends[0].name)*(1-float(args.unknown_estimate))

    # get the sample id
    sample_id=os.path.basename(task.depends[0].name).replace(args.input_extension,"")

    # open the file of species to write
    write_file_handle=open(task.targets[0].name,"w") 

    # get the total reads and lengths for this sample
    sample_id=os.path.basename(task.depends[0].name).replace(args.input_extension,"")
    for species, value in taxonomy_by_sample.get(sample_id,[]):
        if value > float(args.min_abundance):
            # now check if this meets the coverage requirement
            if species in species_genome_sizes:
                size=species_genome_sizes[species]
                coverage=(total_sample_nt*value)/size

                if coverage > args.min_coverage:
                    write_file_handle.write("\t".join([species,str(value),str(coverage)])+"\n")

    write_file_handle.close()

def compute_final_strains(task):
    # read in all the species from all the samples with min abundance and coverage
    species={}
    for infile in task.depends:
        with open(infile.name) as file_handle:
            for line in file_handle:
                data=line.rstrip().split("\t")
                species[data[0]]=species.get(data[0],0)+1

    # find all strains that pass the prevalence threshold
    total_samples = len(task.depends)
    min_samples = total_samples * float(args.min_prevalence)

    with open(task.targets[0].name,"w") as file_handle:
        for name, count in species.items():
            if count > min_samples:
                file_handle.write("\t".join([name,str(count)])+"\n")

# for each of the input files get a list of strains of interest
all_strain_files=[]
for infile in workflow.get_input_files(extension=args.input_extension):
    outfile = workflow.name_output_files(infile).replace(args.input_extension,"_strains_of_interest.tsv")
    workflow.add_task(
        get_strains,
        depends=[infile,args.taxonomic_profile,args.pkl_database],
        targets=outfile)
    all_strain_files.append(outfile)

# compile the list of strains to incldue those that meet the min prevalence overall
final_file=workflow.name_output_files("all_strains_of_interest.tsv")
workflow.add_task(
    compute_final_strains,
    depends=all_strain_files,
    targets=final_file)

workflow.go()
