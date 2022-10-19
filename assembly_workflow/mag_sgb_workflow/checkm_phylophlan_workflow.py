from anadama2 import Workflow
import glob
import os
import shutil
from anadama2.tracked import TrackedDirectory
from pathlib import Path

workflow = Workflow(version="0.4", description="MAG and SGB workflow")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run any individual command", type=int, default=180000)
workflow.add_argument("checkm-path", desc='path to checkm2 directory with bin/', default="")
workflow.add_argument("checkm-predict-options", desc='checkm2 predict options as a text string with quotes', default="")
workflow.add_argument("completeness", desc='completeness threshold for retaining bins', default=50)
workflow.add_argument("contamination", desc='contamination threshold for retaining bins', default=10)
workflow.add_argument("time", desc="The maximum time in minutes allocated to run any individual command", type=int, default=10000)
workflow.add_argument("n", desc='contamination threshold for retaining bins', default=300, type=int)
workflow.add_argument("phylophlan-database", desc="Database name", default="SGB.Jul20")
workflow.add_argument("phylophlan-database-folder", desc="Folder with the phylophlan database files such as SGB.Jul20/, SGB.Jul20.txt.bz2, etc.")
workflow.add_argument("checkm-data-path", desc='CHECKM_DATA_PATH with selected_marker_sets.tsv, taxon_marker_sets.tsv, etc. if not already set')
workflow.add_argument("phylophlan-metagenomic-options", desc='PhyloPhlAn metagenomic options as a text string with quotes', default="")
args = workflow.parse_args()

# Set CHECKM_DATA_PATH
if args.checkm_data_path:
	os.environ["CHECKM_DATA_PATH"] = args.checkm_data_path

try:
	os.environ["CHECKM_DATA_PATH"]
except:
	raise ValueError("CHECKM_DATA_PATH not provided or set")

database = args.phylophlan_database
database_folder = args.phylophlan_database_folder
cores = args.cores
partition = args.grid_partition

checkm_path = "/" + args.checkm_path.strip("/") + "/"

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"
max_time = args.time

def make_directory(path):
	if not os.path.isdir(path):
		os.makedirs(path)

output = "/" + args.output.strip("/") + "/"
assembly_output = output.split("checkm_phylophlan_steps/")[0]
scratch = "/" + args.grid_scratch.strip("/") + "/"
bins_dir = assembly_output + "bins/"
checkm_scratch = scratch + "checkm/"
phylophlan_scratch = scratch + "phylophlan/"
phylophlan_dir = output + "phylophlan/"
checkm_dir = output + "checkm/"
qa_dir = checkm_dir + "qa/"
qa_unmerged_dir = checkm_dir + "qa_unmerged/"
phylophlan_unmerged_dir = phylophlan_dir + "unmerged/"
checkm_n50_dir = checkm_dir + "n50/"
checkm_bins_dir_scratch = checkm_scratch + "bins/"
checkm_qa_scratch = checkm_scratch + "qa_unmerged/"
phylophlan_unmerged_scratch = phylophlan_scratch + "unmerged/"
make_directory(checkm_qa_scratch)
make_directory(phylophlan_scratch)
make_directory(qa_dir)

if not os.path.isdir(checkm_bins_dir_scratch):
    os.makedirs(checkm_bins_dir_scratch)

    all_paths = Path(bins_dir).rglob('*.bin.*')
    new_files = []
    for path in all_paths:
        file = path.as_posix()
        if not "unbinned" in file and not "tooShort" in file and not "lowDepth" in file and not "final.contigs" in file and not "qc_bins" in file:
            new_files.append(file)

    n = len(new_files)
    size = args.n
    subsets = [new_files[i:i+size] for i in range(0, n, size)]
    for i in range(len(subsets)):
        sub_dir = checkm_bins_dir_scratch + "small_dir_" + str(i) + "/"
        if not os.path.isdir(sub_dir):
            os.makedirs(sub_dir)
        sub_files = subsets[i]
        for file in sub_files:
            shutil.copy(file, sub_dir)

small_dirs = sorted([f.path + "/" for f in os.scandir(checkm_bins_dir_scratch) if f.is_dir()])

# Calculate n50 of MAGs
if not os.path.isfile(checkm_n50_dir + "mags_n50.tsv"):
	n50 = "python " + this_folder + "mag_n50_calc.py -i " + checkm_bins_dir_scratch + " -o " + checkm_n50_dir + " -t " + str(cores)
	workflow.add_task(n50, targets=checkm_n50_dir + "mags_n50.tsv", depends=checkm_bins_dir_scratch)

# Run CheckM
if not os.path.isfile(qa_dir + "quality_report.tsv"):
    for small_dir in small_dirs:
        make_directory(checkm_qa_scratch + small_dir.split("/")[-2] + "/")
        make_directory(qa_unmerged_dir + small_dir.split("/")[-2] + "/")
        command = checkm_path + "bin/checkm2 predict -x fa --input " + small_dir + " --output-directory " + checkm_qa_scratch + small_dir.split("/")[-2] + "/" + " --threads " + str(cores) + " " + args.checkm_predict_options
        workflow.add_task_gridable(actions=command,
            depends=checkm_bins_dir_scratch,
            targets=qa_unmerged_dir + small_dir.split("/")[-2] + "/quality_report.tsv",
            time=min(1500, max_time),
            mem=args.mem,
            cores=args.cores,
            partition=partition)

    for count, small_dir in enumerate(small_dirs):
        if count == 0:
            workflow.add_task("cat " + qa_unmerged_dir + small_dir.split("/")[-2] + "/quality_report.tsv > " + qa_dir + "quality_report_tmp_0.tsv",
        		depends = qa_unmerged_dir + small_dir.split("/")[-2] + "/quality_report.tsv",
        		targets = qa_dir + "quality_report_tmp_0.tsv",
        		)
        else:
            workflow.add_task("cat " + qa_dir + "quality_report_tmp_" + str(count - 1) + ".tsv > " + qa_dir + "quality_report_tmp_" + str(count) + ".tsv && tail -n +2 " + qa_unmerged_dir + small_dir.split("/")[-2] + "/quality_report.tsv >> " + qa_dir + "quality_report_tmp_" + str(count) + ".tsv && rm " + qa_dir + "quality_report_tmp_" + str(count - 1) + ".tsv",
        		depends = [qa_dir + "quality_report_tmp_" + str(count - 1) + ".tsv", qa_unmerged_dir + small_dir.split("/")[-2] + "/quality_report.tsv"],
        		targets = qa_dir + "quality_report_tmp_" + str(count) + ".tsv",
        		)
    workflow.add_task("mv " + qa_dir + "quality_report_tmp_" + str(len(small_dirs) - 1) + ".tsv " + qa_dir + "quality_report.tsv",
        depends = qa_dir + "quality_report_tmp_" + str(len(small_dirs) - 1) + ".tsv",
        targets = qa_dir + "quality_report.tsv",
        )

if not os.path.isfile(qa_dir + "checkm_qa_and_n50.tsv"):
	workflow.add_task("python " + this_folder + "checkm_wrangling.py --checkm-qa " + qa_dir + "quality_report.tsv" + " --n50 " + checkm_n50_dir + "mags_n50.tsv" + " --out_file " + qa_dir + "checkm_qa_and_n50.tsv" + " --completeness " + str(args.completeness) + " --contamination " + str(args.contamination),
		depends = [checkm_n50_dir + "mags_n50.tsv", qa_dir + "quality_report.tsv"],
		targets = qa_dir + "checkm_qa_and_n50.tsv",
		)

if not os.path.isfile(assembly_output + "checkm/qa/" + "checkm_qa_and_n50.tsv"):
	workflow.add_task("mkdir -p " + assembly_output + "checkm/qa/" + " && mv " + qa_dir + "checkm_qa_and_n50.tsv " + assembly_output + "checkm/qa/" + "checkm_qa_and_n50.tsv",
		depends = [qa_dir + "checkm_qa_and_n50.tsv"],
		targets = assembly_output + "checkm/qa/" + "checkm_qa_and_n50.tsv",
		)

##################
# run PhyloPhlAn #
##################

if not os.path.isfile(phylophlan_dir + "phylophlan_out.tsv"):
    for small_dir in small_dirs:
        if not os.path.isfile(phylophlan_unmerged_dir + small_dir.split("/")[-2] + "/" + "phylophlan_out.tsv"):
            make_directory(phylophlan_unmerged_scratch + small_dir.split("/")[-2] + "/")
            make_directory(phylophlan_unmerged_dir + small_dir.split("/")[-2] + "/")
            command = "phylophlan_metagenomic -i " + small_dir + " -n 1 --add_ggb --add_fgb -d " + database + " -o " + phylophlan_unmerged_scratch + small_dir.split("/")[-2] + "/" + "phylophlan_out --nproc " + str(cores) + " --verbose -e fa --database_folder " + database_folder + " " + args.phylophlan_metagenomic_options
            workflow.add_task_gridable(actions=command,
            depends=checkm_bins_dir_scratch,
                targets=phylophlan_unmerged_dir + small_dir.split("/")[-2] + "/" + "phylophlan_out.tsv",
                time=min(2000, max_time),
                mem=args.mem,
                cores=args.cores,
                partition=partition)

    for count, small_dir in enumerate(small_dirs):
        if count == 0:
            workflow.add_task("cat " + phylophlan_unmerged_dir + small_dir.split("/")[-2] + "/" + "phylophlan_out.tsv > " + phylophlan_dir + "phylophlan_out_tmp_0.tsv",
        		depends = phylophlan_unmerged_dir + small_dir.split("/")[-2] + "/" + "phylophlan_out.tsv",
        		targets = phylophlan_dir + "phylophlan_out_tmp_0.tsv",
        		)
        else:
            workflow.add_task("cat " + phylophlan_dir + "phylophlan_out_tmp_" + str(count - 1) + ".tsv > " + phylophlan_dir + "phylophlan_out_tmp_" + str(count) + ".tsv && tail -n +5 " + phylophlan_unmerged_dir + small_dir.split("/")[-2] + "/" + "phylophlan_out.tsv >> " + phylophlan_dir + "phylophlan_out_tmp_" + str(count) + ".tsv && rm " + phylophlan_dir + "phylophlan_out_tmp_" + str(count - 1) + ".tsv",
        		depends = [phylophlan_dir + "phylophlan_out_tmp_" + str(count - 1) + ".tsv", phylophlan_unmerged_dir + small_dir.split("/")[-2] + "/" + "phylophlan_out.tsv"],
        		targets = phylophlan_dir + "phylophlan_out_tmp_" + str(count) + ".tsv",
        		)
    workflow.add_task("mv " + phylophlan_dir + "phylophlan_out_tmp_" + str(len(small_dirs) - 1) + ".tsv " + phylophlan_dir + "phylophlan_out.tsv",
        depends = phylophlan_dir + "phylophlan_out_tmp_" + str(len(small_dirs) - 1) + ".tsv",
        targets = phylophlan_dir + "phylophlan_out.tsv",
        )

if not os.path.isfile(assembly_output + "phylophlan/" + "phylophlan_out.tsv"):
    workflow.add_task("mkdir -p " + assembly_output + "phylophlan/" + " && mv " + phylophlan_dir + "phylophlan_out.tsv " + assembly_output + "phylophlan/" + "phylophlan_out.tsv",
		depends = phylophlan_dir + "phylophlan_out.tsv",
		targets = assembly_output + "phylophlan/" + "phylophlan_out.tsv",
		)

workflow.go()
