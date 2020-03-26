
import os

from anadama2 import Workflow

# to run provide the new workflow run input and output folders
# $ python anadama2_add_files_to_database.py --input $NEW_INPUT_FOLDER --output $NEW_OUTPUT_FOLDER

workflow = Workflow()

# add the list of possible file extensions
workflow.add_argument("input-extensions", desc="the comma-delimited list of extensions of the input files", default="txt,tsv,fastq,fastq.gz,log,sam")
args = workflow.parse_args()

# get all of the files in the input folder with the extensions provided
def get_files_to_add(input_folder):
    posible_extensions=set(args.input_extensions.split(","))
    input_files=[]
    for folder, directories, files in os.walk(input_folder):
        if not ".anadama" in folder:
            for filename in files:
                if any(map(lambda ext: filename.endswith(ext), posible_extensions)):
                    input_files.append(os.path.join(folder,filename))
    return input_files

# get the files to add from the input and output folder
input_files = get_files_to_add(args.input)
output_files = get_files_to_add(args.output)

for filename in input_files+output_files:
    workflow.add_task("echo 'Adding file [depends[0]]'", depends=filename, targets=filename)

workflow.go()

