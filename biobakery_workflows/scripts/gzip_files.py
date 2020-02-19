# To run: $ python gzip_files.gz --input $FOLDER --output $FOLDER_OUT --grid-jobs 100

from anadama2 import Workflow

workflow = Workflow()

in_files = workflow.get_input_files(extension=".fastq")
out_files =[filename+".gz" for filename in workflow.name_output_files(in_files)]

workflow.add_task_group_gridable(
    "gzip -c [depends[0]] > [targets[0]]",
    depends=in_files,
    targets=out_files,
    docker_image="biobakery/kneaddata:0.7.2_cloud_v4",
    mem=5*1024,
    time=3*60,
    cores=1)

workflow.go()

