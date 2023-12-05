
import os
from anadama2 import Workflow

## Run this workflow to subsample fastq files to "--total-reads <2500000>".
## python use_seqtk_to_downsample.py --input <input_folder> --output <output_folder> --grid-jobs 100

workflow = Workflow()
workflow.add_argument("input-extension", desc="the extensions of the input files", default="fastq.gz")
workflow.add_argument("seqtk-path", desc="the path to the seqtk executable", default="./seqtk/")
workflow.add_argument("total-reads", desc="the number of reads to downsample to", default="2500000")
workflow.add_argument("seed", desc="the seed to use for setk", default="100")
args = workflow.parse_args()


in_files = workflow.get_input_files(extension=args.input_extension)
out_files = workflow.name_output_files(name=in_files)

for infile, outfile in zip(in_files, out_files):
    workflow.add_task_gridable(
        "rm -f [targets[0]] && [args[0]] sample -s [args[1]] [depends[0]] [args[2]] > [args[3]] && gzip [args[3]]",
        depends=infile,
        targets=outfile,
        args=[os.path.join(args.seqtk_path,"seqtk"),args.seed,args.total_reads,outfile.replace(".gz","")],
        cores=1,
        mem=5*1024,
        time=15)

workflow.go()
