
from anadama2 import Workflow

from biobakery_workflows.tasks import dadatwo

workflow = Workflow()
workflow.add_argument("fwd-primer", desc="forward primer, required for its workflow",required=True)
workflow.add_argument("rev-primer", desc="reverse primer, required for its workflow",required=True)
workflow.add_argument("pair-identifier", desc="the string to identify the first file in a pair", default="_R1_001")
workflow.add_argument("threads", desc="number of threads/cores for each task to use", default=1)
args = workflow.parse_args()

dadatwo.remove_primers(workflow,args.fwd_primer,args.rev_primer,args.input,args.output,args.pair_identifier,args.threads)

workflow.go()
