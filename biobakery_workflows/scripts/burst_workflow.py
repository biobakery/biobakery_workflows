#!/usr/bin/env python

import os
import argparse
import subprocess
import time
import shlex

GRID_JOBS_OPTION="--grid-jobs"

STARTUP_TIME = 5 * 60
MAX_SLEEP = 40 * 60

def process_ps_stdout(stdout):
    """ Process the stdout of the ps command """
    return [i.split()[0] for i in filter(lambda x: x, stdout.decode("utf-8").split("\n")[1:])]

def get_children(pid):
    """ Get the pids of the children of this pid """
    try:
        stdout=subprocess.check_output(["ps","--ppid",pid,"-o","pid"])
    except subprocess.CalledProcessError:
        stdout=[]

    pids=[]
    if stdout:
        pids=process_ps_stdout(stdout)

    return pids

def traverse_tree(pid,nodes):
    """ Get all of the children of the pid by walking the tree """

    for child in get_children(pid):
        nodes.update(traverse_tree(child,nodes))
    nodes.add(pid)

    return nodes

def get_pids(pid):
    """ Get the children and their children from parent pid """

    pids=set([pid])
    for child in get_children(pid):
        pids.update(traverse_tree(child,pids))
    
    return list(pids)

def main():

    # gather the arguments
    parser = argparse.ArgumentParser(
        description= "workflow burst wrapper",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--workflow-command", 
        help="the workflow command to run \n[REQUIRED]", 
        metavar="<biobakery_workflows wmgx>", 
        required=True)
    parser.add_argument(
        "--max-sleep", 
        help="the total seconds to wait before terminating the workflow",
        type=int,
        default=MAX_SLEEP)

    args = parser.parse_args()

    # capture the total number of grid jobs
    command=shlex.split(args.workflow_command)
    try:
        grid_jobs = command[command.index(GRID_JOBS_OPTION)+1]
        sleep_time = int(grid_jobs) * 5 + STARTUP_TIME
        if sleep_time > MAX_SLEEP :
            sleep_time = MAX_SLEEP
    except ValueError, IndexError:
        print("Unable to determine grid jobs using default sleep setting")
        sleep_time = MAX_SLEEP

    # set max sleep if decided by the user
    if args.max_sleep != MAX_SLEEP:
        sleep_time = args.max_sleep

    # run the workflow and get the pid
    print("Starting workflow and then will wait for {} seconds".format(sleep_time))
    print(args.workflow_command)
    process=subprocess.Popen(args.workflow_command,shell=True)
    time.sleep(sleep_time)
    print("Stopping workflow. Please run '$ sacct' to track your submitted jobs.")
    # get all pids of children
    pids = get_pids(str(process.pid))
    process.kill()
    for id in pids:
        subprocess.call(["kill","-9",str(id)])

if __name__ == "__main__":
    main()
