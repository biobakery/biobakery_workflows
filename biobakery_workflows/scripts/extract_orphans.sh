#!/bin/bash
# Usage: extract_orphans.sh <RAW INTERLEAVED FASTQ> <BALANCED INTERLEAVED FASTQ> <OUTPUT DIRECTORY>
#
# Given a raw interleaved fastq file and a balanced interleaved fastq (containing no orpahns)
# generate the read ID lists necessary to extract all orphans from the original file.
#

RAW_SEQS=$1
BALANCED_SEQS=$2
OUT_DIR=$3

SAMPLE_NAME=`echo $RAW_SEQS | sed -e 's/\..*//'`

## Grab read IDs from both the raw sequences and the "balanced" sequences
RAW_IDS=${OUT_DIR}/${SAMPLE_NAME}.raw_ids.txt
BALANCED_IDS=${OUT_DIR}/${SAMPLE_NAME}.matched_ids.txt
ORPHAN_IDS=${OUT_DIR}/${SAMPLE_NAME}.orphan_ids.txt

ORPHAN_SEQS=${OUT_DIR}/${SAMPLE_NAME}_orphans.fastq

grep -e '^@.*/[1|2]$' $RAW_SEQS | sed -e 's/^@//' | sort > ${RAW_IDS}
grep -e '^@.*/[1|2]$' $BALANCED_SEQS | sed -e 's/^@//' | sort > ${BALANCED_IDS}

## Produce a list of read IDs in the raw ID list that isn't in the balanced ID list
comm -23 $RAW_IDS $BALANCED_IDS > $ORPHAN_IDS

## And use seqtk to give us an orphans file
seqtk subseq $RAW_SEQS $ORPHAN_IDS > $ORPHAN_SEQS

## Get rid of the ID files we generated
rm $RAW_IDS
rm $BALANCED_IDS
rm $ORPHAN_IDS