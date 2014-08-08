#!/bin/bash

# This script will run both versions of isofinder and compare their output.

set -e

#
# Analysis configuration
#
significance=0.95
method=p1
winlen=3000

#
# Input files
#
if [ -z $1 ]; then
    # See if we're on Sean's machine.
    if [ -e /media/student/poly/genome/fa/chr1.fa ]; then
        fastas=$(echo /media/student/poly/genome/fa/chr{{1..22},X,Y}.fa)
    else
        echo "usage: $(basename $0) fasta_file..."
        exit 1
    fi
fi

#
# Perform analysis
#
outdir=/tmp/isofinder-compare
rm -rf $outdir
mkdir -p $outdir

(
    cd version1
    for fasta in $fastas; do
        out=$outdir/$(basename $fasta).v1
        ./isofinder $fasta $significance $method $winlen $out
        bounds=${out}.bounds
        cat $out | awk '{print $2 " " $3 " " $9}' > $bounds
    done
)

(
    cd version2
    for fasta in $fastas; do
        out=$outdir/$(basename $fasta).v2
        ./isofinder $fasta $significance $method $winlen $out
        bounds=${out}.bounds
        cat $out | awk '{print $1 " " $2 " " $3}' > $bounds
    done
)

#
# Compare results
#
diffs=false

for fasta in $fastas; do
    bounds1=$outdir/$(basename $fasta).v1.bounds
    bounds2=$outdir/$(basename $fasta).v2.bounds
    if ! diff -q $bounds1 $bounds2; then
        diffs=true
    fi
done

if ! $diffs; then
    echo "Comparison successful"
fi
