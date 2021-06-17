import glob
from collections import defaultdict
import os
import itertools
import sys
#import pandas as pd
from pprint import pprint

minRQccs="0.9"
minRQpostCcs="0.99"
minPasses="1"

RUNIDS,= glob_wildcards( "raw/{runId}.subreads.bam")

shell.prefix("source ~/.bashrc; set +eu; conda deactivate;  set -euo pipefail; ")

pprint(RUNIDS)

rule all:
	input:
		expand("ccs/fastq/{runId}.min-rq" + minRQpostCcs + ".min-passes" + minPasses + ".fastq.gz", runId=RUNIDS)



rule makeCCS:
	input:
		bam="raw/{runId}.subreads.bam",
	output:
		"ccs/{runId}.min-rq" + minRQccs + ".min-passes" + minPasses +".ccs.bam"
	threads: 8
	shell:
		'''
uuid=$(uuidgen)
set +eu
conda activate pacBio.env
set -eu

ccs --minLength=10 --min-rq={minRQccs} --min-passes={minPasses} --max-length 1000000 -j {threads} --report-file ccs/{wildcards.runId}.min-rq{minRQccs}.min-passes{minPasses}.ccs_report.txt {input.bam} {config[TMPDIR]}/$uuid.bam
mv {config[TMPDIR]}/$uuid.bam {output}

		'''

rule rqFilterCcs:
	input: "ccs/{runId}.min-rq" + minRQccs + ".min-passes" + minPasses + ".ccs.bam"
	output: temp("ccs/rqFilter/{runId}.min-rq" + minRQpostCcs + ".min-passes" + minPasses + ".ccs.bam")
	shell:
		'''
uuid=$(uuidgen)
set +eu
conda activate pacBio.env
set -eu

bamtools filter -in {input} -out {config[TMPDIR]}/$uuid -tag "rq":">={minRQpostCcs}"
mv {config[TMPDIR]}/$uuid {output}

		'''

rule removePrimers:
	input: 
		bam="ccs/rqFilter/{runId}.min-rq" + minRQpostCcs + ".min-passes" + minPasses + ".ccs.bam",
		pb_adapters=config["PB_ADAPT"]
	output: 
		temp("ccs/clipped/{runId}.min-rq" + minRQpostCcs + ".min-passes" + minPasses + ".clipped.primer_5p--primer_3p.bam")
	threads: 8
	shell:
		'''
uuid=$(uuidgen)
set +eu
conda activate pacBio.env
set -eu

wd=$PWD
mkdir {config[TMPDIR]}/$uuid
lima --isoseq --min-score-lead 0 -j {threads} {input.bam} {input.pb_adapters} {config[TMPDIR]}/$uuid/$(basename {output[0]} .primer_5p--primer_3p.bam).bam
cd {config[TMPDIR]}/$uuid/
gzip *.report
gzip *.clips
cd $wd
mv {config[TMPDIR]}/$uuid/* $(dirname {output[0]})

		'''

rule removeChimeras:
	input: 
		bam="ccs/clipped/{runId}.min-rq" + minRQpostCcs + ".min-passes" + minPasses + ".clipped.primer_5p--primer_3p.bam",
		pb_adapters=config["PB_ADAPT"]
	output: temp("ccs/nc/{runId}.min-rq" + minRQpostCcs + ".min-passes" + minPasses + ".nc.ccs.bam")
	threads: 8
	shell:
		'''
uuid=$(uuidgen)
set +eu
conda activate pacBio.env
set -eu

wd=$PWD
mkdir {config[TMPDIR]}/$uuid
isoseq3 refine -j {threads}  {input.bam} {input.pb_adapters} {config[TMPDIR]}/$uuid/$(basename {output})

cd {config[TMPDIR]}/$uuid/
gzip *report.csv
cd $wd
mv {config[TMPDIR]}/$uuid/* $(dirname {output})

		'''

rule bamToFastq:
	input: "ccs/nc/{runId}.min-rq" + minRQpostCcs + ".min-passes" + minPasses + ".nc.ccs.bam"
	output: "ccs/fastq/{runId}.min-rq" + minRQpostCcs + ".min-passes" + minPasses + ".fastq.gz"
	shell:
		'''
uuid=$(uuidgen)
set +eu
conda activate pacBio.env
set -eu

bamtools convert -format fastq -in {input} | gzip > {config[TMPDIR]}/$uuid
mv {config[TMPDIR]}/$uuid {output}

		'''

