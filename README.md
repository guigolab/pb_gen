# Description
`pg_gen` is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow that takes PacBio subread BAM files as input, and generates FASTQ files containing CCS (Circular Consensus Sequencing) reads. It uses standard PacBio software from the [pbbioconda](https://github.com/PacificBiosciences/pbbioconda) / [IsoSeq](https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda) software suite. `pb_gen` runs only on Linux.



# Installation and dependencies

Move to the directory where you intend to run `pb_gen` and clone the source from GitHub:

`git clone https://github.com/julienlag/pb_gen.git`

This workflow uses Snakemake and Conda environments. Therefore, a working installation of [conda](https://docs.anaconda.com/anaconda/install/linux/) and [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) is a prerequisite. The Snakemake/Conda setup ensures that all software `pg_gen` depends on will be installed automatically by Snakemake in an isolated environment the first time you run the workflow.

If using [DRMAA](https://en.wikipedia.org/wiki/DRMAA) (Distributed Resource Management Application API), you should install its Python bindings on your system (in the same environment as the one in which you're running Snakemake):

`conda install -c anaconda drmaa`

# Usage

## Execution

As with any other Snakemake workflow, the command issued will depend on your computing environment. Here's an example which runs well in a UGE/DRMAA cluster environment, in my experience (`{variables}` in curly braces are expanded by Snakemake):

```bash

snakemake -p --reason --latency-wait 100 --use-conda -s master.smk -j 4500 \
--configfile config_pacBio_AlzheimerBrain_ccs.json `# change according to your needs` \
--jobname {rulename}.{jobid}._pb_gen `# job naming rule for the HPC scheduler` \
--cluster-config cluster_config.json `# cluster configuration (queue names, job resource requirements etc.). This file is not provided in the current repo` \
--max-jobs-per-second 10 \
--drmaa " -V -q {cluster.queue} -l disk={cluster.disk} -l virtual_free={cluster.virtual_free} -l h_rt={cluster.h_rt}  -o {cluster.out} -e {cluster.err} {cluster.threads} -P {cluster.project}" \
--rerun-incomplete --keep-going  --show-failed-logs

```


## Input

Mandatory input files include:

- **Subread BAM files** should be placed in a subdirectory named `raw/` and be named according to the following scheme: 

   **`<runId>.subreads.bam`**
   
   where `<runId>` is an arbitrary unique string identifying each sequencing run.

- **FASTA file of adapter sequences**: see `PB_ADAPT` config variable in [Configuration variables](#configuration-variables) below).

## Output 

### CCS FASTQ files

One gzipped FASTQ file is produced per input BAM file. The output FASTQ files are written in the **`ccs/fastq/`** subdirectory and are named according to the following scheme:

`<runId>.min-rq<minRQpostCcs>.min-passes<minPasses>.fastq.gz`

`<runId>` is the basename of the corresponding input subread BAM, while `<minRQpostCcs>` and `<minPasses>` represent the `min-rq` (Minimum predicted accuracy, default **0.99**) and `min-passes` (Minimum number of full-length subreads required to generate a CCS read, default **1**) options passed to the IsoSeq software. These last two variables are currently hardcoded.


### Intermediate files

The output files generated by each intermediate rule are removed after all rules that use it as input are completed (see relevant Snakemake documentation [here](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#protected-and-temporary-files)). Most report and summary files generated along the way by IsoSeq programs are kept, though.

## Configuration variables

These can be specified via the `--config` or `--config-file` Snakemake options (see Snakemake CLI documentation [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html)).

- `PB_ADAPT`: Path to the FASTA file containing the 5' and 3' adapter sequences used in the cDNA library preparation. `PB_ADAPT` is used in rules `removePrimers` (`lima` execution) and `removeChimeras` (`isoseq refine` execution). These primer sequences should be identified as `primer_5p` and `primer_3p` in the FASTA file, respectively. See example file named `pb_adapters.fa` (corresponding to [PacBio's Alzheimer Brain dataset](https://downloads.pacbcloud.com/public/dataset/Alzheimer2019_IsoSeq/)) in the current repository.

- `TMPDIR`: directory to write temporary files to during job executions. If in doubt, set it to `/tmp/`. In some cluster environments, it's useful to set it to `$TMPDIR`, which is dynamically assigned to each job by the scheduler.

# Notes

## Differences with the standard IsoSeq pipeline

The `pb_gen` processing pipeline is pretty close to PacBio's recommendations, with some tweaks. The goal of these modifications is to make the output CCS FASTQs compliant with [LyRic](https://github.com/julienlag/LyRic), our long-read transcriptome analysis workflow.

The most notable differences between `pb_gen` and the [standard IsoSeq pipeline](https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda) are:

- `pb_gen` does not trim poly(A) tails from the reads. 
- `pb_gen` skips the clustering step (`isoseq3 cluster`). Inter-read clustering is considered optional. Clustering reads into non-redundant transcript models is performed by LyRic at a later stage, based on their alignment to the reference genome.

## Online resources

[PacBio glossary of terms](https://www.pacb.com/documentation/pacific-biosciences-glossary-of-terms/)

[PacBio Iso-Seq google group](https://groups.google.com/g/smrt_IsoSeq)
