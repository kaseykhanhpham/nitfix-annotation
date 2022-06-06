#!/usr/bin/env python

########
# NITFIX PROJECT -- LONG READ GENOME ANNOTATION -- 02. MAKER
# This is a script to automate creation of jobs to run the MAKER pipeline on UF's HiperGator cluster
# Change the values of variables in the first section of this script labeled 'VARIABLES'
########

'''
# VARIABLES
## ACCOUNT INFO
ACCOUNT_NAME = "soltis"
EMAIL = "kasey.pham@ufl.edu"
## FILES
TAXON = ""
REF_GENOME = ""
BASE_DIR = ""
TRANSCRIPTOME = ""
PROTEOMES = "/blue/soltis/kasey.pham/nitfix/ref_seqs/ortho/selected_sequences/medicago_truncatula.faa,/blue/soltis/kasey.pham/nitfix/ref_seqs/ortho/selected_sequences/arachis_hypogaea.faa,/blue/soltis/kasey.pham/nitfix/ref_seqs/ortho/selected_sequences/glycine_max.faa"
PIPELINE_DIR = "/blue/soltis/kasey.pham/nitfix/pipeline_scripts"
AUGUSTUS_CONFIG = "/blue/soltis/kasey.pham/nitfix/augustus_config"
BUSCO_DOWNLOADS = "/blue/soltis/kasey.pham/nitfix/busco_downloads"
## PROGRAMS 
INTEL = "intel/2020.0.166"
OPENMPI = "openmpi/4.0.5"
MAKER = "maker/3.01.03"
SNAP = "snap/20100728"
BEDTOOLS = "bedtools/2.30.0"
BUSCO = "busco/4.1.4"
TRNASCAN = "trnascan-se/1.23"
## PROGRAM LOCATIONS

## DERIVED VARIABLES
MAKER_VERSION = MAKER.split("/")[1]
REPEAT_PROTEIN = "/apps/maker/{maker_version}/maker/data/te_proteins.fasta".format(maker_version = MAKER_VERSION)
REPEAT_GFF = "{base_dir}/repeat_lib/full_mask/{taxon}_repeat_fullmask.out.reformat.gff3".format(base_dir = BASE_DIR, taxon = TAXON)
BUSCO_VERSION = BUSCO.split("/")[1]
SNAP_VERSION = SNAP.split("/")[1]
TRNASCAN_VERSION = TRNASCAN.split("/")[1]
'''

# IMPORT LIBRARIES
import os
import sys

# IMPORT VARIABLES (not secure lmao)
variable_file = sys.argv[1]
exec(open(variable_file).read())

# 1. Create directory structure in base directory
os.system("mkdir {base_dir}/maker".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/maker/maker01".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/maker/snap".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/maker/augustus".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/maker/maker02".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/maker/busco".format(base_dir = BASE_DIR))

# 2. MAKER Round 1
# a. Format control files
# opts file
maker01opts_in = open("{pipeline_dir}/maker01_ctl/maker01_opts.ctl".format(pipeline_dir = PIPELINE_DIR), "r")
maker01opts_out = open("{base_dir}/maker/maker01/maker01_opts.ctl".format(base_dir = BASE_DIR), "w")

maker01opts_text = maker01opts_in.readlines()

for line in maker01opts_text:
    line_new = line.replace("REF_GENOME_HERE", REF_GENOME)
    line_new = line_new.replace("TRANSCRIPTOME_HERE", TRANSCRIPTOME)
    line_new = line_new.replace("PROTEOMES_HERE", PROTEOMES)
    line_new = line_new.replace("REPEAT_PROTEIN_HERE", REPEAT_PROTEIN)
    line_new = line_new.replace("REPEAT_GFF_HERE", REPEAT_GFF)
    maker01opts_out.write(line_new)

maker01opts_in.close()
maker01opts_out.close()

# exe file
maker01exe_in = open("{pipeline_dir}/maker01_ctl/maker01_exe.ctl".format(pipeline_dir = PIPELINE_DIR), "r")
maker01exe_out = open("{base_dir}/maker/maker01/maker01_exe.ctl".format(base_dir = BASE_DIR), "w")

maker01exe_text = maker01exe_in.readlines()

for line in maker01exe_text:
    line_new = line.replace("VERSION_HERE", MAKER_VERSION)
    maker01exe_out.write(line_new)

maker01exe_in.close()
maker01exe_out.close()

# evm file
os.system("cp {pipeline_dir}/maker01_ctl/maker01_evm.ctl {base_dir}/maker/maker01/maker01_evm.ctl".format(pipeline_dir = PIPELINE_DIR, base_dir = BASE_DIR, taxon = TAXON))
# bopts file
os.system("cp {pipeline_dir}/maker01_ctl/maker01_bopts.ctl {base_dir}/maker/maker01/maker01_bopts.ctl".format(pipeline_dir = PIPELINE_DIR, base_dir = BASE_DIR, taxon = TAXON))

# b. Run MAKER
maker01job_text='''#!/bin/bash
#SBATCH --job-name=maker01
#SBATCH --output=maker01_%j.out
#SBATCH --error=maker01_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --qos={account_name}
#SBATCH --account={account_name}
#SBATCH --nodes=6                       # Number of nodes
#SBATCH --ntasks=12                     # Number of MPI ranks
#SBATCH --cpus-per-task=1               # Number of cores per MPI rank
#SBATCH --ntasks-per-node=2             # How many tasks on each node
#SBATCH --ntasks-per-socket=1           # How many tasks on each CPU or socket
#SBATCH --distribution=cyclic:cyclic    # Distribute tasks cyclically on nodes and sockets
#SBATCH --mem-per-cpu=4gb               # Memory per processor
#SBATCH -t 14-00:00:00

pwd; hostname; date
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

module purge
module load {intel}
module load {openmpi}
module load {maker}

srun --mpi=pmix maker -base {taxon}_maker01 maker01_opts.ctl maker01_bopts.ctl maker01_exe.ctl
'''.format(email = EMAIL, account_name = ACCOUNT_NAME, intel = INTEL, openmpi = OPENMPI, maker = MAKER, taxon = TAXON)

maker01job_out = open("{base_dir}/maker/maker01/01-maker01.job".format(base_dir = BASE_DIR), "w")
maker01job_out.write(maker01job_text)
maker01job_out.close()

# c. Process MAKER output
processmakerjob_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=procmak
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=2gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=process_maker_%j.out
#SBATCH --error=process_maker_%j.err

module load {maker}

cd {base_dir}/maker/maker01/{taxon}_maker01.maker.output
gff3_merge -s -d {taxon}_maker01_master_datastore_index.log > {taxon}.all.maker.gff
fasta_merge -d {taxon}_maker01_master_datastore_index.log
gff3_merge -n -s -d {taxon}_maker01_master_datastore_index.log > {taxon}.all.maker.noseq.gff
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, maker = MAKER, base_dir = BASE_DIR, taxon = TAXON)

processmakerjob_out = open("{base_dir}/maker/maker01/02-process_maker.job".format(base_dir = BASE_DIR), "w")
processmakerjob_out.write(processmakerjob_text)
processmakerjob_out.close()

# 3. SNAP
# a. Retrieve SNAP-relevant information from MAKER output
snapjob_text='''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=snap
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=1gb
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=snap_%j.out
#SBATCH --error=snap_%j.err

module purge
module load {maker}

maker2zff -x 0.25 -l 50 -d {base_dir}/maker/maker01/{taxon}_maker01.maker.output/{taxon}_maker01_master_datastore_index.log
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, maker = MAKER, base_dir = BASE_DIR, taxon = TAXON)

snapjob_out = open("{base_dir}/maker/snap/03-snap.job".format(base_dir = BASE_DIR), "w")
snapjob_out.write(snapjob_text)
snapjob_out.close()

# b. Get stats on gene models from maker and compile into a training set for SNAP
#     sources:
#     https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2
#     https://reslp.github.io/blog/My-MAKER-Pipeline/
trainsnapjob_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=snap
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=2gb
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=train_snap_%j.out
#SBATCH --error=train_snap_%j.err

module purge
module load {snap}

rename genome {taxon}_maker01.zff.length50_aed0.25 *

fathom {taxon}_maker01.zff.length50_aed0.25.ann {taxon}_maker01.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom {taxon}_maker01.zff.length50_aed0.25.ann {taxon}_maker01.zff.length50_aed0.25.dna -validate > validate.log 2>&1

grep "error" validate.log > error.log
# remove error models with custom python parsing script
python {pipeline_dir}/parse_snap_errors.py {taxon} error.log

# Export gene models which passed error checks plus 1000 basepairs surrounding
fathom {taxon}_maker01.zff.length50_aed0.25.noerr.ann {taxon}_maker01.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

# create training parameters and HMM model
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl {taxon}_maker01.zff.length50_aed0.25 params > {taxon}_maker01.zff.length50_aed0.25.hmm
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, snap = SNAP, taxon = TAXON, pipeline_dir = PIPELINE_DIR)

trainsnapjob_out = open("{base_dir}/maker/snap/04-train_snap.job".format(base_dir = BASE_DIR), "w")
trainsnapjob_out.write(trainsnapjob_text)
trainsnapjob_out.close()

# 4. AUGUSTUS
# a. Export training models from MAKER Round 1 to FASTA sequences
exportfastash_text='''
module load {bedtools}
awk -v OFS="\\t" '{{ if ($3 == "mRNA") print $1, $4, $5 }}' {base_dir}/maker/maker01/{taxon}_maker01.maker.output/{taxon}.all.maker.noseq.gff | \
awk -v OFS="\\t" '{{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }}' | \
bedtools getfasta -fi {ref_genome} -bed - -fo {taxon}_maker01.all.maker.transcripts1000.fasta &> bedtools.err
'''.format(bedtools = BEDTOOLS, base_dir = BASE_DIR, taxon = TAXON, ref_genome = REF_GENOME)

exportfastash_out = open("{base_dir}/maker/augustus/05-export_fasta.sh".format(base_dir = BASE_DIR), "w")
exportfastash_out.write(exportfastash_text)
exportfastash_out.close()

# b. Run BUSCO on MAKER Round 1 models
# Prepare reference files for BUSCO
buscoconfig_in = open("{pipeline_dir}/busco_config/config01.ini".format(pipeline_dir = PIPELINE_DIR), "r")
buscoconfig_out = open("{base_dir}/maker/augustus/config.ini".format(base_dir = BASE_DIR), "w")

buscoconfig_text = buscoconfig_in.readlines()

for line in buscoconfig_text:
    line_new = line.replace("BASE_DIR_HERE", BASE_DIR)
    line_new = line_new.replace("TAXON_HERE", TAXON)
    line_new = line_new.replace("BUSCO_DOWNLOADS_DIR_HERE", BUSCO_DOWNLOADS)
    line_new = line_new.replace("VERSION_HERE", BUSCO_VERSION)
    buscoconfig_out.write(line_new)

buscoconfig_in.close()
buscoconfig_out.close()

# job: busco_aug.job
buscoaugjob_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=busco
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=5gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
#SBATCH --output=busco_aug_%j.out
#SBATCH --error=busco_aug_%j.err

module load {busco}
export BUSCO_CONFIG_FILE="{base_dir}/maker/augustus/config.ini"
export AUGUSTUS_CONFIG_PATH="{augustus_config}"

busco
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, busco = BUSCO, base_dir = BASE_DIR, augustus_config = AUGUSTUS_CONFIG)

buscoaugjob_out = open("{base_dir}/maker/augustus/06-busco_aug.job".format(base_dir = BASE_DIR), "w")
buscoaugjob_out.write(buscoaugjob_text)
buscoaugjob_out.close()

# c. Process BUSCO results
processbuscosh_text= '''#!/bin/bash
cd {base_dir}/maker/augustus/{taxon}_augustus/run_embryophyta_odb10/augustus_output/retraining_parameters/BUSCO_{taxon}_augustus

rename BUSCO_{taxon}_augustus {taxon} *
sed -i 's/BUSCO_{taxon}_augustus/{taxon}/g' {taxon}_parameters.cfg
sed -i 's/BUSCO_{taxon}_augustus/{taxon}/g' {taxon}_parameters.cfg.orig1

mkdir {augustus_config}/species/{taxon}
cp {taxon}* {augustus_config}/species/{taxon}
'''.format(base_dir = BASE_DIR, taxon = TAXON, augustus_config = AUGUSTUS_CONFIG)

processbuscosh_out = open("{base_dir}/maker/augustus/07-process_busco.sh".format(base_dir = BASE_DIR), "w")
processbuscosh_out.write(processbuscosh_text)
processbuscosh_out.close()

# 5. MAKER Round 2
# a. Copy over mapped transcript alignments
maker02prepsh_text='''#!/bin/bash
awk '{{ if ($2 == "est2genome") print $0 }}' {base_dir}/maker/maker01/{taxon}_maker01.maker.output/{taxon}.all.maker.noseq.gff > {base_dir}/maker/maker02/{taxon}.all.maker.noseq.est2genome.gff
# copy over mapped protein alignments
awk '{{ if ($2 == "protein2genome") print $0 }}' {base_dir}/maker/maker01/{taxon}_maker01.maker.output/{taxon}.all.maker.noseq.gff > {base_dir}/maker/maker02/{taxon}.all.maker.noseq.protein2genome.gff
# copy over mapped repeat alignments
awk '{{ if ($2 ~ "repeat") print $0 }}' {base_dir}/maker/maker01/{taxon}_maker01.maker.output/{taxon}.all.maker.noseq.gff > {base_dir}/maker/maker02/{taxon}.all.maker.noseq.repeats.gff
'''.format(base_dir = BASE_DIR, taxon = TAXON)

maker02prepsh_out = open("{base_dir}/maker/maker02/08-maker02_prep.sh".format(base_dir = BASE_DIR), "w")
maker02prepsh_out.write(maker02prepsh_text)
maker02prepsh_out.close()

# b. Edit config files
maker02opts_in = open("{pipeline_dir}/maker02_ctl/maker02_opts.ctl".format(pipeline_dir = PIPELINE_DIR), "r")
maker02opts_out = open("{base_dir}/maker/maker02/maker02_opts.ctl".format(base_dir = BASE_DIR), "w")

maker02opts_text = maker02opts_in.readlines()

for line in maker02opts_text:
    line_new = line.replace("REF_GENOME_HERE", REF_GENOME)
    line_new = line_new.replace("TRANSCRIPTOME_GFF_HERE", "{base_dir}/maker/maker02/{taxon}.all.maker.noseq.est2genome.gff".format(base_dir = BASE_DIR, taxon = TAXON))
    line_new = line_new.replace("PROTEOMES_GFF_HERE", "{base_dir}/maker/maker02/{taxon}.all.maker.noseq.protein2genome.gff".format(base_dir = BASE_DIR, taxon = TAXON))
    line_new = line_new.replace("REPEAT_MAPPED_GFF_HERE", "{base_dir}/maker/maker02/{taxon}.all.maker.noseq.repeats.gff".format(base_dir = BASE_DIR, taxon = TAXON))
    line_new = line_new.replace("SNAP_HMM_HERE", "{base_dir}/maker/snap/{taxon}_maker01.zff.length50_aed0.25.hmm".format(base_dir = BASE_DIR, taxon = TAXON))
    line_new = line_new.replace("TAXON_HERE", TAXON)
    maker02opts_out.write(line_new)

maker02opts_in.close()
maker02opts_out.close()

# exe file
maker02exe_in = open("{pipeline_dir}/maker02_ctl/maker02_exe.ctl".format(pipeline_dir = PIPELINE_DIR), "r")
maker02exe_out = open("{base_dir}/maker/maker02/maker02_exe.ctl".format(base_dir = BASE_DIR), "w")

maker02exe_text = maker02exe_in.readlines()

for line in maker02exe_text:
    line_new = line.replace("MAKER_VERSION_HERE", MAKER_VERSION)
    line_new = line_new.replace("SNAP_VERSION_HERE", SNAP_VERSION)
    line_new = line_new.replace("TRNASCAN_VERSION_HERE", TRNASCAN_VERSION)
    maker02exe_out.write(line_new)

maker02exe_in.close()
maker02exe_out.close()

# evm file
os.system("cp {pipeline_dir}/maker02_ctl/maker02_evm.ctl {base_dir}/maker/maker02/maker02_evm.ctl".format(pipeline_dir = PIPELINE_DIR, base_dir = BASE_DIR, taxon = TAXON))
# bopts file
os.system("cp {pipeline_dir}/maker02_ctl/maker02_bopts.ctl {base_dir}/maker/maker02/maker02_bopts.ctl".format(pipeline_dir = PIPELINE_DIR, base_dir = BASE_DIR, taxon = TAXON))

# c. Run MAKER
maker02job_text = '''#!/bin/bash
#SBATCH --job-name=maker02
#SBATCH --output=maker02_%j.out
#SBATCH --error=maker02_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --qos={account_name}
#SBATCH --account={account_name}
#SBATCH --nodes=6                       # Number of nodes
#SBATCH --ntasks=12                     # Number of MPI ranks
#SBATCH --cpus-per-task=1               # Number of cores per MPI rank
#SBATCH --ntasks-per-node=2             # How many tasks on each node
#SBATCH --ntasks-per-socket=1           # How many tasks on each CPU or socket
#SBATCH --distribution=cyclic:cyclic    # Distribute tasks cyclically on nodes and sockets
#SBATCH --mem-per-cpu=1gb               # Memory per processor
#SBATCH -t 10:00:00

pwd; hostname; date
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

module purge
module load {intel}
module load {openmpi}
module load {maker}
module load {snap}

export AUGUSTUS_CONFIG_PATH={augustus_config}

srun --mpi=pmix maker -base {taxon}_maker02 maker02_opts.ctl maker02_bopts.ctl maker02_exe.ctl
'''.format(email = EMAIL, account_name = ACCOUNT_NAME, intel = INTEL, openmpi = OPENMPI, maker = MAKER, snap = SNAP, augustus_config = AUGUSTUS_CONFIG, taxon = TAXON)

maker02job_out = open("{base_dir}/maker/maker02/09-maker02.job".format(base_dir = BASE_DIR), "w")
maker02job_out.write(maker02job_text)
maker02job_out.close()
#
# b. Process output
processmaker02job_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=procmak
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=3gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=process_maker_%j.out
#SBATCH --error=process_maker_%j.err

module load {maker}

cd {base_dir}/maker/maker02/{taxon}_maker02.maker.output
gff3_merge -s -d  {taxon}_maker02_master_datastore_index.log > {taxon}_maker02.all.maker.gff
fasta_merge -d {taxon}_maker02_master_datastore_index.log
# GFF w/o the sequences
gff3_merge -n -s -d {taxon}_maker02_master_datastore_index.log > {taxon}_maker02.all.maker.noseq.gff

# get stats
# gene num and length
cat {taxon}_maker02.all.maker.gff | awk '{{ if ($3 == "gene") print $0 }}' | awk '{{ sum += ($5 - $4) }} END {{ print NR, sum / NR }}' > gene_counts.txt
perl $HPC_MAKER_DIR/maker/bin/AED_cdf_generator.pl -b 0.025 {taxon}_maker02.all.maker.gff > AED_summary.txt 
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, maker = MAKER, base_dir = BASE_DIR, taxon = TAXON)

processmaker02job_out = open("{base_dir}/maker/maker02/10-process_maker02.job".format(base_dir = BASE_DIR), "w")
processmaker02job_out.write(processmaker02job_text)
processmaker02job_out.close()

# 5. Calculate BUSCO score for final annotations
# a. Format config file
buscoconfig_in = open("{pipeline_dir}/busco_config/config02.ini".format(pipeline_dir = PIPELINE_DIR), "r")
buscoconfig_out = open("{base_dir}/maker/busco/config.ini".format(base_dir = BASE_DIR), "w")

buscoconfig_text = buscoconfig_in.readlines()

for line in buscoconfig_text:
    line_new = line.replace("BASE_DIR_HERE", BASE_DIR)
    line_new = line.replace("BUSCO_DOWNLOADS_DIR_HERE", BUSCO_DOWNLOADS)
    line_new = line_new.replace("VERSION_HERE", BUSCO_VERSION)
    buscoconfig_out.write(line_new)

buscoconfig_in.close()
buscoconfig_out.close()

# b. Run BUSCO
buscojob_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=busco
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=2gb
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
#SBATCH --output=busco_%j.out
#SBATCH --error=busco_%j.err

module load {busco}

IDIR="{base_dir}/maker/maker02/{taxon}_maker02.maker.output"
export BUSCO_CONFIG_FILE="{base_dir}/maker/busco/config.ini"
export AUGUSTUS_CONFIG_PATH="{augustus_config}"

busco -f -i "$IDIR"/{taxon}_maker02.all.maker.transcripts.fasta -o maker_busco
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, busco = BUSCO, base_dir = BASE_DIR, taxon = TAXON, augustus_config = AUGUSTUS_CONFIG)

buscojob_out = open("{base_dir}/maker/busco/11-busco.job".format(base_dir = BASE_DIR), "w")
buscojob_out.write(buscojob_text)
buscojob_out.close()