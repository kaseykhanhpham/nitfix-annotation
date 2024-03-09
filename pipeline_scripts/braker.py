#!/usr/bin/env python

########
# NITFIX PROJECT -- LONG READ GENOME ANNOTATION -- 03. BRAKER
# This is a script to automate creation of jobs to run the BRAKER pipeline on UF's HiperGator cluster
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
PROTEOMES_AGGR = "/blue/soltis/kasey.pham/nitfix/ref_seqs/ortho/annot_test_protein_set.faa"
PIPELINE_DIR = "/blue/soltis/kasey.pham/nitfix/pipeline_scripts"
AUGUSTUS_CONFIG = "/blue/soltis/kasey.pham/nitfix/augustus_config"
BUSCO_DOWNLOADS = "/blue/soltis/kasey.pham/nitfix/busco_downloads"

## PROGRAMS 
BEDTOOLS = "bedtools/2.30.0"
GENEMARK = "genemark_es/4.65"
PYTHON = "python/3.8"
PROTHINT = "prothint/2.6.0"
BRAKER = "braker/2.1.6"
AUGUSTUS = "augustus/3.4.0"
PERL = "perl/5.24.1"
CUFFLINKS = "cufflinks/2.2.1.1"
BUSCO = "busco/4.0.6"

## DERIVED VARIABLES
PROTHINT_VERSION = PROTHINT.split("/")[1]
AUGUSTUS_VERSION = AUGUSTUS.split("/")[1]
GENEMARK_VERSION = GENEMARK.split("/")[1]
BRAKER_VERSION = BRAKER.split("/")[1]
PYTHON_VERSION = PYTHON.split("/")[1]
'''

# IMPORT LIBRARIES
import os
import sys

# IMPORT VARIABLES (not secure lmao)
variable_file = sys.argv[1]
exec(open(variable_file).read())

# 1. Create directory structure for BRAKER files 
os.system("mkdir {base_dir}/braker".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/braker/braker".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/braker/busco".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/braker/genemark_es".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/braker/prothint".format(base_dir = BASE_DIR))

# 2. Soft-mask repeats
maskgenomesh_text = '''#!/bin/bash
module load {bedtools}
cd {base_dir}/repeat_lib/full_mask
bedtools maskfasta -soft -fi {ref_genome} -bed {base_dir}/repeat_lib/full_mask/{taxon}_repeat_fullmask.out.reformat.gff3 -fo {taxon}_genome_fullmask_softmasked.fa
'''.format(bedtools = BEDTOOLS, base_dir = BASE_DIR, ref_genome = REF_GENOME, taxon = TAXON)

maskgenomesh_out = open("{base_dir}/braker/01-mask_genome.sh".format(base_dir = BASE_DIR), "w")
maskgenomesh_out.write(maskgenomesh_text)
maskgenomesh_out.close()

# 2. Run Genemark-ES
genemarkesjob_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=genemark
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=3gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=genemark_%j.out
#SBATCH --error=genemark_%j.err

module purge
module load {genemark}

gmes_petap.pl --sequence {base_dir}/repeat_lib/full_mask/{taxon}_genome_fullmask_softmasked.fa --ES
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, genemark = GENEMARK, base_dir = BASE_DIR, taxon = TAXON)

genemarkesjob_out = open("{base_dir}/braker/genemark_es/02-genemark_es.job".format(base_dir = BASE_DIR), "w")
genemarkesjob_out.write(genemarkesjob_text)
genemarkesjob_out.close()

# 4. Run ProtHint
prothintjob_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=prothint
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=3gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
#SBATCH --output=prothint_%j.out
#SBATCH --error=prothint_%j.err

module purge
module load {python}
module load {prothint}

export PROTHINT_PATH="/apps/prothint/{prothint_version}/prothint/bin"

prothint.py {base_dir}/repeat_lib/full_mask/{taxon}_genome_fullmask_softmasked.fa {proteomes_aggr} --geneMarkGtf {base_dir}/braker/genemark_es/genemark.gtf
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, python = PYTHON, prothint = PROTHINT, prothint_version = PROTHINT_VERSION, base_dir = BASE_DIR, taxon = TAXON, proteomes_aggr = PROTEOMES_AGGR)

prothintjob_out = open("{base_dir}/braker/prothint/03-prothint.job".format(base_dir = BASE_DIR), "w")
prothintjob_out.write(prothintjob_text)
prothintjob_out.close()

# 5. Run BRAKER
# a. Directory and configuration prep
brakerprepsh_text = '''#!/bin/bash
module load {braker}
module load {augustus}
mkdir augustus
ln -s {augustus_config} config
ln -s /apps/augustus/{augustus_version}/bin augustus/bin
mkdir augustus/scripts
cd $HPC_BRAKER_DIR/bin
ls *.pl | while read NAME; do ln -s apps/braker/{braker_version}/bin/"$NAME" {base_dir}/braker/braker/augustus/scripts/"$NAME"; done
'''.format(braker = BRAKER, augustus = AUGUSTUS, augustus_config = AUGUSTUS_CONFIG, augustus_version = AUGUSTUS_VERSION, braker_version = BRAKER_VERSION, base_dir = BASE_DIR)

brakerprepsh_out = open("{base_dir}/braker/braker/04-braker_prep.sh".format(base_dir = BASE_DIR), "w")
brakerprepsh_out.write(brakerprepsh_text)
brakerprepsh_out.close()

# b. Run BRAKER
brakerjob_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=braker
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=6gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=9
#SBATCH --nodes=1
#SBATCH --output=braker_%j.out
#SBATCH --error=braker_%j.err

module purge
module load {braker}
module load {augustus}
module load {perl}

export GENEMARK_PATH="/apps/genemark/genemark-es/{genemark_version}"
export BAMTOOLS_PATH="/apps/braker/{braker_version}/bin/"
export DIAMOND_PATH="/apps/braker/{braker_version}/bin/"
export AUGUSTUS_CONFIG_PATH="{augustus_config}"
export AUGUSTUS_BIN_PATH="/apps/augustus/{augustus_version}/bin"
export AUGUSTUS_SCRIPTS_PATH="/apps/braker/{braker_version}/bin"
export PYTHON3_PATH="/apps/python/{python_version}/bin"
export CDBTOOLS_PATH="/apps/braker/{braker_version}/bin/"

braker.pl --epmode --species {taxon}_noRNA --genome={base_dir}/repeat_lib/full_mask/{taxon}_genome_fullmask_softmasked.fa --hints={base_dir}/braker/prothint/prothint_augustus.gff --softmasking --cores=8
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, braker = BRAKER, augustus = AUGUSTUS, perl = PERL, genemark_version = GENEMARK_VERSION, braker_version = BRAKER_VERSION, augustus_config = AUGUSTUS_CONFIG, augustus_version = AUGUSTUS_VERSION, python_version = PYTHON_VERSION, taxon = TAXON, base_dir = BASE_DIR)

brakerjob_out = open("{base_dir}/braker/braker/05-braker.job".format(base_dir = BASE_DIR), "w")
brakerjob_out.write(brakerjob_text)
brakerjob_out.close()

# 7. Calculate BUSCO score for annotations
# a. Get FASTA of BRAKER annotations and configure BUSCO
# location: /blue/soltis/kasey.pham/nitfix/bulnesia/braker/braker/braker
getbrakerfassh_text = '''
module load {cufflinks}

cd {base_dir}/braker/braker/braker
gffread -x braker_models.fasta -g {base_dir}/repeat_lib/full_mask/{taxon}_genome_fullmask_softmasked.fa braker.gtf
'''.format(cufflinks = CUFFLINKS, base_dir = BASE_DIR, taxon = TAXON)

getbrakerfassh_out = open("{base_dir}/braker/braker/06-get_braker_fas.sh".format(base_dir = BASE_DIR), "w")
getbrakerfassh_out.write(getbrakerfassh_text)
getbrakerfassh_out.close()

# b. Format BUSCO config file
buscoconfig_in = open("{pipeline_dir}/busco_config/config02.ini".format(pipeline_dir = PIPELINE_DIR), "r")
buscoconfig_out = open("{base_dir}/braker/busco/config.ini".format(base_dir = BASE_DIR), "w")

buscoconfig_text = buscoconfig_in.readlines()

for line in buscoconfig_text:
    line_new = line.replace("BASE_DIR_HERE", BASE_DIR)
    line_new = line.replace("BUSCO_DOWNLOADS_DIR_HERE", BUSCO_DOWNLOADS)
    line_new = line_new.replace("VERSION_HERE", BUSCO_VERSION)
    buscoconfig_out.write(line_new)

buscoconfig_in.close()
buscoconfig_out.close()

# c. Run BUSCO
buscojob_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=busco
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=2gb
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --output=busco_%j.out
#SBATCH --error=busco_%j.err

module load {busco}
IDIR="{base_dir}/braker/braker/braker"
export BUSCO_CONFIG_FILE="{base_dir}/braker/busco/config.ini"
export AUGUSTUS_CONFIG_PATH="{augustus_config}"

busco -f -i "$IDIR"/braker_models.fasta -o braker_busco
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, busco = BUSCO, base_dir = BASE_DIR, augustus_config = AUGUSTUS_CONFIG)


buscojob_out = open("{base_dir}/braker/busco/07-busco.job".format(base_dir = BASE_DIR), "w")
buscojob_out.write(buscojob_text)
buscojob_out.close()