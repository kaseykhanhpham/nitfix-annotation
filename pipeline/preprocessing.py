#!/usr/bin/env python

########
# NITFIX PROJECT -- LONG READ GENOME ANNOTATION -- 01. PREPROCESSING
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
## PROGRAMS 
BBMAP = "bbmap/38.44"
REPEATMODELER = "repeatmodeler/2.0"
REPEATMASKER = "repeatmasker/4.0.9"
## DERIVED VARIABLES
REF_GENOME_DIR = "/".join(REF_GENOME.split("/")[0:len(REF_GENOME.split("/")) - 1]) # directory where genome file is saved
REF_GENOME_FILENAME = REF_GENOME.split("/")[len(REF_GENOME.split("/")) - 1] # name of reference genome file
REF_GENOME_BASE = ".".join(REF_GENOME_FILENAME.split(".")[0:len(REF_GENOME_FILENAME.split(".")) - 1]) # name of reference genome without file extension
'''

# IMPORT LIBRARIES
import os
import sys

# IMPORT VARIABLES (not secure lmao)
variable_file = sys.argv[1]
exec(open(variable_file).read())

# 1. Create dir structure for repeats
os.system("mkdir {base_dir}/repeat_lib".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/repeat_lib/denovo".format(base_dir = BASE_DIR))
os.system("mkdir {base_dir}/repeat_lib/full_mask".format(base_dir = BASE_DIR))

# 2. Get assembly statistics
assembstatssh_text ='''#!/bin/bash
module load {bbmap}
$HPC_BBMAP_DIR/bin/statswrapper.sh in={ref_genome} format=2 > {ref_genome_dir}/{basename}_stats.txt
'''.format(bbmap = BBMAP, ref_genome = REF_GENOME, ref_genome_dir = REF_GENOME_DIR, basename = REF_GENOME_BASE)

assembstatssh_out = open("{base_dir}/repeat_lib/01-assembstats.sh".format(base_dir = BASE_DIR), "w")
assembstatssh_out.write(assembstatssh_text)
assembstatssh_out.close()

# 3. RepeatModeler
repeat01job_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=repeat1
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=12gb
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --output=repeat01_%j.out
#SBATCH --error=repeat01_%j.err

pwd; hostname; date
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

module purge
module load {repeatmodeler}

BuildDatabase -name {taxon} -engine ncbi {ref_genome}
RepeatModeler -pa 11 -engine ncbi -database {taxon}
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, repeatmodeler = REPEATMODELER, taxon = TAXON, ref_genome = REF_GENOME)

repeat01job_out = open("{base_dir}/repeat_lib/denovo/02-repeat01.job".format(base_dir = BASE_DIR), "w")
repeat01job_out.write(repeat01job_text)
repeat01job_out.close()

# 4. RepeatMasker with RepBase
repeat02job_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=repeat2
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=20gb
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --output=repeat02_%j.out
#SBATCH --error=repeat02_%j.err

pwd; hostname; date
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

module purge
module load {repeatmasker}

RepeatMasker -pa 12 -e ncbi -species viridiplantae -dir repbase_mask {ref_genome}
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, repeatmasker = REPEATMASKER, ref_genome = REF_GENOME)

repeat02job_out = open("{base_dir}/repeat_lib/03-repeat02.job".format(base_dir = BASE_DIR), "w")
repeat02job_out.write(repeat02job_text)
repeat02job_out.close()

# 5. RepeatMasker with de novo library
repeat03job_text = '''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=repeat3
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=15gb
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --output=repeat03_%j.out
#SBATCH --error=repeat03_%j.err

pwd; hostname; date
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

module purge
module load {repeatmasker}

RepeatMasker -pa 12 -e ncbi -lib {base_dir}/repeat_lib/denovo/{taxon}-families.fa -dir denovo_mask {base_dir}/repeat_lib/repbase_mask/{ref_genome_filename}.masked
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, repeatmasker = REPEATMASKER, base_dir = BASE_DIR, taxon = TAXON, ref_genome_filename = REF_GENOME_FILENAME)

repeat03job_out = open("{base_dir}/repeat_lib/04-repeat03.job".format(base_dir = BASE_DIR), "w")
repeat03job_out.write(repeat03job_text)
repeat03job_out.close()

# 6. Analyze all repeats together
# a. consolidate repeat masking results
consrepssh_text='''#!/bin/bash
ln -s {base_dir}/repeat_lib/denovo_mask/{ref_genome_filename}.masked.masked {base_dir}/repeat_lib/full_mask/{taxon}_repeat_fullmask.fa
ln -s {base_dir}/repeat_lib/denovo_mask/{ref_genome_filename}.masked.out {base_dir}/repeat_lib/full_mask/{taxon}_repeat_fullmask.out
gunzip {base_dir}/repeat_lib/repbase_mask/{ref_genome_filename}.cat.gz
gunzip {base_dir}/repeat_lib/denovo_mask/{ref_genome_filename}.masked.cat.gz
cat {base_dir}/repeat_lib/repbase_mask/{ref_genome_filename}.cat {base_dir}/repeat_lib/denovo_mask/{ref_genome_filename}.masked.cat > {taxon}_repeat_fullmask.cat
'''.format(base_dir = BASE_DIR, ref_genome_filename = REF_GENOME_FILENAME, taxon = TAXON)

consrepssh_out = open("{base_dir}/repeat_lib/full_mask/05-consreps.sh".format(base_dir = BASE_DIR), "w")
consrepssh_out.write(consrepssh_text)
consrepssh_out.close()

# b. Get final masked genome
procrepsjob_text ='''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=procreps
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=5gb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=procreps_%j.out
#SBATCH --error=procreps_%j.err

pwd; hostname; date
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

module purge
module load {repeatmasker}

ProcessRepeats -species viridiplantae {taxon}_repeat_fullmask.cat
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, repeatmasker = REPEATMASKER, taxon = TAXON)

procrepsjob_out = open("{base_dir}/repeat_lib/full_mask/06-procreps.job".format(base_dir = BASE_DIR), "w")
procrepsjob_out.write(procrepsjob_text)
procrepsjob_out.close()

# c. Separate complex repeats
comprepsjob_text ='''#!/bin/bash
#SBATCH --account={account_name}
#SBATCH --qos={account_name}-b
#SBATCH --job-name=compreps
#SBATCH --mail-type=ALL
#SBATCH --mail-user={email}
#SBATCH --mem=50gb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=compreps_%j.out
#SBATCH --error=compreps_%j.err

pwd; hostname; date
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

module purge
module load {repeatmasker}

$HPC_REPEATMASKER_DIR/share/RepeatMasker/util/rmOutToGFF3.pl {taxon}_repeat_fullmask.out > {taxon}_repeat_fullmask.out.gff3
cat {taxon}_repeat_fullmask.out.gff3 | \
  perl -ane '$id; if(!/^\#/){{@F = split(/\\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\\t", @F)."\\n"}} print $_' \
  > {taxon}_repeat_fullmask.out.reformat.gff3
grep -v -e "Satellite" -e ")n" -e "-rich" {taxon}_repeat_fullmask.out.gff3 > {taxon}_repeat_fullmask.out.complex.gff3
cat {taxon}_repeat_fullmask.out.complex.gff3 | \
  perl -ane '$id; if(!/^\#/){{@F = split(/\\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\\t", @F)."\\n"}} print $_' \
  > {taxon}_repeat_fullmask.out.complex_reformat.gff3
'''.format(account_name = ACCOUNT_NAME, email = EMAIL, repeatmasker = REPEATMASKER, taxon = TAXON)

comprepsjob_out = open("{base_dir}/repeat_lib/full_mask/07-compreps.job".format(base_dir = BASE_DIR), "w")
comprepsjob_out.write(comprepsjob_text)
comprepsjob_out.close()

