#!/usr/bin/env python

########
# NITFIX PROJECT -- LONG READ GENOME ANNOTATION -- 04. CLEAN OUTPUT
# python script to clean intermediate files from MAKER and BRAKER annotation pipeline. Do not run this until you are done with running both pipelines and inspecting outputs!
########

# IMPORT LIBRARIES
import os
import sys

# IMPORT VARIABLES (not secure lmao)
variable_file = sys.argv[1]
exec(open(variable_file).read())

# CLEAN REPEAT LIBRARY FILES
# RepeatModeler
os.system("rm -r {base_dir}/repeat_lib/denovo/RM_*".format(base_dir = BASE_DIR))
os.system("rm -r {base_dir}/repeat_lib/denovo_mask".format(base_dir = BASE_DIR))
# Repbase/RepeatMasker
os.system("rm -r {base_dir}/repeat_lib/repbase_mask".format(base_dir = BASE_DIR))
# Consolidated repeats
os.system("rm {base_dir}/repeat_lib/full_mask/{taxon}_repeat_fullmask.cat".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/repeat_lib/full_mask/{taxon}_repeat_fullmask.fa".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/repeat_lib/full_mask/{taxon}_repeat_fullmask.out".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/repeat_lib/full_mask/{taxon}_repeat_fullmask.out.complex.gff3".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/repeat_lib/full_mask/{taxon}_repeat_fullmask.out.gff3".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/repeat_lib/full_mask/{taxon}_repeat_fullmask.tbl".format(base_dir = BASE_DIR, taxon = TAXON))

# CLEAN MAKER FILES
# MAKER round 1
os.system("rm -r {base_dir}/maker/maker01/{taxon}_maker01.maker.output/{taxon}_maker01_datastore".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm -r {base_dir}/maker/maker01/{taxon}_maker01.maker.output/mpi_blastdb".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker01/{taxon}_maker01.maker.output/{taxon}.all.maker.gff".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker01/{taxon}_maker01.maker.output/{taxon}_maker01.all.maker.proteins.fasta".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker01/{taxon}_maker01.maker.output/{taxon}_maker01.all.maker.transcripts.fasta".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker01/{taxon}_maker01.maker.output/{taxon}_maker01.db".format(base_dir = BASE_DIR, taxon = TAXON))
# SNAP
os.system("rm -r {base_dir}/maker/snap/params".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/{taxon}_maker01.zff.length50_aed0.25.ann".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/snap/{taxon}_maker01.zff.length50_aed0.25.dna".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/snap/alt.ann".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/alt.dna".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/err.ann".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/err.dna".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/export.aa".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/export.ann".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/export.dna".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/export.tx".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/olp.ann".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/olp.dna".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/uni.ann".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/uni.dna".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/wrn.ann".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/snap/wrn.dna".format(base_dir = BASE_DIR))
# AUGUSTUS
os.system("rm -r {base_dir}/maker/augustus/{taxon}_augustus".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm -r {base_dir}/maker/augustus/tmp_opt_BUSCO_{taxon}_augustus".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/augustus/bedtools.err".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/maker/augustus/{taxon}_maker01.all.maker.transcripts1000.fasta".format(base_dir = BASE_DIR, taxon = TAXON))
# MAKER round 2
os.system("rm -r {base_dir}/maker/maker02/{taxon}_maker02.maker.output/{taxon}_maker02_datastore".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm -r {base_dir}/maker/maker02/{taxon}_maker02.maker.output/mpi_blastdb".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}.all.maker.noseq.est2genome.gff".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}.all.maker.noseq.protein2genome.gff".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}.all.maker.noseq.repeats.gff".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}_maker02.maker.output/{taxon}_maker02.all.maker.augustus_masked.proteins.fasta".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}_maker02.maker.output/{taxon}_maker02.all.maker.augustus_masked.transcripts.fasta".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}_maker02.maker.output/{taxon}_maker02.all.maker.gff".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}_maker02.maker.output/{taxon}_maker02.all.maker.non_overlapping_ab_initio.proteins.fasta".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}_maker02.maker.output/{taxon}_maker02.all.maker.non_overlapping_ab_initio.transcripts.fasta".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}_maker02.maker.output/{taxon}_maker02.db".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}_maker02.maker.output/{taxon}_maker02_master_datastore_index.log".format(base_dir = BASE_DIR, taxon = TAXON))
os.system("rm {base_dir}/maker/maker02/{taxon}_maker02.maker.output/seen.dbm".format(base_dir = BASE_DIR, taxon = TAXON))
# BUSCO
os.system("mv {base_dir}/maker/busco/maker_busco/short_summary.specific.embryophyta_odb10.maker_busco.txt {base_dir}/maker/busco".format(base_dir = BASE_DIR))
os.system("rm -r {base_dir}/maker/busco/maker_busco".format(base_dir = BASE_DIR))

# CLEAN BRAKER FILES
# GeneMarkES
os.system("rm -r {base_dir}/braker/genemark_es/data".format(base_dir = BASE_DIR))
os.system("rm -r {base_dir}/braker/genemark_es/info".format(base_dir = BASE_DIR))
os.system("rm -r {base_dir}/braker/genemark_es/output".format(base_dir = BASE_DIR))
os.system("rm -r {base_dir}/braker/genemark_es/run".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/genemark_es/gmhmm.mod".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/genemark_es/run.cfg".format(base_dir = BASE_DIR))
# ProtHint
os.system("rm -r {base_dir}/braker/prothint/diamond".format(base_dir = BASE_DIR))
os.system("rm -r {base_dir}/braker/prothint/Spaln".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/prothint/evidence.gff".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/prothint/gene_stat.yaml".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/prothint/nuc.fasta".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/prothint/prothint.gff".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/prothint/seed_proteins.faa".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/prothint/top_chains.gff".format(base_dir = BASE_DIR))
# BRAKER
os.system("rm -r {base_dir}/braker/braker/braker/GeneMark-EP".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/braker/braker/augustus.hints.aa".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/braker/braker/augustus.hints.codingseq".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/braker/braker/augustus.hints.gtf".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/braker/braker/genemark_evidence.gff".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/braker/braker/genemark_hintsfile.gff".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/braker/braker/genome_header.map".format(base_dir = BASE_DIR))
os.system("rm {base_dir}/braker/braker/braker/hintsfile.gff".format(base_dir = BASE_DIR))
# BUSCO
os.system("mv {base_dir}/braker/busco/braker_busco/short_summary.specific.embryophyta_odb10.braker_busco.txt {base_dir}/braker/busco".format(base_dir = BASE_DIR))
os.system("rm -r {base_dir}/braker/busco/braker_busco".format(base_dir = BASE_DIR))
