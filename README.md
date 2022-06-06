# Nitfix Project Annotation
Documentation for the pipeline used to annotate genes for long read genomes assembled as part of the NSF Nitfix project. I used two different annotation pipelines: [`MAKER`](http://www.yandell-lab.org/software/maker.html) ([Cantarel et al. 2008 _Genome Research_](https://doi.org/10.1101/gr.6743907)) and [`BRAKER2`](https://github.com/Gaius-Augustus/BRAKER) ([Bruna et al. 2021 _NAR Genomics and Bioinformatics_](https://doi.org/10.1093/nargab/lqaa108)).

## Documentation and Tutorials
I relied heavily on the following resources:
* [Daren Card's annotation pipeline for the boa constrictor genome](https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2)
* [The MAKER tutorial for WGS Assembly and Annotation Winter School 2018](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018)
* [BRAKER2 documentation](https://github.com/Gaius-Augustus/BRAKER#running-braker)

## Genome Assemblies
The following taxa were extracted and sequenced on Nanopore PromethION (long read) and Illumina (short read) platforms by Chris Dervinish and assembled in [`HASLR`](https://github.com/vpc-ccg/haslr) by Neeka Sewnath before being passed to me for annotation.
* _Ceanothus americanus_
* _Coriaria nepalensis_
* _Corynocarpus laevigatus_
* _Euonymus americanus_
* _Frangula americana_
* _Gonopterodendron arborea_ (labeled as _Bulnesia arborea_ in files)
* _Myrica cerifera_
* _Physocarpus capitata_
* _Physocarpus opulifolius_

`MAKER` uses transcript and protein data; its configuration allows for the user to provide transcripts (ESTs) and proteins from the genome taxon and from related taxa. The mode of `BRAKER2` I ran uses only protein data; it allows for proteins from taxa that are distantly related to the genome taxon.

## Reference Sequences
Transcript assemblies from the [1000 Plant Transcriptomes Project](https://db.cngb.org/datamart/plant/DATApla4/) ([One Thousand Plant Transcriptomes Initiative 2019](https://doi.org/10.1038/s41586-019-1693-2), [Carpenter et al. 2019](https://doi.org/10.1093/gigascience/giz126)) were used as input RNA data for `MAKER`. Transcriptomes were chosen based on taxonomic distance to the genome taxon and ploidy level. If a transcriptome for the genome taxon was available, that was chosen. If one wasn't available, the most closely-related taxon with a similar _n_ value based on the [Chromosome Counts Database](http://ccdb.tau.ac.il/) ([Rice et al. 2014](https://doi.org/10.1111/nph.13191)) to the genome taxon was selected. In addition, for a few genome taxa, publically available genomes/transcriptomes were also available and closer taxonomically than the 1KP data. In this case, both the 1KP and other published transcript data were used.

| Genome Taxon               | Transcriptome Taxon               | 1KP Code | Other transcript source                 |
| -------------------------- | --------------------------------- | -------- | --------------------------------------- |
| _Ceanothus americanus_ | _Frangula carolinana_, _Ceanothus thyrsiflorus_ | WVEF | [Salgado et al. 2018](https://doi.org/10.3389/fpls.2018.01629) |
| _Coriaria nepalensis_      | _Coriaria nepalensis_             | NNGU     | NA                                      |
| _Corynocarpus laevigatus_ | _Coriaria nepalensis_, _Datisca glomerata_ | NNGU | [Salgado et al. 2018](https://doi.org/10.3389/fpls.2018.01629) |
| _Euonymus americanus_ | _Crossopetalum rhacoma_, _Tripterygium  wilfordii_ | IHCQ | [Tu et al. 2020](https://doi.org/10.1038/s41467-020-14776-1) |
| _Frangula americana_       | _Frangula caroliniana_            | WVEF      | NA                                    |
| _Gonopterodendron arborea_ | _Tribulus eichleriana_            | KVAY      | NA                                    |
| _Myrica cerifera_          | _Myrica cerifera_                 | INSP      | NA                                    |
| _Physocarpus capitata_     | _Physocarpus opulifolius_         | SXCE      | NA                                    |
| _Physocarpus opulifolius_  | _Physocarpus opulifolius_         | SXCE      | NA                                    |

The same protein data was used as input for all annotation runs. These were the translated CDS from _Medicago trunculata_, _Arachis hypogaea_, and _Glycine max_ genome assemblies, with the longest isoforms selected. The dataset was received from Sara Knaack, who curated it for use in a parallel project. 

## Results
[`BUSCO`](https://busco-archive.ezlab.org/v3/) scores of annotated gene models for both `MAKER` and `BRAKER` are listed in the file `sample_annotation_busco.xlsx`. Genomes and annotations have been uploaded to [`CoGe`](https://genomevolution.org/coge/) for browsing. For now they are private, but I can grant access to individuals if you contact me at kasey.pham@ufl.edu.

## Pipeline Overview

### Repeat Library Construction
This step must be done first, before running `MAKER` or `BRAKER`.

1. `01-assembstats.sh`: Get statistics on genome assembly.
2. `02-repeat01.job`: Model repeats _de novo_ using [`RepeatModeler`](http://www.repeatmasker.org/RepeatModeler/).
3. `03-repeat02.job`: Identify plant repeats from [`Repbase`](https://www.girinst.org/repbase/).
4. `04-repeat03.job`: Mask _de novo_ and plant repeats using [`RepeatMasker`](http://repeatmasker.org/).
5. `05-consreps.sh`: Consolidate _de novo_ and plant repeats.
6. `06-procreps.job`: Process repeats using `RepeatMasker` tools.
7. `07-compreps.job`: Generate GFF files with repeats for annotation tools.

### Run `MAKER`
This step can be run before, after, or concurrent to `BRAKER`.

1. `01-maker01.job`: Run `MAKER` round 1 using annotated genome, transcriptomes, proteomes, and repeat library as input.
2. `02-process_maker.job`: Process and format `MAKER` round 1 output to prepare for running [`SNAP`](https://github.com/KorfLab/SNAP).
3. `03-snap.job`: Get `SNAP` models.
4. `04-train_snap.job`: Generate training sets from `SNAP` models.
5. `05-export_fasta.sh`: Process `MAKER` output for running [`AUGUSTUS`](https://github.com/Gaius-Augustus/Augustus).
6. `06-busco_aug.job`: Run `AUGUSTUS`.
7. `07-process_busco.sh`: Add `AUGUSTUS` output to reference species in config folder.
8. `08-maker02_prep.sh`: Process `AUGUSTUS` output for `MAKER` round 2.
9. `09-maker02.job`: Run `MAKER` round 2.
10. `10-process_maker02.job`: Reformat `MAKER` round 2 output and get annotation statistics.
11. `11-busco.job`: Run `BUSCO` on annotated gene models.

### Run `BRAKER`
This step can be run before, after, or concurrent to `MAKER`.

1. `01-mask_genome.sh`: Mask repeats in genome based on repeat library generated.
2. `02-genemark_es.job`: Run [`GenemarkES`](http://exon.gatech.edu/GeneMark/) to identify gene models.
3. `03-prothint.job`: Run `ProtHint` to generate hints file based on `GeneMark` output to provide to `BRAKER`.
4. `04-braker_prep.sh`: Prepare for running `BRAKER` by creating symlinks to the configuration directory and `BRAKER` scripts in local directory.
5. `05-braker.job`: Run `BRAKER`.
6. `06-get_braker_fas.sh`: Process `BRAKER` output and get `FASTA` file of gene models.
7. `07-busco.job`: Run `BUSCO` on `BRAKER` output gene models.
