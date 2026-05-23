# Code and intermediate data for the paper "Genomic architecture of the developing forelegs of Drosophila prolongata; an exaggerated weapon and ornament" 

If you use this data and or code, please cite the forthcoming paper.

Tyler Audet, Jhoniel Perdigón Ferreira, Abhishek Meena, Mariam Abass, Arteen Torabi-Marashi, John Yeom, Fatima Zaghloul, Nour Zaghloul, Arie Mizrahi, Julia Novikov, Emma Xi, Stefan Lüpold, Ian Dworkin. (2026).Genomic architecture of the developing forelegs of *Drosophila prolongata*; an exaggerated weapon and ornament. [Evolution](https://doi.org/10.1093/evolut/qpag096). 

Raw sequence data is available via [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) via the following accession GSE330706.

Static version of this repository, plus additional intermediate data too large for github are available here
10.6084/m9.figshare.32311218

## Overview

Beginning with the RNAseq data in a folder called 'RNA', and the reference genomes for carrolli, prolongata, and melanogaster in a folder called 'genome' and sub-directories named after each species. Bash scripts contain path structures for all scripts and are formatted for use with slurm schedulers.

Scripts:

-R_scripts

  -DEG_analysis.Rmd: This is the script for all differential gene expression analyses for the paper 'Genomic architecture of the developing forelegs of Drosophila prolongata; an exaggerated weapon and ornament'. This. script reads in the modelled emmeans estimates from `Model.R` which it uses to generate shrunken estimates and outputs these estimates as `[species]_norms` which are used in this and other scripts as our final differential gene expression estimates. This script generates Figures 4 and 5, as well as supplemental figures S2, S3, S7, S8, S9.

  -PCA_work.Rmd: This is the script for the PCA analyses for the paper 'Genomic architecture of the developing forelegs of Drosophila prolongata; an exaggerated weapon and ornament'. This script reads in raw count data generated from `STAR_[species].sh` and gemerates Figures 2 and 3, as well as supplemental figure S6.

  -PenGAL4_cross.Rmd: This is the script for the function pendulin GAL4 cross analyses for the paper 'Genomic architecture of the developing forelegs of Drosophila prolongata; an exaggerated weapon and ornament'. This script reads in our phenotypic measurements `pengal4_data.csv` and outputs figure 7 and supplemental figure S10.

  -DeSeq_analysis.Rmd: This is the script used to verify our methods against standard DeSeq methods for the paper 'Genomic architecture of the developing forelegs of Drosophila prolongata; an exaggerated weapon and ornament'. This script reads in raw counts generated from the script `STAR_[species].sh` to run the DeSeq analyses as well as our shrunked estimates as the files `[species]_norms` generated from the script `DGE_analysis.Rmd`. It generates the supplemental figures S4 and S5.

  -Vector_Correlations.Rmd: This is the script for all vector correlation analyses for the paper 'Genomic architecture of the developing forelegs of Drosophila prolongata; an exaggerated weapon and ornament'. It reads in the data generated from the `DEG_analysis.Rmd` script and generates Figures 8, 9, and 10 as well as supplemental figures S11-S29. It also uses publically available file from Flybase.org `pathway_group_data_fb_2024_04.tsv` available from [https://s3ftp.flybase.org/releases/FB2024_04/precomputed_files/genes/index.html], as well as the publically available mod-encode file for putative grn targets available at [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8109144] and converted to gene names, provided here as file `grn_targetGenes.txt`.

  -Luecke_cell_data_analysis_ID.Rmd: This is the script for the reanalysis of cell numbers for the paper 'Genomic architecture of the developing forelegs of Drosophila prolongata; an exaggerated weapon and ornament'. This script will take cell number and size data from Luecke & Kopp 2019 available via Dryad at [https://datadryad.org/dataset/doi:10.25338/B8303J] and model it to demonstrate that cell number is not a confounding factor for our gene expression results. This script outputs supplemental figure S1.


-bash_scripts

  -QC_trimming: These scripts are associated with step 1 below and were used to trim adapters, and verify that sequencing was sufficient, adapters were removed, and no issues were present in the original sequencing data prior to mapping.

  -merging: These are the scripts associasted with step 2 below and were used to merge two distinct sequenciong runs used to achieve sufficient coverage.

  -blasting_between_spp: These scripts are all the scripts associated with step 3 below and were used to look for contamination and verification that our sequencing was correct.

  -mapping: These are the scripts associated with step 4 below and are used for indexing genomes for STAR and mapping/counting for our sequencing. It outputs count data files that are then read in to the scripts `[species]_model.R` to output modelled count estimates.

  -modelling: These scripts are associated with step 6 below and use the output from the STAR mapping scripts to model and output emmeans estimates which are then read in to `DGE_analysis.Rmd` for estimate shrinkage and downstream analyses.

data:
This is data used to generate results and outputs.

  -PenGAL4 files are phenotypic data for the PenGAL4 .Rmd

  -Count files are raw count data (reads mapped to the genome using STAR) used in DEG analysis .Rmd to generate the ashr corrected files, generated from the STAR mapping shell scripts

  -norms files are ashr regularized measurement used for analysis

For complete analyses, the output files from star ending with 'ReadsPerGene.out.tab' should be splaced in a folder for each species here, however due to the size of these files they are not included here but generated by running the `star_[species].sh` scripts in `/scripts/bash_scripts/mapping/`.

## 1. QC_Trimming

As we used paired-end sequencing, we have 226 R1 and 226 R2 files. First I run fastqc pre-trim - things look fine, adapters found.

I trim with `bbduk.sh`.

A second pass of fastqc also looks good. 

## 2. merging

Next I split all reads from the first run in to their own directory called run1 and the second run in to a directory called run2, while removing the beginning of the name added by the sequencing centre. I then merge these in to a single file using: `merge.sh`

## 3. blasting_between_spp

Before anything, I did a quick check by blasting some sequences from each fastq using the online blast tool. I identified two samples that appeared to be blasting to incorrect species (a melanogaster to a carrolli and the reverse), so I check all samples to make sure they blast best to the species which they are labelled as.

`make_blast_db.sh`

`make_fasta_subset.sh`

`blastin.sh`

I can go through the outputted 'counts' file to identify samples that are blasting to the wrong species.

mel_M_mid_early_3 blasts to carrolli and prolongata equally and much better than melanogaster. I am going to leave this sample aside for the time being.
car_M_mid_early_6 blasts best to melanogaster so I will also leave this sample aside for now.

* A first run of the analysis identified pro_M_mid_late_2 and car_M_mid_late_2 as being outliers in PCR and appear to have strange extreme gene expression that doesn't match any other sample suggesting possible tissue contamination, so I set these aside as well.

## 4. mapping
scaffold names don't match between the genomes and annotations for prolongata and carrolli, so I need to fix that and make names match for mapping.

`sed 's/.*Scaffold/\>Scaffold/g' GCA_036346975.1_ASM3634697v1_genomic.fna | sed 's/\,\ whole\ genome\ shotgun\ sequence//g' > prolongata_genome.fa`

`sed 's/.*contig/\>contig/g' GCA_018152295.1_ASM1815229v1_genomic.fna | sed 's/\,\ whole\ genome\ shotgun\ sequence//g' > carrolli_genome.fa`

I need to convert my gff to a gtf so that STAR maps properly. Previous attempts at mapping to the gff directly resulted in very few hits, and genes not mapping at all that map when a gtf is used. The writer of star has suggested conversion when issues have been brought up [https://github.com/alexdobin/STAR/issues/901]

I attempted this step with both gffread and aagat, gffread appears to not convert the D. prolongata gff properly resulting in very few genes mapping, agat however results in similar mapping numbers as D. melanogaster to the well established FlyBase gtf annotation. So going forward I use the agat converted gff -> gtf

`agat_gff2gtf_pro.sh`

`agat_gff2gtf_car.sh`

and then I index each species genome and annotation.

`star_index_pro.sh`

`star_index_car.sh`

`star_index_mel.sh`

mapping was run in the 2-pass method recomended in the manual. The scripts were run once without splice site annotations, and then ran a second time with `--sjdbFileChrStartEnd /path/to/sjdbFile.txt.`

```
star_pro.sh

star_car.sh

star_mel.sh  

```

The 2-pass mapping really didn't seem to make much of a difference. But it didn't hurt anything, so on we go.

## 5. converting_geneNames

I need a way to convert from the gene numbers in the prolongata/carrolli annotations, to the geneIDs that are orthologous with melanogaster, so I extract the information I need from the annotation for converting between all of these relevant things.

I run:

`orthoswapper_pro.r` via `orthoswapper_pro.sh`

`orthoswapper_car.r` via `orthoswapper_car.sh`

## 6. modelling

Next, I model each species separately and output emmeans contrasts and model estimates for each gene.

`pro_modelling.R` via `pro_model.sh`

`car_modelling.R` via `car_model.sh`

`mel_modelling.R` via `mel_model.sh`

# Rscripts for analyses

## 1. PCA
This script generates Figures 2 and 3, and corresponds to the results section "Expression profiles show species and sex-specific clustering, with little evidence for clustering based on specific leg imaginal disc, or developmental stage"

`PCA_work.Rmd`

## 2. Magnitude and gene number relationship with SSD
This script generates Figure 4 and corresponds to the results section  "Sex-biased gene expression in D. prolongata forelegs is species-specific"

`MagnitudesDistancesBetweenVectors.Rmd`

## 3. Differential expression analysis
This script generates Figure 5 and corresponds to the results section "Sex-biased gene expression in D. prolongata forelegs is species-specific"

`DEG_analysis.Rmd`

## 4. Functional analysis
This scripts generates Figure 7 and corresponds to the results section "RNAi knockdowns of a subset of our candidates shows changes to femur width and length as well as a qualitatively D. prolongata-like phenotype in D. melanogaster femurs"

`PenGAL4_cross.Rmd`

## 5. Vector correlation analyses
This script generates Figures 8, 9, and 10 and corresponds to the results section "Signalling pathways and candidate sex-biased genes appear to be expressed in similar direction and magnitude between all three species" and "Sex-biased expression in key signalling pathways are broadly similar in direction and magnitude between all three species"

`Vector_Correlations.Rmd`

## 6. DeSeq checks
This script corresponds with Supplemental Figures 1 and 2 and is a sanity check to ensure that our methods are roughly in line with the standard DESeq pipeline methods.
`DeSeq_analysis.Rmd`
