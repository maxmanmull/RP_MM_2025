#!/bin/bash

mkdir -p singularity_containers
cd singularity_containers

singularity pull --name trim-galore_0.6.10--hdfd78af_0.sif docker://quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0
singularity pull --name spades_3.15.5--h95f258a_0.sif docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0
singularity pull --name checkm-genome_1.2.4--pyhdfd78af_2.sif docker://quay.io/biocontainers/checkm-genome:1.2.4--pyhdfd78af_2
singularity pull --name bwa_0.7.17--hed695b0_7.sif docker://quay.io/biocontainers/bwa:0.7.17--hed695b0_7
singularity pull --name samtools_1.16.1--h6899075_1.sif docker://quay.io/biocontainers/samtools:1.16.1--h6899075_1
singularity pull --name concoct_1.1.0--py312h71dcd68_7.sif docker://quay.io/biocontainers/concoct:1.1.0--py312h71dcd68_7
singularity pull --name gtdbtk_2.4.1--pyhdfd78af_1.sif docker://quay.io/biocontainers/gtdbtk:2.4.1--pyhdfd78af_1
singularity pull --name quast_5.2.0--py39pl5321h2add14b_1.sif docker://quay.io/biocontainers/quast:5.2.0--py39pl5321h2add14b_1
singularity pull --name busco_v6.0.0_cv1.sif docker://quay.io/biocontainers/busco:6.0.0--pyhdfd78af_0
singularity pull --name dbcan_5.1.2--pyhdfd78af_0.sif docker://quay.io/biocontainers/dbcan:5.1.2--pyhdfd78af_0
singularity pull --name gapseq-1.4.0.sif docker://quay.io/biocontainers/gapseq:1.4.0--h9ee0642_1
singularity pull --name standalone_8.0.2.sif docker://antismash/standalone:8.0.2
singularity pull --name iqtree_3.0.1.sif docker://quay.io/biocontainers/iqtree:3.0.1--h503566f_0
