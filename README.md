Metabolic Modeling Pipeline
Simple setup guide for running the metagenomic analysis pipeline.

1.
Setup Environment:
conda env create -f conda_env.yaml
conda activate conda_metabolic_pipeline

2. Edit config.yaml
Change these paths to match your system:

directories:
  singularity: "/path/to/your/singularity/containers"  # Container images location
  data: "/path/to/your/fastq/files"                    # Raw FASTQ files directory
  results: "results"                                    # Output directory
  media: "/path/to/gapseq/media/files"                 # Gapseq media files
  input: "input"                                        # Input metadata directory

databases:
  gtdbtk: "/path/to/gtdbtk/database"                   # GTDB-Tk database
  checkm: "/path/to/checkm/database"                   # CheckM database
  busco: "/path/to/busco/database"                     # BUSCO database
  dbcan: "/path/to/dbcan/database"                     # dbCAN database

input_files:
  metadata: "input/meta_input_pipeline.csv"            # sample metadata file

3. Prepare Input Files
Create input/meta_input_pipeline.csv with samples:
sample_id,Taurine,Creatinine,Carnitine,Xylan,Chitin
Sample_001,Sample_001,,Sample_001,,
Sample_002,,Sample_002,,Sample_002,Sample_002
Reference_MAG.fa,Reference_MAG.fa,Reference_MAG.fa,Reference_MAG.fa,Reference_MAG.fa,Reference_MAG.fa

Use sample IDs that match FASTQ file names
For reference genomes, add .fa extension
Leave empty cells for media conditions not tested

4. Run the pipeline
snakemake -n

#full run
snakemake --cores 8 --use-singularity

#specific parts
snakemake --cores 8 --use-singularity results/gtdbtk/done.txt  # Only taxonomy
snakemake --cores 8 --use-singularity results/gapseq_models/multi_media_complete.txt  # Only metabolic models

5. Python Analysis Scripts
These are way less modular and req more changing in the hardcode if needed (specificly look at the mappings)
in general these filepaths have to be edited as well
PATHWAY_FILE = "./results/gapseq_pathways/meta_pwy.tbl"  # From gapseq database
GAPSEQ_MODELS_DIR = "./results/gapseq_models"
QUALITY_FILE = "./results/qc_summary/mag_quality_master.tsv"
TAXONOMY_FILE = "./results/gtdbtk/gtdbtk.bac120.summary.tsv"

to run just use python{scriptname}.py

#####
Required External Files
Download these from the gapseq GitHub repository:
meta_pwy.tbl - Pathway metadata
meta_rea_pwy-gapseq.tbl - Reaction-pathway mapping
+*All the databases as above 

Get from: https://github.com/jotech/gapseq/tree/master/dat


FASTQ files should be named: {sample}_R1_001.fastq.gz and {sample}_R2_001.fastq.gz
Reference genomes go in results/all_genomes/ as .fa files
All paths in Python scripts are relative to where you run them from
Run Python scripts from the main project directory



