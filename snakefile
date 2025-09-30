import pandas as pd
from pathlib import Path
import glob
import os

# ============================================================================
# LOAD CONFIGURATION
# ============================================================================

configfile: "config.yaml"

# ============================================================================
# EXTRACT CONFIGURATION VALUES
# ============================================================================

# Directories
SINGULARITY_DIR = config["directories"]["singularity"]
DATA_DIR = config["directories"]["data"]
RESULTS_DIR = config["directories"]["results"]
MEDIA_DIR = config["directories"]["media"]
INPUT_DIR = config["directories"]["input"]

# Databases
GTDBTK_DB = config["databases"]["gtdbtk"]
DBCAN_DB = config["databases"]["dbcan"]
BUSCO_DB = config["databases"]["busco"]
CHECKM_DB = config["databases"]["checkm"]

# Media types
MEDIA_COLS = config["media_types"]

# ============================================================================
# SAMPLE PARSING
# ============================================================================

# Parse metadata to get FASTQ and FASTA samples
metadata = pd.read_csv(config["input_files"]["metadata"])
FASTQ_SAMPLES = []
FASTA_SAMPLES = []

# Extract unique sample IDs from all media columns
for col in MEDIA_COLS:
    if col in metadata.columns:
        samples = metadata[col].dropna()
        for s in samples:
            s_str = str(s)
            if s_str.endswith('.fa'):
                if s_str not in FASTA_SAMPLES:
                    FASTA_SAMPLES.append(s_str)
            else:
                # Handle float conversion (459.0 -> 459)
                if '.' in s_str:
                    s_str = str(int(float(s_str)))
                if s_str not in FASTQ_SAMPLES:
                    FASTQ_SAMPLES.append(s_str)

print(f"Found {len(FASTQ_SAMPLES)} FASTQ samples: {FASTQ_SAMPLES}")
print(f"Found {len(FASTA_SAMPLES)} FASTA samples: {FASTA_SAMPLES}")

# ============================================================================
# MASTER RULE
# ============================================================================

rule all:
    input:
        # Stage 1: Read processing and assembly
        expand(f"{RESULTS_DIR}/trimmed/{{sample}}_R1_001_val_1.fq.gz", sample=FASTQ_SAMPLES),
        expand(f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta", sample=FASTQ_SAMPLES),

        # Stage 2: Contamination assessment and conditional binning
        expand(f"{RESULTS_DIR}/contamination_check/{{sample}}/binning_decision.txt", sample=FASTQ_SAMPLES),
        expand(f"{RESULTS_DIR}/final_genomes/{{sample}}/genome_list.txt", sample=FASTQ_SAMPLES),

        # Stage 3: Genome collection and analysis
        f"{RESULTS_DIR}/all_genomes/collected.txt",
        f"{RESULTS_DIR}/gtdbtk/done.txt",
        f"{RESULTS_DIR}/gtdbtk/phylogenetic_tree/tree.nwk",

        # Stage 4: Functional analyses
        f"{RESULTS_DIR}/gapseq_models/multi_media_complete.txt",
        f"{RESULTS_DIR}/antismash/all_complete.txt",
        f"{RESULTS_DIR}/cazymes/all_done.txt",
        f"{RESULTS_DIR}/busco/busco_summary.tsv",

        # Stage 5: Final report
        f"{RESULTS_DIR}/final_report/genome_quality_report.tsv",

# ============================================================================
# STAGE 1: READ PROCESSING AND ASSEMBLY
# ============================================================================

rule trim_reads:
    input:
        r1 = f"{DATA_DIR}/{{sample}}_R1_001.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        r1 = f"{RESULTS_DIR}/trimmed/{{sample}}_R1_001_val_1.fq.gz",
        r2 = f"{RESULTS_DIR}/trimmed/{{sample}}_R2_001_val_2.fq.gz"
    params:
        outdir = f"{RESULTS_DIR}/trimmed",
        quality = config["parameters"]["trim_galore"]["quality"],
        length = config["parameters"]["trim_galore"]["length"]
    threads: config["parameters"]["trim_galore"]["threads"]
    container: os.path.join(SINGULARITY_DIR, config["containers"]["trim_galore"])
    shell:
        """
        mkdir -p {params.outdir}

        trim_galore \
            --quality {params.quality} \
            --length {params.length} \
            --paired \
            --gzip \
            --cores {threads} \
            --output_dir {params.outdir} \
            {input.r1} {input.r2}
        """

rule assemble_spades:
    input:
        r1 = f"{RESULTS_DIR}/trimmed/{{sample}}_R1_001_val_1.fq.gz",
        r2 = f"{RESULTS_DIR}/trimmed/{{sample}}_R2_001_val_2.fq.gz"
    output:
        contigs = f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta",
        scaffolds = f"{RESULTS_DIR}/assembly/{{sample}}/scaffolds.fasta"
    params:
        outdir = f"{RESULTS_DIR}/assembly/{{sample}}",
        kmers = config["parameters"]["spades"]["kmers"],
        memory = config["parameters"]["spades"]["memory"]
    threads: config["parameters"]["spades"]["threads"]
    container: os.path.join(SINGULARITY_DIR, config["containers"]["spades"])
    shell:
        """
        spades.py \
            --isolate \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            -k {params.kmers} \
            -t {threads} \
            -m {params.memory}
        """

# ============================================================================
# STAGE 2: CONTAMINATION ASSESSMENT
# ============================================================================

rule checkm_assembly_contamination:
    input:
        contigs = f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta"
    output:
        contamination_check = f"{RESULTS_DIR}/contamination_check/{{sample}}/contamination_assessment.txt",
        checkm_results = f"{RESULTS_DIR}/contamination_check/{{sample}}/checkm_results.txt",
        checkm_table = f"{RESULTS_DIR}/contamination_check/{{sample}}/checkm_table.tsv",
        decision = f"{RESULTS_DIR}/contamination_check/{{sample}}/binning_decision.txt"
    params:
        outdir = f"{RESULTS_DIR}/contamination_check/{{sample}}",
        checkm_container = os.path.join(SINGULARITY_DIR, config["containers"]["checkm"]),
        checkm_data = CHECKM_DB,
        contamination_threshold = config["parameters"]["checkm"]["contamination_threshold"],
        completeness_threshold = config["parameters"]["checkm"]["completeness_threshold"]
    threads: config["parameters"]["checkm"]["threads"]
    shell:
        """
        export CHECKM_DATA_PATH={params.checkm_data}

        mkdir -p {params.outdir}/temp_genome

        # Create temporary directory with assembly as single "bin"
        cp {input.contigs} {params.outdir}/temp_genome/assembly.fa

        # Run CheckM on the whole assembly
        singularity exec -B {params.checkm_data} {params.checkm_container} \
            checkm lineage_wf \
            -t {threads} \
            -x fa \
            {params.outdir}/temp_genome \
            {params.outdir}/checkm_output \
            > {output.checkm_results}

        # Generate table format
        singularity exec -B {params.checkm_data} {params.checkm_container} \
            checkm qa \
            -o 2 \
            -f {output.checkm_table} \
            --tab_table \
            {params.outdir}/checkm_output/lineage.ms \
            {params.outdir}/checkm_output

        # Extract completeness and contamination
        comp=$(tail -n 1 {output.checkm_table} | cut -f6)
        cont=$(tail -n 1 {output.checkm_table} | cut -f7)

        echo "Assembly Contamination Assessment for {wildcards.sample}" > {output.contamination_check}
        echo "Completeness: $comp%" >> {output.contamination_check}
        echo "Contamination: $cont%" >> {output.contamination_check}

        # Decision logic
        if awk -v cont="$cont" -v thresh="{params.contamination_threshold}" 'BEGIN {{exit !(cont > thresh)}}'; then
            echo "DECISION: RUN_BINNING" > {output.decision}
            echo "Reason: Contamination ($cont%) exceeds threshold ({params.contamination_threshold}%)" >> {output.decision}
        elif awk -v comp="$comp" 'BEGIN {{exit !(comp < 50)}}'; then
            echo "DECISION: RUN_BINNING" > {output.decision}
            echo "Reason: Low completeness ($comp%) suggests multiple partial genomes" >> {output.decision}
        else
            echo "DECISION: USE_ASSEMBLY" > {output.decision}
            echo "Reason: Clean assembly (Contamination: $cont%, Completeness: $comp%)" >> {output.decision}
        fi

        cat {output.decision} >> {output.contamination_check}
        """

# ============================================================================
# STAGE 3: CONDITIONAL BINNING
# ============================================================================

rule map_reads_for_binning:
    input:
        contigs = f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta",
        r1 = f"{RESULTS_DIR}/trimmed/{{sample}}_R1_001_val_1.fq.gz",
        r2 = f"{RESULTS_DIR}/trimmed/{{sample}}_R2_001_val_2.fq.gz",
        decision = f"{RESULTS_DIR}/contamination_check/{{sample}}/binning_decision.txt"
    output:
        bam = f"{RESULTS_DIR}/mapping/{{sample}}/{{sample}}.sorted.bam",
        bai = f"{RESULTS_DIR}/mapping/{{sample}}/{{sample}}.sorted.bam.bai"
    params:
        outdir = f"{RESULTS_DIR}/mapping/{{sample}}",
        bwa_container = os.path.join(SINGULARITY_DIR, config["containers"]["bwa"]),
        samtools_container = os.path.join(SINGULARITY_DIR, config["containers"]["samtools"])
    threads: config["parameters"]["mapping"]["threads"]
    shell:
        """
        # Only run mapping if binning is needed
        if grep -q "RUN_BINNING" {input.decision}; then
            mkdir -p {params.outdir}

            # Index contigs
            singularity exec {params.bwa_container} \
                bwa index {input.contigs}

            # Map reads
            singularity exec {params.bwa_container} \
                bwa mem -t {threads} {input.contigs} {input.r1} {input.r2} | \
            singularity exec {params.samtools_container} \
                samtools sort -@ {threads} -o {output.bam}

            # Index BAM
            singularity exec {params.samtools_container} \
                samtools index -@ {threads} {output.bam}
        else
            # Create empty files if not binning
            mkdir -p {params.outdir}
            touch {output.bam}
            touch {output.bai}
        fi
        """

rule process_based_on_contamination:
    input:
        contigs = f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta",
        decision = f"{RESULTS_DIR}/contamination_check/{{sample}}/binning_decision.txt",
        bam = f"{RESULTS_DIR}/mapping/{{sample}}/{{sample}}.sorted.bam",
        contamination_check = f"{RESULTS_DIR}/contamination_check/{{sample}}/contamination_assessment.txt"
    output:
        genome_dir = directory(f"{RESULTS_DIR}/final_genomes/{{sample}}"),
        genome_list = f"{RESULTS_DIR}/final_genomes/{{sample}}/genome_list.txt",
        summary = f"{RESULTS_DIR}/final_genomes/{{sample}}/processing_summary.txt"
    params:
        concoct_container = os.path.join(SINGULARITY_DIR, config["containers"]["concoct"]),
        checkm_container = os.path.join(SINGULARITY_DIR, config["containers"]["checkm"]),
        chunk_size = config["parameters"]["concoct"]["chunk_size"],
        overlap = config["parameters"]["concoct"]["overlap"]
    threads: config["parameters"]["concoct"]["threads"]
    shell:
        """
        mkdir -p {output.genome_dir}

        # Read decision
        if grep -q "USE_ASSEMBLY" {input.decision}; then
            echo "=== Using full assembly as single genome ===" > {output.summary}

            # Copy assembly as single genome
            cp {input.contigs} {output.genome_dir}/{wildcards.sample}_genome.fa
            echo "{wildcards.sample}_genome.fa" > {output.genome_list}

            echo "" >> {output.summary}
            cat {input.contamination_check} >> {output.summary}
            echo "Processing: Skipped binning - clean assembly" >> {output.summary}

        else
            echo "=== Running binning to separate contaminating organisms ===" > {output.summary}

            # Run CONCOCT binning
            mkdir -p {output.genome_dir}/binning_temp

            # Cut up contigs
            singularity exec {params.concoct_container} \
                cut_up_fasta.py {input.contigs} -c {params.chunk_size} -o {params.overlap} --merge_last \
                -b {output.genome_dir}/binning_temp/contigs_10K.bed \
                > {output.genome_dir}/binning_temp/contigs_10K.fa

            # Generate coverage table
            singularity exec {params.concoct_container} \
                concoct_coverage_table.py {output.genome_dir}/binning_temp/contigs_10K.bed \
                {input.bam} > {output.genome_dir}/binning_temp/coverage_table.tsv

            # Run CONCOCT
            singularity exec {params.concoct_container} \
                concoct --composition_file {output.genome_dir}/binning_temp/contigs_10K.fa \
                --coverage_file {output.genome_dir}/binning_temp/coverage_table.tsv \
                -b {output.genome_dir}/binning_temp/ \
                -t {threads}

            # Extract bins
            singularity exec {params.concoct_container} \
                merge_cutup_clustering.py {output.genome_dir}/binning_temp/clustering_gt1000.csv \
                > {output.genome_dir}/binning_temp/clustering_merged.csv
            mkdir -p {output.genome_dir}/binning_temp/bins
            singularity exec {params.concoct_container} \
                extract_fasta_bins.py {input.contigs} \
                {output.genome_dir}/binning_temp/clustering_merged.csv \
                --output_path {output.genome_dir}/binning_temp/bins

            # Run CheckM on bins
            mkdir -p {output.genome_dir}/binning_temp/checkm_output
            singularity exec {params.checkm_container} \
                checkm lineage_wf \
                -t {threads} \
                -x fa \
                {output.genome_dir}/binning_temp/bins \
                {output.genome_dir}/binning_temp/checkm_output \
                > {output.genome_dir}/binning_temp/checkm_results.txt

            # Generate CheckM table
            singularity exec {params.checkm_container} \
                checkm qa \
                -o 2 \
                -f {output.genome_dir}/binning_temp/checkm_table.tsv \
                --tab_table \
                {output.genome_dir}/binning_temp/checkm_output/lineage.ms \
                {output.genome_dir}/binning_temp/checkm_output

            # Select high-quality bins
            echo "Selecting high-quality bins..." >> {output.summary}
            > {output.genome_list}

            # Parse CheckM table (skip header)
            tail -n +2 {output.genome_dir}/binning_temp/checkm_table.tsv | while IFS=$'\t' read -r bin_id marker lineage genomes markers0 comp cont strain_het rest; do
                if [ ! -z "$bin_id" ]; then
                    # Check if bin passes thresholds (>90% complete, <5% contamination)
                    if awk -v comp="$comp" -v cont="$cont" 'BEGIN {{exit !(comp > 90 && cont < 5)}}'; then
                        cp {output.genome_dir}/binning_temp/bins/${{bin_id}}.fa \
                           {output.genome_dir}/{wildcards.sample}_${{bin_id}}.fa
                        echo "{wildcards.sample}_${{bin_id}}.fa" >> {output.genome_list}
                        echo "  Selected: $bin_id (Comp: $comp%, Cont: $cont%)" >> {output.summary}
                    else
                        echo "  Rejected: $bin_id (Comp: $comp%, Cont: $cont%)" >> {output.summary}
                    fi
                fi
            done

            # If no good bins, use original assembly
            if [ ! -s {output.genome_list} ]; then
                echo "WARNING: No high-quality bins recovered. Using original assembly." >> {output.summary}
                cp {input.contigs} {output.genome_dir}/{wildcards.sample}_genome.fa
                echo "{wildcards.sample}_genome.fa" > {output.genome_list}
            fi

            # Clean up
            rm -rf {output.genome_dir}/binning_temp
        fi

        echo "" >> {output.summary}
        echo "Final genomes in directory:" >> {output.summary}
        ls -lh {output.genome_dir}/*.fa 2>/dev/null | awk '{{print $NF}}' | xargs -n1 basename >> {output.summary}
        """

# ============================================================================
# STAGE 4: COLLECT ALL GENOMES
# ============================================================================

rule collect_all_genomes:
    input:
        genome_lists = expand(f"{RESULTS_DIR}/final_genomes/{{sample}}/genome_list.txt", sample=FASTQ_SAMPLES),
        processing_summaries = expand(f"{RESULTS_DIR}/final_genomes/{{sample}}/processing_summary.txt", sample=FASTQ_SAMPLES)
    output:
        collected = f"{RESULTS_DIR}/all_genomes/collected.txt",
        summary = f"{RESULTS_DIR}/all_genomes/collection_summary.txt",
        genome_list = f"{RESULTS_DIR}/all_genomes/genome_names.txt"
    params:
        output_dir = f"{RESULTS_DIR}/all_genomes"
    run:
        import os
        import shutil

        os.makedirs(params.output_dir, exist_ok=True)

        genome_count = 0
        all_genomes = []
        binned_samples = []
        unbinned_samples = []

        # Collect genomes from FASTQ samples
        for sample in FASTQ_SAMPLES:
            genome_dir = f"{RESULTS_DIR}/final_genomes/{sample}"
            genome_list_file = f"{genome_dir}/genome_list.txt"
            processing_summary = f"{genome_dir}/processing_summary.txt"

            # Check if sample was binned
            if os.path.exists(processing_summary):
                with open(processing_summary, 'r') as f:
                    content = f.read()
                    if "Skipped binning" in content:
                        unbinned_samples.append(sample)
                    else:
                        binned_samples.append(sample)

            # Copy genomes
            if os.path.exists(genome_list_file):
                with open(genome_list_file, 'r') as f:
                    for genome in f:
                        genome = genome.strip()
                        if genome:
                            src = f"{genome_dir}/{genome}"
                            dst = os.path.join(params.output_dir, genome)
                            if os.path.exists(src):
                                shutil.copy(src, dst)
                                genome_count += 1
                                genome_base = genome.replace('.fa', '').replace('.fasta', '')
                                all_genomes.append(genome_base)

        # Copy reference genomes
        for fasta in FASTA_SAMPLES:
            src = os.path.join(INPUT_DIR, fasta)
            dst = os.path.join(params.output_dir, fasta)
            if os.path.exists(src):
                shutil.copy(src, dst)
                genome_count += 1
                genome_base = fasta.replace('.fa', '').replace('.fasta', '')
                all_genomes.append(genome_base)

        # Write genome names for checkpoint
        with open(output.genome_list, 'w') as f:
            for genome in all_genomes:
                f.write(f"{genome}\n")

        # Write summary
        with open(output.summary, 'w') as f:
            f.write("Genome Collection Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total genomes collected: {genome_count}\n\n")
            f.write(f"Samples processed without binning (clean): {len(unbinned_samples)}\n")
            for s in unbinned_samples:
                f.write(f"  - {s}\n")
            f.write(f"\nSamples that required binning (contaminated): {len(binned_samples)}\n")
            for s in binned_samples:
                f.write(f"  - {s}\n")

        with open(output.collected, 'w') as f:
            f.write(f"Collection complete: {genome_count} genomes\n")

# Checkpoint to get genome list dynamically
checkpoint get_genome_list:
    input:
        collected = f"{RESULTS_DIR}/all_genomes/collected.txt"
    output:
        genome_list = f"{RESULTS_DIR}/all_genomes/checkpoint_genome_names.txt"
    shell:
        """
        ls {RESULTS_DIR}/all_genomes/*.fa 2>/dev/null | xargs -n1 basename | sed 's/.fa$//' > {output.genome_list}
        """

# ============================================================================
# STAGE 5: GTDB-Tk CLASSIFICATION
# ============================================================================

rule gtdbtk_classify:
    input:
        collected = f"{RESULTS_DIR}/all_genomes/collected.txt"
    output:
        summary = f"{RESULTS_DIR}/gtdbtk/gtdbtk.bac120.summary.tsv",
        done = f"{RESULTS_DIR}/gtdbtk/done.txt"
    params:
        genome_dir = f"{RESULTS_DIR}/all_genomes",
        outdir = f"{RESULTS_DIR}/gtdbtk",
        gtdbtk_db = GTDBTK_DB,
        gtdbtk_container = os.path.join(SINGULARITY_DIR, config["containers"]["gtdbtk"]),
        pplacer_cpus = config["parameters"]["gtdbtk"]["pplacer_cpus"]
    threads: config["parameters"]["gtdbtk"]["threads"]
    shell:
        """
        export GTDBTK_DATA_PATH={params.gtdbtk_db}

        genome_count=$(ls {params.genome_dir}/*.fa 2>/dev/null | wc -l)

        if [ "$genome_count" -gt 0 ]; then
            singularity exec -B {params.gtdbtk_db} {params.gtdbtk_container} \
                gtdbtk classify_wf \
                --genome_dir {params.genome_dir} \
                --extension fa \
                --out_dir {params.outdir} \
                --cpus {threads} \
                --skip_ani_screen \
                --pplacer_cpus {params.pplacer_cpus}
        else
            echo "No genomes to classify" > {output.summary}
        fi

        touch {output.done}
        """

rule build_phylogenetic_tree:
    input:
        done = f"{RESULTS_DIR}/gtdbtk/done.txt"  # Wait for GTDB-Tk to finish
    output:
        tree = f"{RESULTS_DIR}/gtdbtk/phylogenetic_tree/tree.nwk",
        tree_ml = f"{RESULTS_DIR}/gtdbtk/phylogenetic_tree/tree.treefile",
        log = f"{RESULTS_DIR}/gtdbtk/phylogenetic_tree/tree.log"
    params:
        msa = f"{RESULTS_DIR}/gtdbtk/align/gtdbtk.bac120.user_msa.fasta.gz",
        outdir = f"{RESULTS_DIR}/gtdbtk/phylogenetic_tree",
        iqtree_container = os.path.join(SINGULARITY_DIR, config["containers"]["iqtree"]),
        bootstrap = config["parameters"]["iqtree"]["bootstrap"],
        alrt = config["parameters"]["iqtree"]["alrt"]
    threads: config["parameters"]["iqtree"]["threads"]
    shell:
        """
        mkdir -p {params.outdir}

        # Check if MSA exists
        if [ -f {params.msa} ] && [ -s {params.msa} ]; then
            gunzip -c {params.msa} > {params.outdir}/msa.fasta

            singularity exec {params.iqtree_container} \
                iqtree \
                -s {params.outdir}/msa.fasta \
                -pre {params.outdir}/tree \
                -m MFP \
                -bb {params.bootstrap} \
                -alrt {params.alrt} \
                -nt {threads}

            if [ -f {params.outdir}/tree.treefile ]; then
                cp {params.outdir}/tree.treefile {output.tree}
            else
                echo "Tree construction failed" > {output.tree}
                touch {output.tree_ml}
                echo "Tree construction failed" > {output.log}
            fi

            rm -f {params.outdir}/msa.fasta
        else
            echo "No MSA available - GTDB-Tk may have found only one genome or failed" > {output.tree}
            echo "No MSA available" > {output.tree_ml}
            echo "GTDB-Tk did not produce an MSA file - possibly only one genome or classification failed" > {output.log}
        fi
        """
# ============================================================================
# STAGE 6: GAPSEQ METABOLIC RECONSTRUCTION
# ============================================================================

rule gapseq_pathways_single:
    input:
        genome = f"{RESULTS_DIR}/all_genomes/{{genome}}.fa"
    output:
        pathways = f"{RESULTS_DIR}/gapseq_pathways/{{genome}}-all-Pathways.tbl",
        reactions = f"{RESULTS_DIR}/gapseq_pathways/{{genome}}-all-Reactions.tbl"
    params:
        outdir = f"{RESULTS_DIR}/gapseq_pathways",
        gapseq_container = os.path.join(SINGULARITY_DIR, config["containers"]["gapseq"])
    threads: config["parameters"]["gapseq"]["threads"]
    shell:
        """
        cd {params.outdir}
        singularity exec {params.gapseq_container} \
            gapseq find -p all ../../{input.genome}
        """

rule gapseq_transporters_single:
    input:
        genome = f"{RESULTS_DIR}/all_genomes/{{genome}}.fa"
    output:
        transporters = f"{RESULTS_DIR}/gapseq_transporters/{{genome}}-Transporter.tbl"
    params:
        outdir = f"{RESULTS_DIR}/gapseq_transporters",
        gapseq_container = os.path.join(SINGULARITY_DIR, config["containers"]["gapseq"])
    threads: config["parameters"]["gapseq"]["threads"]
    shell:
        """
        cd {params.outdir}
        singularity exec {params.gapseq_container} \
            gapseq find-transport ../../{input.genome}
        """

rule gapseq_draft_single:
    input:
        genome = f"{RESULTS_DIR}/all_genomes/{{genome}}.fa",
        pathways = f"{RESULTS_DIR}/gapseq_pathways/{{genome}}-all-Pathways.tbl",
        reactions = f"{RESULTS_DIR}/gapseq_pathways/{{genome}}-all-Reactions.tbl",
        transporters = f"{RESULTS_DIR}/gapseq_transporters/{{genome}}-Transporter.tbl"
    output:
        draft = f"{RESULTS_DIR}/gapseq_draft/{{genome}}-draft.RDS",
        weights = f"{RESULTS_DIR}/gapseq_draft/{{genome}}-rxnWeights.RDS",
        genes = f"{RESULTS_DIR}/gapseq_draft/{{genome}}-rxnXgenes.RDS"
    params:
        outdir = f"{RESULTS_DIR}/gapseq_draft",
        gapseq_container = os.path.join(SINGULARITY_DIR, config["containers"]["gapseq"])
    threads: config["parameters"]["gapseq"]["threads"]
    shell:
        """
        cd {params.outdir}
        singularity exec {params.gapseq_container} \
            gapseq draft \
            -r ../../{input.reactions} \
            -t ../../{input.transporters} \
            -p ../../{input.pathways} \
            -c ../../{input.genome}
        """

# Checkpoint to determine genome-media combinations
checkpoint get_genome_media_map:
    input:
        collected = f"{RESULTS_DIR}/all_genomes/collected.txt"
    output:
        genome_media_map = f"{RESULTS_DIR}/gapseq_models/genome_media_map.txt"
    run:
        import os
        import glob

        metadata = pd.read_csv(config["input_files"]["metadata"])
        genome_media_pairs = []

        genome_dir = f"{RESULTS_DIR}/all_genomes"
        existing_genomes = [os.path.basename(f).replace('.fa', '')
                           for f in glob.glob(f"{genome_dir}/*.fa")]

        for genome in existing_genomes:
            genome_media = set()

            if genome == "Marinacidobacteraceae":
                for media in MEDIA_COLS:
                    if media in metadata.columns:
                        samples = metadata[media].dropna()
                        if any(str(s).endswith('.fa') for s in samples):
                            genome_media.add(media)
            else:
                sample_id = genome.split('_')[0]

                for media in MEDIA_COLS:
                    if media in metadata.columns:
                        samples = metadata[media].dropna()
                        for s in samples:
                            s_str = str(s)
                            if not s_str.endswith('.fa'):
                                if '.' in s_str:
                                    s_str = str(int(float(s_str)))
                                if s_str == sample_id:
                                    genome_media.add(media)

            for media in genome_media:
                genome_media_pairs.append(f"{genome}\t{media}")

        os.makedirs(f"{RESULTS_DIR}/gapseq_models", exist_ok=True)
        with open(output.genome_media_map, 'w') as f:
            if genome_media_pairs:
                for pair in genome_media_pairs:
                    f.write(f"{pair}\n")
            else:
                f.write("# No genome-media combinations found\n")

rule create_dummy_model:
    output:
        f"{RESULTS_DIR}/gapseq_models/dummy.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/gapseq_models
        echo "No models to create" > {output}
        """

rule gapseq_fill_multi_media:
    input:
        draft = f"{RESULTS_DIR}/gapseq_draft/{{genome}}-draft.RDS",
        weights = f"{RESULTS_DIR}/gapseq_draft/{{genome}}-rxnWeights.RDS",
        genes = f"{RESULTS_DIR}/gapseq_draft/{{genome}}-rxnXgenes.RDS",
        media = f"{MEDIA_DIR}/{{media}}_media.csv"
    output:
        model = f"{RESULTS_DIR}/gapseq_models/{{genome}}_{{media}}.RDS"
    params:
        outdir = f"{RESULTS_DIR}/gapseq_models",
        gapseq_container = os.path.join(SINGULARITY_DIR, config["containers"]["gapseq"])
    threads: config["parameters"]["gapseq"]["threads"]
    shell:
        """
        cd {params.outdir}

        singularity exec {params.gapseq_container} \
            gapseq fill \
            -m ../../{input.draft} \
            -c ../../{input.weights} \
            -g ../../{input.genes} \
            -n {input.media} \
            -o {wildcards.genome}

        if [ -f "{wildcards.genome}.RDS" ]; then
            mv {wildcards.genome}.RDS {wildcards.genome}_{wildcards.media}.RDS
        else
            echo "ERROR: Expected output file {wildcards.genome}.RDS not found"
            exit 1
        fi

        if [ -f "{wildcards.genome}.xml" ]; then
            mv {wildcards.genome}.xml {wildcards.genome}_{wildcards.media}.xml
        fi
        """

def get_all_genome_media_models(wildcards):
    checkpoint_output = checkpoints.get_genome_media_map.get(**wildcards).output[0]

    models = []
    if os.path.exists(checkpoint_output):
        with open(checkpoint_output) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        genome, media = parts
                        models.append(f"{RESULTS_DIR}/gapseq_models/{genome}_{media}.RDS")

    return models if models else [f"{RESULTS_DIR}/gapseq_models/dummy.txt"]

rule gapseq_multi_media_complete:
    input:
        models = get_all_genome_media_models,
        genome_map = f"{RESULTS_DIR}/gapseq_models/genome_media_map.txt"
    output:
        done = f"{RESULTS_DIR}/gapseq_models/multi_media_complete.txt",
        summary = f"{RESULTS_DIR}/gapseq_models/multi_media_summary.txt"
    run:
        import os

        media_counts = {media: 0 for media in MEDIA_COLS}
        genome_counts = {}
        total_models = 0

        genome_media_pairs = []
        if os.path.exists(input.genome_map):
            with open(input.genome_map, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) == 2:
                            genome, media = parts
                            genome_media_pairs.append((genome, media))

        for genome, media in genome_media_pairs:
            model_path = f"{RESULTS_DIR}/gapseq_models/{genome}_{media}.RDS"
            if os.path.exists(model_path):
                total_models += 1

                if media in media_counts:
                    media_counts[media] += 1

                if genome not in genome_counts:
                    genome_counts[genome] = []
                genome_counts[genome].append(media)

        with open(output.summary, 'w') as f:
            f.write(f"Multi-Media Gap-Filling Summary\n")
            f.write(f"================================\n\n")
            f.write(f"Total models created: {total_models}\n")
            f.write(f"Total genomes: {len(genome_counts)}\n\n")

            f.write("Models per media type:\n")
            for media, count in sorted(media_counts.items()):
                f.write(f"  {media}: {count} models\n")

            f.write(f"\nGenomes with multiple media:\n")
            for genome, media_list in sorted(genome_counts.items()):
                if len(media_list) > 1:
                    f.write(f"  {genome}: {', '.join(media_list)}\n")

        with open(output.done, 'w') as f:
            f.write(f"Completed {total_models} multi-media models\n")

# ============================================================================
# STAGE 7: ANTISMASH ANALYSIS
# ============================================================================

rule antismash_single:
    input:
        genome = f"{RESULTS_DIR}/all_genomes/{{genome}}.fa"
    output:
        done = f"{RESULTS_DIR}/antismash/{{genome}}/done.txt",
        json = f"{RESULTS_DIR}/antismash/{{genome}}/{{genome}}.json",
        gbk = f"{RESULTS_DIR}/antismash/{{genome}}/{{genome}}.gbk"
    params:
        outdir = f"{RESULTS_DIR}/antismash/{{genome}}",
        antismash_container = os.path.join(SINGULARITY_DIR, config["containers"]["antismash"])
    threads: config["parameters"]["antismash"]["threads"]
    shell:
        """
        mkdir -p {params.outdir}

        singularity exec {params.antismash_container} \
            antismash \
            --taxon bacteria \
            --genefinding-tool prodigal \
            --output-dir {params.outdir} \
            --output-basename {wildcards.genome} \
            --cpus {threads} \
            --minimal \
            --enable-genefunctions \
            --enable-lanthipeptides \
            --enable-lassopeptides \
            --enable-nrps-pks \
            --enable-sactipeptides \
            --enable-t2pks \
            --enable-thiopeptides \
            --enable-tta \
            {input.genome}

        touch {output.done}
        """

def aggregate_antismash_outputs(wildcards):
    checkpoint_output = checkpoints.get_genome_list.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        genomes = [line.strip() for line in f if line.strip()]
    return expand(f"{RESULTS_DIR}/antismash/{{genome}}/done.txt", genome=genomes)

rule antismash_summary:
    input:
        antismash_done = aggregate_antismash_outputs,
        gapseq_done = f"{RESULTS_DIR}/gapseq_models/multi_media_complete.txt"
    output:
        summary = f"{RESULTS_DIR}/antismash/bgc_summary.txt",
        done = f"{RESULTS_DIR}/antismash/all_complete.txt"
    params:
        antismash_dir = f"{RESULTS_DIR}/antismash"
    run:
        import json
        import os

        bgc_summary = []
        genome_bgc_counts = {}
        bgc_types = {}

        for done_file in input.antismash_done:
            genome_dir = os.path.dirname(done_file)
            genome_name = os.path.basename(genome_dir)
            json_file = os.path.join(genome_dir, f"{genome_name}.json")

            if os.path.exists(json_file):
                try:
                    with open(json_file, 'r') as f:
                        data = json.load(f)

                    records = data.get('records', [])
                    total_bgcs = 0
                    genome_bgc_types = []

                    for record in records:
                        areas = record.get('areas', [])
                        for area in areas:
                            total_bgcs += 1
                            products = area.get('products', [])
                            for product in products:
                                genome_bgc_types.append(product)
                                if product not in bgc_types:
                                    bgc_types[product] = 0
                                bgc_types[product] += 1

                    genome_bgc_counts[genome_name] = {
                        'total': total_bgcs,
                        'types': genome_bgc_types
                    }

                except Exception as e:
                    print(f"Error parsing {json_file}: {e}")
                    genome_bgc_counts[genome_name] = {'total': 0, 'types': []}
            else:
                genome_bgc_counts[genome_name] = {'total': 0, 'types': []}

        with open(output.summary, 'w') as f:
            f.write("antiSMASH Biosynthetic Gene Cluster Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total genomes analyzed: {len(genome_bgc_counts)}\n")
            f.write(f"Total BGCs found: {sum(g['total'] for g in genome_bgc_counts.values())}\n\n")

            if bgc_types:
                f.write("BGC types found:\n")
                for bgc_type, count in sorted(bgc_types.items(), key=lambda x: x[1], reverse=True):
                    f.write(f"  {bgc_type}: {count}\n")

            f.write("\nBGCs per genome:\n")
            for genome, info in sorted(genome_bgc_counts.items(), key=lambda x: x[1]['total'], reverse=True):
                if info['total'] > 0:
                    f.write(f"  {genome}: {info['total']} BGCs\n")

        with open(output.done, 'w') as f:
            f.write(f"antiSMASH analysis complete for {len(genome_bgc_counts)} genomes\n")

# ============================================================================
# STAGE 8: CAZYME ANALYSIS
# ============================================================================

rule run_dbcan:
    input:
        genome = f"{RESULTS_DIR}/all_genomes/{{genome}}.fa"
    output:
        overview = f"{RESULTS_DIR}/cazymes/{{genome}}/overview.txt",
        done = f"{RESULTS_DIR}/cazymes/{{genome}}/done.txt"
    params:
        outdir = f"{RESULTS_DIR}/cazymes/{{genome}}",
        dbcan_container = os.path.join(SINGULARITY_DIR, config["containers"]["dbcan"]),
        db_dir = DBCAN_DB
    threads: config["parameters"]["dbcan"]["threads"]
    shell:
        """
        mkdir -p {params.outdir}

        singularity exec {params.dbcan_container} \
            run_dbcan \
            {input.genome} \
            prok \
            --out_dir {params.outdir} \
            --db_dir {params.db_dir} \
            --use_signalP=FALSE \
            --dia_cpu {threads} \
            --hmm_cpu {threads}

        touch {output.done}
        """

def aggregate_cazyme_outputs(wildcards):
    checkpoint_output = checkpoints.get_genome_list.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        genomes = [line.strip() for line in f if line.strip()]
    return expand(f"{RESULTS_DIR}/cazymes/{{genome}}/done.txt", genome=genomes)

rule cazymes_complete:
    input:
        aggregate_cazyme_outputs
    output:
        f"{RESULTS_DIR}/cazymes/all_done.txt"
    shell:
        """
        echo "CAZyme analysis complete" > {output}
        echo "Analyzed genomes:" >> {output}
        ls -d {RESULTS_DIR}/cazymes/*/ | wc -l >> {output}
        """

# ============================================================================
# STAGE 9: BUSCO ASSESSMENT
# ============================================================================

rule busco_genome_assessment:
    input:
        genome = f"{RESULTS_DIR}/all_genomes/{{genome}}.fa"
    output:
        done = f"{RESULTS_DIR}/busco/{{genome}}/done.txt",
        summary = f"{RESULTS_DIR}/busco/{{genome}}/{{genome}}/short_summary.specific.{config['parameters']['busco']['lineage']}.{{genome}}.txt",
        full_table = f"{RESULTS_DIR}/busco/{{genome}}/{{genome}}/run_{config['parameters']['busco']['lineage']}/full_table.tsv"
    params:
        outdir = f"{RESULTS_DIR}/busco/{{genome}}",
        busco_container = os.path.join(SINGULARITY_DIR, config["containers"]["busco"]),
        busco_downloads = BUSCO_DB,
        lineage = config["parameters"]["busco"]["lineage"]
    threads: config["parameters"]["busco"]["threads"]
    shell:
        """
        mkdir -p {params.outdir}

        export BUSCO_DOWNLOADS_PATH={params.busco_downloads}
        export BUSCO_CONFIG_FILE={params.busco_downloads}/config.ini

        echo "[busco]" > {params.busco_downloads}/config.ini
        echo "offline = True" >> {params.busco_downloads}/config.ini
        echo "download_path = {params.busco_downloads}" >> {params.busco_downloads}/config.ini

        singularity exec -B {params.busco_downloads} {params.busco_container} \
            busco \
            -i {input.genome} \
            -o {wildcards.genome} \
            --out_path {params.outdir} \
            -m genome \
            -l {params.busco_downloads}/lineages/{params.lineage} \
            -c {threads} \
            -f \
            --offline

        touch {output.done}
        """

def aggregate_busco_outputs(wildcards):
    checkpoint_output = checkpoints.get_genome_list.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        genomes = [line.strip() for line in f if line.strip()]
    return expand(f"{RESULTS_DIR}/busco/{{genome}}/done.txt", genome=genomes)

rule busco_summary:
    input:
        busco_done = aggregate_busco_outputs
    output:
        summary_table = f"{RESULTS_DIR}/busco/busco_summary.tsv"
    params:
        busco_dir = f"{RESULTS_DIR}/busco"
    run:
        import os
        import re
        import glob

        busco_results = []

        for done_file in input.busco_done:
            genome_dir = os.path.dirname(done_file)
            genome_name = os.path.basename(genome_dir)

            summary_files = glob.glob(f"{genome_dir}/short_summary.*.txt")
            if summary_files:
                summary_file = summary_files[0]

                with open(summary_file, 'r') as f:
                    content = f.read()

                    complete_match = re.search(r'(\d+\.?\d*)\s*%.*Complete BUSCOs', content)
                    single_match = re.search(r'(\d+\.?\d*)\s*%.*Complete and single-copy', content)
                    duplicated_match = re.search(r'(\d+\.?\d*)\s*%.*Complete and duplicated', content)
                    fragmented_match = re.search(r'(\d+\.?\d*)\s*%.*Fragmented BUSCOs', content)
                    missing_match = re.search(r'(\d+\.?\d*)\s*%.*Missing BUSCOs', content)
                    total_match = re.search(r'Total BUSCO groups searched\s+(\d+)', content)

                    busco_results.append({
                        'Genome': genome_name,
                        'Complete_pct': float(complete_match.group(1)) if complete_match else 0,
                        'Single_copy_pct': float(single_match.group(1)) if single_match else 0,
                        'Duplicated_pct': float(duplicated_match.group(1)) if duplicated_match else 0,
                        'Fragmented_pct': float(fragmented_match.group(1)) if fragmented_match else 0,
                        'Missing_pct': float(missing_match.group(1)) if missing_match else 0,
                        'Total_BUSCOs': int(total_match.group(1)) if total_match else 0
                    })
            else:
                busco_results.append({
                    'Genome': genome_name,
                    'Complete_pct': 0,
                    'Single_copy_pct': 0,
                    'Duplicated_pct': 0,
                    'Fragmented_pct': 0,
                    'Missing_pct': 0,
                    'Total_BUSCOs': 0
                })

        if busco_results:
            import pandas as pd
            busco_df = pd.DataFrame(busco_results)
            busco_df = busco_df.sort_values('Complete_pct', ascending=False)
            busco_df.to_csv(output.summary_table, sep='\t', index=False)
        else:
            with open(output.summary_table, 'w') as f:
                f.write("Genome\tComplete_pct\tSingle_copy_pct\tDuplicated_pct\tFragmented_pct\tMissing_pct\tTotal_BUSCOs\n")

# ============================================================================
# STAGE 10: QUAST ASSEMBLY STATISTICS (SINGLE RUN)
# ============================================================================

rule quast_all_genomes:
    input:
        genome = f"{RESULTS_DIR}/all_genomes/{{genome}}.fa"
    output:
        report = f"{RESULTS_DIR}/quast/{{genome}}/report.tsv"
    params:
        outdir = f"{RESULTS_DIR}/quast/{{genome}}",
        quast_container = os.path.join(SINGULARITY_DIR, config["containers"]["quast"])
    threads: config["parameters"]["quast"]["threads"]
    shell:
        """
        singularity exec {params.quast_container} \
            quast.py \
            {input.genome} \
            -o {params.outdir} \
            --threads {threads} \
            --min-contig 500 \
            --no-plots \
            --no-html \
            --no-icarus \
            --gene-finding \
            --rna-finding \
            --conserved-genes-finding
        """

def get_all_quast_reports(wildcards):
    checkpoint_output = checkpoints.get_genome_list.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        genomes = [line.strip() for line in f if line.strip()]
    return expand(f"{RESULTS_DIR}/quast/{{genome}}/report.tsv", genome=genomes)

# ============================================================================
# STAGE 11: FINAL REPORTING
# ============================================================================
rule create_final_genome_report:
    input:
        quast_reports = get_all_quast_reports,
        contamination_checks = expand(f"{RESULTS_DIR}/contamination_check/{{sample}}/checkm_table.tsv", sample=FASTQ_SAMPLES),
        processing_summaries = expand(f"{RESULTS_DIR}/final_genomes/{{sample}}/processing_summary.txt", sample=FASTQ_SAMPLES),
        gtdbtk = f"{RESULTS_DIR}/gtdbtk/gtdbtk.bac120.summary.tsv",
        busco_summary = f"{RESULTS_DIR}/busco/busco_summary.tsv"
    output:
        final_report = f"{RESULTS_DIR}/final_report/genome_quality_report.tsv"
    run:
        import pandas as pd
        import os

        # Parse QUAST reports
        quast_data = {}
        for report_file in input.quast_reports:
            genome = report_file.split('/')[-2]
            try:
                df = pd.read_csv(report_file, sep='\t')
                data = {}
                for idx, row in df.iterrows():
                    metric = row.iloc[0].replace('# ', '').replace(' (>= 0 bp)', '').replace(' ', '_')
                    data[metric] = row.iloc[1]
                quast_data[genome] = data
            except Exception as e:
                print(f"Could not read {report_file}: {e}")

        # Parse CheckM data
        checkm_data = {}
        for sample in FASTQ_SAMPLES:
            checkm_file = f"{RESULTS_DIR}/contamination_check/{sample}/checkm_table.tsv"
            if os.path.exists(checkm_file):
                df = pd.read_csv(checkm_file, sep='\t')
                if not df.empty:
                    checkm_data[f"{sample}_genome"] = {
                        'completeness': df.iloc[0]['Completeness'],
                        'contamination': df.iloc[0]['Contamination']
                    }

            # Check for bins
            processing_file = f"{RESULTS_DIR}/final_genomes/{sample}/processing_summary.txt"
            if os.path.exists(processing_file):
                with open(processing_file, 'r') as f:
                    content = f.read()
                    for line in content.split('\n'):
                        if "Selected:" in line:
                            parts = line.split()
                            if len(parts) >= 5:
                                bin_name = parts[1]
                                comp = float(parts[3].replace('(Comp:', '').replace('%,', ''))
                                cont = float(parts[5].replace('Cont:', '').replace('%)', ''))
                                checkm_data[f"{sample}_{bin_name}"] = {
                                    'completeness': comp,
                                    'contamination': cont
                                }

        # Add reference genomes
        for fasta in FASTA_SAMPLES:
            genome_name = fasta.replace('.fa', '')
            checkm_data[genome_name] = {
                'completeness': 100.0,
                'contamination': 0.0
            }

        # Parse GTDB-Tk
        gtdbtk_data = {}
        if os.path.exists(input.gtdbtk):
            df = pd.read_csv(input.gtdbtk, sep='\t')
            for _, row in df.iterrows():
                gtdbtk_data[row['user_genome']] = {
                    'classification': row.get('classification', 'Unknown'),
                    'fastani_ani': row.get('fastani_ani', 'N/A')
                }

        # Parse BUSCO
        busco_data = {}
        if os.path.exists(input.busco_summary):
            df = pd.read_csv(input.busco_summary, sep='\t')
            for _, row in df.iterrows():
                busco_data[row['Genome']] = {
                    'complete_pct': row['Complete_pct'],
                    'single_copy_pct': row['Single_copy_pct'],
                    'duplicated_pct': row['Duplicated_pct']
                }

        # Compile final report
        final_data = []
        for genome, quast in quast_data.items():
            checkm = checkm_data.get(genome, {})
            gtdb = gtdbtk_data.get(genome, {})
            busco = busco_data.get(genome, {})

            completeness = checkm.get('completeness', 0)
            contamination = checkm.get('contamination', 100)

            # Determine quality tier
            quality_tiers = config["quality_tiers"]
            if completeness >= quality_tiers["high"]["completeness_min"] and contamination <= quality_tiers["high"]["contamination_max"]:
                quality_tier = "High"
            elif completeness >= quality_tiers["medium"]["completeness_min"] and contamination <= quality_tiers["medium"]["contamination_max"]:
                quality_tier = "Medium"
            elif completeness >= quality_tiers["low"]["completeness_min"] and contamination <= quality_tiers["low"]["contamination_max"]:
                quality_tier = "Low"
            else:
                quality_tier = "Failed"

            final_data.append({
                'Genome': genome,
                'Quality_Tier': quality_tier,
                'CheckM_Completeness': completeness,
                'CheckM_Contamination': contamination,
                'BUSCO_Complete': busco.get('complete_pct', 'N/A'),
                'BUSCO_Single_Copy': busco.get('single_copy_pct', 'N/A'),
                'BUSCO_Duplicated': busco.get('duplicated_pct', 'N/A'),
                'Size_Mb': round(quast.get('Total_length', 0)/1000000, 2),
                'Contigs': quast.get('contigs', 0),
                'N50': quast.get('N50', 0),
                'GC%': quast.get('GC_(%)', 0),
                'GTDB_Classification': gtdb.get('classification', 'Not classified'),
                'ANI%': gtdb.get('fastani_ani', 'N/A')
            })

        if final_data:
            final_df = pd.DataFrame(final_data)
            final_df = final_df.sort_values(['Quality_Tier', 'CheckM_Completeness'], ascending=[True, False])
            final_df.to_csv(output.final_report, sep='\t', index=False)
        else:
            pd.DataFrame().to_csv(output.final_report, sep='\t', index=False)


