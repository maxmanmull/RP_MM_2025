#!/usr/bin/env python3

import pandas as pd
import re
import json
from collections import defaultdict
from pathlib import Path

# ==================== FILE PATHS CONFIGURATION ====================
RESULTS_DIR = Path("./results")

# Input files
QUALITY_FILE = RESULTS_DIR / "qc_summary/mag_quality_master.tsv"
GTDB_FILE = RESULTS_DIR / "gtdbtk/gtdbtk.bac120.summary.tsv"
META_FILE = RESULTS_DIR / "./input/meta_input_pipeline.csv"
TREE_FILE = RESULTS_DIR / "gtdbtk/phylogenetic_tree/tree.nwk"
ANTISMASH_DIR = RESULTS_DIR / "antismash"

# Output directory and files
OUTPUT_DIR = RESULTS_DIR / "./tree_analysis/iTOL"
OUTPUT_TREE = OUTPUT_DIR / "tree.nwk"
OUTPUT_PHYLUM = OUTPUT_DIR / "phylum.txt"
OUTPUT_BGC = OUTPUT_DIR / "bgc_classes.txt"
OUTPUT_MEDIA = OUTPUT_DIR / "media.txt"
OUTPUT_RED = OUTPUT_DIR / "red.txt"
OUTPUT_ANI = OUTPUT_DIR / "ani.txt"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ==================== BGC MAPPING  ====================
BGC_CLASS_MAPPING = {
    # PKS types
    'T1PKS': 'PKS',
    'T3PKS': 'PKS',
    'arylpolyene': 'PKS',

    # NRPS types
    'NRPS': 'NRPS',
    'NRPS-like': 'NRPS',
    'NAPAA': 'NRPS',

    # RiPPs
    'RiPP-like': 'RiPPs',
    'lanthipeptide-class-ii': 'RiPPs',
    'thioamitides': 'RiPPs',
    'azole-containing-RiPP': 'RiPPs',
    'ranthipeptide': 'RiPPs',

    # Terpenes
    'terpene': 'Terpenes',
    'terpene-precursor': 'Terpenes',

    # Siderophores/Metallophores
    'NI-siderophore': 'Siderophores',
    'NRP-metallophore': 'Siderophores',

    # Others
    'betalactone': 'Beta-lactones',
    'ectoine': 'Osmoprotectants',
    'redox-cofactor': 'Cofactors',
    'NAGGN': 'Others',
    'hydrogen-cyanide': 'Others',
    'hserlactone': 'Signaling',
    'acyl_amino_acids': 'Others',
    'cytokinin': 'Phytohormones'
}

# Simple color scheme
BGC_CLASS_COLORS = {
    'PKS': '#1f77b4',
    'NRPS': '#ff7f0e',
    'RiPPs': '#2ca02c',
    'Terpenes': '#d62728',
    'Siderophores': '#9467bd',
    'Beta-lactones': '#8c564b',
    'Osmoprotectants': '#e377c2',
    'Cofactors': '#7f7f7f',
    'Signaling': '#bcbd22',
    'Phytohormones': '#17becf',
    'Others': '#cccccc'
}


# ==================== DATA LOADING ====================
def load_data():
    """Load and merge all required datasets"""
    quality = pd.read_csv(QUALITY_FILE, sep="\t")
    gtdb = pd.read_csv(GTDB_FILE, sep="\t")
    meta = pd.read_csv(META_FILE)

    data = quality.merge(gtdb, left_on="Genome", right_on="user_genome", how="left")

    # Extract taxonomy
    data['phylum'] = data['classification'].str.extract(r'p__([^;]+)')[0].fillna('Unknown')
    data['genus'] = data['classification'].str.extract(r'g__([^;]+)')[0].fillna('')
    data['species'] = data['classification'].str.extract(r's__([^;]+)')[0].fillna('')

    # Create clean names
    data['name'] = data.apply(lambda row: create_name(row), axis=1)
    data['red'] = pd.to_numeric(data['red_value'], errors='coerce')
    data['ani'] = pd.to_numeric(data['closest_genome_ani'], errors='coerce')

    return data, meta


def create_name(row):
    """Create clean name for tree visualization"""
    if row['genus']:
        if row['species']:
            sp = row['species'].replace(f"{row['genus']}_", "").replace(f"{row['genus']} ", "")
            return f"{row['genus']}_{sp}"
        return f"{row['genus']}_sp"
    return row['Genome']


# ==================== BGC PARSING ====================
def parse_bgc_data(data):
    """Parse antiSMASH results - directly using product names"""
    bgc_class_data = {}
    all_bgc_classes = set()

    for _, row in data.iterrows():
        genome = row['Genome']
        name = row['name']
        json_path = ANTISMASH_DIR / genome / f"{genome}.json"

        if json_path.exists():
            try:
                with open(json_path) as f:
                    antismash_data = json.load(f)

                bgc_class_counts = defaultdict(int)

                if 'records' in antismash_data:
                    for record in antismash_data['records']:
                        for feature in record.get('features', []):
                            if feature.get('type') == 'region':
                                products = feature.get('qualifiers', {}).get('product', [])
                                for product in products:
                                    # Use product directly as it comes from antiSMASH
                                    bgc_class = BGC_CLASS_MAPPING.get(product, 'Others')
                                    bgc_class_counts[bgc_class] += 1
                                    all_bgc_classes.add(bgc_class)

                bgc_class_data[name] = bgc_class_counts

                if bgc_class_counts:
                    total = sum(bgc_class_counts.values())
                    print(f"{genome}: {total} BGCs - {dict(bgc_class_counts)}")

            except Exception as e:
                print(f"Error parsing {json_path}: {e}")
                bgc_class_data[name] = {}
        else:
            bgc_class_data[name] = {}

    return bgc_class_data, sorted(all_bgc_classes)


# ==================== OTHER FUNCTIONS ====================
def parse_media_data(data, meta):
    """Parse media growth data"""
    media_cols = ["Taurine", "Carnitine", "Xylan", "Chitin"]
    media_results = {}

    for _, row in data.iterrows():
        genome = row['Genome']
        name = row['name']
        genome_id = re.match(r'^(\d+)', genome)
        genome_id = genome_id.group(1) if genome_id else genome

        media_results[name] = []
        for col in media_cols:
            found = 0
            if col in meta.columns:
                for val in meta[col].dropna().astype(str):
                    val = val.replace('.fa', '').strip()
                    if val in [genome_id, genome]:
                        found = 1
                        break
            media_results[name].append(found)

    return media_results, media_cols


def modify_tree(data):
    """Update tree with clean names"""
    with open(TREE_FILE) as f:
        tree = f.read()

    for _, row in data.iterrows():
        old = row['Genome']
        new = row['name']
        patterns = [f"'{old}'", f'"{old}"', f"({old}:", f",{old}:",
                    f"({old},", f",{old},", f",{old})", f"({old})"]
        for p in patterns:
            if p in tree:
                tree = tree.replace(p, p.replace(old, new))

    with open(OUTPUT_TREE, 'w') as f:
        f.write(tree)


def write_itol_files(data, bgc_data, media_results, media_cols, bgc_classes):
    """Generate all iTOL annotation files"""

    # 1. Phylum colors
    phyla = data['phylum'].unique()
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']
    phylum_colors = {p: colors[i % len(colors)] for i, p in enumerate(phyla)}

    with open(OUTPUT_PHYLUM, 'w') as f:
        f.write("TREE_COLORS\nSEPARATOR TAB\nLEGEND_TITLE\tPhylum\n")
        f.write("LEGEND_SHAPES" + "\t1" * len(phyla) + "\n")
        f.write("LEGEND_COLORS" + "".join(f"\t{phylum_colors[p]}" for p in phyla) + "\n")
        f.write("LEGEND_LABELS" + "".join(f"\t{p}" for p in phyla) + "\n")
        f.write("DATA\n")
        for _, row in data.iterrows():
            f.write(f"{row['name']}\trange\t{phylum_colors[row['phylum']]}\t{row['phylum']}\n")

    # 2. BGC classes
    if bgc_classes:
        with open(OUTPUT_BGC, 'w') as f:
            f.write("DATASET_MULTIBAR\nSEPARATOR TAB\n")
            f.write("DATASET_LABEL\tBGC Classes\nCOLOR\t#000000\nWIDTH\t200\n")
            f.write("LEGEND_TITLE\tBiosynthetic Gene Cluster Classes\n")
            f.write("LEGEND_SHAPES" + "\t1" * len(bgc_classes) + "\n")
            f.write("LEGEND_COLORS" + "".join(f"\t{BGC_CLASS_COLORS.get(c, '#999999')}" for c in bgc_classes) + "\n")
            f.write("LEGEND_LABELS" + "".join(f"\t{c}" for c in bgc_classes) + "\n")
            f.write("FIELD_LABELS" + "".join(f"\t{c}" for c in bgc_classes) + "\n")
            f.write("FIELD_COLORS" + "".join(f"\t{BGC_CLASS_COLORS.get(c, '#999999')}" for c in bgc_classes) + "\n")
            f.write("DATA\n")

            for _, row in data.iterrows():
                name = row['name']
                if name in bgc_data:
                    counts = bgc_data[name]
                    values = [str(counts.get(c, 0)) for c in bgc_classes]
                    if any(int(v) > 0 for v in values):
                        f.write(f"{name}\t" + "\t".join(values) + "\n")

    # 3. Media growth
    with open(OUTPUT_MEDIA, 'w') as f:
        f.write("DATASET_BINARY\nSEPARATOR TAB\n")
        f.write("DATASET_LABEL\tMedia Growth\nCOLOR\t#000000\n")
        f.write("FIELD_LABELS" + "".join(f"\t{col}" for col in media_cols) + "\n")
        f.write("FIELD_COLORS\t#2166ac\t#762a83\t#5aae61\t#de77ae\n")
        f.write("FIELD_SHAPES" + "\t1" * len(media_cols) + "\n")
        f.write("DATA\n")
        for name, values in media_results.items():
            f.write(name + "\t" + "\t".join(str(v) for v in values) + "\n")

    # 4. RED values
    with open(OUTPUT_RED, 'w') as f:
        f.write("DATASET_TEXT\nSEPARATOR TAB\n")
        f.write("DATASET_LABEL\tRED\nCOLOR\t#E74C3C\nMARGIN\t20\nDATA\n")
        for _, row in data.iterrows():
            val = f"{row['red']:.2f}" if pd.notna(row['red']) else "-"
            f.write(f"{row['name']}\t{val}\t1\t#E74C3C\tbold\t1\t0\n")

    # 5. ANI values
    with open(OUTPUT_ANI, 'w') as f:
        f.write("DATASET_TEXT\nSEPARATOR TAB\n")
        f.write("DATASET_LABEL\tANI %\nCOLOR\t#000000\nMARGIN\t60\nDATA\n")
        for _, row in data.iterrows():
            val = f"{row['ani']:.1f}" if pd.notna(row['ani']) else "-"
            f.write(f"{row['name']}\t{val}\t1\t#000000\tbold\t1\t0\n")


# ==================== MAIN ====================
def main():
    data, meta = load_data()
    bgc_class_data, bgc_classes = parse_bgc_data(data)
    media_results, media_cols = parse_media_data(data, meta)
    modify_tree(data)
    write_itol_files(data, bgc_class_data, media_results, media_cols, bgc_classes)



    if bgc_classes:
        print(f"\nBGC Classes found:")
        for bgc_class in bgc_classes:
            total = sum(d.get(bgc_class, 0) for d in bgc_class_data.values())
            if total > 0:
                print(f"  {bgc_class}: {total}")



if __name__ == "__main__":
    main()