#!/usr/bin/env python3
"""
UpSet Plot for Reactions Across Organisms -> does enrhichment for missing reactions in Marinacidobacterium_sp
Shows overlaps of reactions between organisms for different quality levels
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_contents
import os
import glob
import warnings
from collections import defaultdict
from scipy.stats import fisher_exact, false_discovery_control
from matplotlib.lines import Line2D

warnings.filterwarnings('ignore')

# ==================== FILE PATHS ====================
REACTIONS_DIR = r".\gapseq_pathways"
OUTPUT_DIR = r".\outdir"
PATHWAY_METADATA_FILE = r".\meta_pwy.tbl"  #get from gapseq git repo

# ==================== DATABASES ====================

# Color Palettes
DEEPSEA_PRIMARY = {
    'abyss': '#0c1b33',
    'midnight': '#1e3a5f',
    'deep_blue': '#234e70',
    'ocean': '#2c5f8e',
    'twilight': '#3875a8',
    'marine': '#4a90c2',
    'reef': '#5fa8d3',
    'lagoon': '#7bbee8',
    'shallow': '#9dd1f1',
    'seafoam': '#bde4f4',
    'biolume_green': '#00ffc3',
    'biolume_cyan': '#00e5ff',
}

DEEPSEA_SECONDARY = {
    'pearl': '#f8f9fa',
    'shell': '#e9ecef',
    'fog': '#dee2e6',
    'stone': '#adb5bd',
    'depth': '#6c757d',
    'trench': '#495057',
    'void': '#343a40',
    'ink': '#212529',
    'jellyfish': '#e4c1f9',
    'urchin': '#d0b8e3',
    'starfish': '#fcb1a6',
    'current': '#a9def9',
    'wave': '#d0f4de',
    'plankton': '#fcf6bd',
}

PHYLUM_COLORS = {
    'Pseudomonadota': DEEPSEA_SECONDARY['current'],
    'Acidobacteriota': DEEPSEA_SECONDARY['wave'],
    'Bacteroidota': DEEPSEA_SECONDARY['jellyfish'],
    'Actinomycetota': DEEPSEA_SECONDARY['starfish']
}

SUBSYSTEM_COLORS = {
    'Biosynthesis': DEEPSEA_PRIMARY['marine'],
    'Degradation': DEEPSEA_PRIMARY['deep_blue'],
    'Energy-Metabolism': DEEPSEA_PRIMARY['twilight'],
    'Detoxification': DEEPSEA_PRIMARY['ocean'],
    'Metabolic-Clusters': DEEPSEA_PRIMARY['reef'],
    'Generation-of-Precursor-Metabolites-and-Energy': DEEPSEA_PRIMARY['lagoon'],
    'Superpathways': DEEPSEA_PRIMARY['midnight'],
    'Other': DEEPSEA_SECONDARY['stone'],
    'Unknown': DEEPSEA_SECONDARY['fog']
}

ORGANISM_PHYLUM = {
    'Vreelandella_maris': 'Pseudomonadota',
    'Alloalcanivorax_venustensis': 'Pseudomonadota',
    'Pseudoalteromonas_tetraodonis': 'Pseudomonadota',
    'Cobetia_amphilecti': 'Pseudomonadota',
    'Pseudomonas_E_sp': 'Pseudomonadota',
    'Marinobacter_sp': 'Pseudomonadota',
    'Polaribacter_sp': 'Bacteroidota',
    'Paraglaciecola_sp': 'Pseudomonadota',
    'Pseudovibrio_sp': 'Pseudomonadota',
    'Halopseudomonas_neustonica': 'Pseudomonadota',
    'Salegentibacter_sp': 'Bacteroidota',
    'Leeuwenhoekiella_aequorea': 'Bacteroidota',
    'Mycobacterium_sp': 'Actinomycetota',
    'Gelidibacter_sp': 'Bacteroidota',
    'Methylophaga_sp002696735': 'Pseudomonadota',
    'Alloalcanivorax_sp002389055': 'Pseudomonadota',
    'Psychrobacter_sp': 'Pseudomonadota',
    'Rhodoglobus_sp963974405': 'Actinomycetota',
    'Marinacidobacterium_sp': 'Acidobacteriota'
}

ORGANISM_GENOME_MAPPING = {
    'Vreelandella_maris': '12_genome',
    'Alloalcanivorax_venustensis': '163_9',
    'Pseudoalteromonas_tetraodonis': '1_genome',
    'Cobetia_amphilecti': '291_genome',
    'Pseudomonas_E_sp': '2_genome',
    'Marinobacter_sp': '336_genome',
    'Polaribacter_sp': '3_genome',
    'Paraglaciecola_sp': '459_genome',
    'Pseudovibrio_sp': '468_genome',
    'Halopseudomonas_neustonica': '482_12',
    'Salegentibacter_sp': '482_9',
    'Leeuwenhoekiella_aequorea': '4_genome',
    'Mycobacterium_sp': '536_0',
    'Gelidibacter_sp': '539_genome',
    'Methylophaga_sp002696735': '581_2',
    'Alloalcanivorax_sp002389055': '581_7',
    'Psychrobacter_sp': '5_genome',
    'Rhodoglobus_sp963974405': '6_genome',
    'Marinacidobacterium_sp': 'Marinacidobacteraceae'
}

ORGANISM_DISPLAY = {
    'Vreelandella_maris': 'Vreelandella maris',
    'Alloalcanivorax_venustensis': 'Alloalcanivorax venustensis',
    'Pseudoalteromonas_tetraodonis': 'Pseudoalteromonas tetraodonis',
    'Cobetia_amphilecti': 'Cobetia amphilecti',
    'Pseudomonas_E_sp': 'Pseudomonas sp.',
    'Marinobacter_sp': 'Marinobacter sp.',
    'Polaribacter_sp': 'Polaribacter sp.',
    'Paraglaciecola_sp': 'Paraglaciecola sp.',
    'Pseudovibrio_sp': 'Pseudovibrio sp.',
    'Halopseudomonas_neustonica': 'Halopseudomonas neustonica',
    'Salegentibacter_sp': 'Salegentibacter sp.',
    'Leeuwenhoekiella_aequorea': 'Leeuwenhoekiella aequorea',
    'Mycobacterium_sp': 'Mycobacterium sp.',
    'Gelidibacter_sp': 'Gelidibacter sp.',
    'Methylophaga_sp002696735': 'Methylophaga sp.',
    'Alloalcanivorax_sp002389055': 'Alloalcanivorax sp.',
    'Psychrobacter_sp': 'Psychrobacter sp.',
    'Rhodoglobus_sp963974405': 'Rhodoglobus sp.',
    'Marinacidobacterium_sp': 'Marinacidobacterium sp.'
}

GROWTH_DATA = {
    'Vreelandella_maris': {'Taurine': 0, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 0},
    'Alloalcanivorax_venustensis': {'Taurine': 1, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 0},
    'Pseudoalteromonas_tetraodonis': {'Taurine': 1, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 0},
    'Cobetia_amphilecti': {'Taurine': 0, 'Creatinine': 0, 'Carnitine': 1, 'Xylan': 0, 'Chitin': 0},
    'Pseudomonas_E_sp': {'Taurine': 1, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 0},
    'Marinobacter_sp': {'Taurine': 0, 'Creatinine': 0, 'Carnitine': 1, 'Xylan': 0, 'Chitin': 0},
    'Polaribacter_sp': {'Taurine': 1, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 0},
    'Paraglaciecola_sp': {'Taurine': 0, 'Creatinine': 0, 'Carnitine': 1, 'Xylan': 0, 'Chitin': 0},
    'Pseudovibrio_sp': {'Taurine': 0, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 1, 'Chitin': 0},
    'Halopseudomonas_neustonica': {'Taurine': 1, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 0},
    'Salegentibacter_sp': {'Taurine': 1, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 0},
    'Leeuwenhoekiella_aequorea': {'Taurine': 1, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 0},
    'Mycobacterium_sp': {'Taurine': 0, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 1},
    'Gelidibacter_sp': {'Taurine': 0, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 1},
    'Methylophaga_sp002696735': {'Taurine': 0, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 1, 'Chitin': 0},
    'Alloalcanivorax_sp002389055': {'Taurine': 0, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 1, 'Chitin': 0},
    'Psychrobacter_sp': {'Taurine': 1, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 0},
    'Rhodoglobus_sp963974405': {'Taurine': 1, 'Creatinine': 0, 'Carnitine': 0, 'Xylan': 0, 'Chitin': 0},
    'Marinacidobacterium_sp': {'Taurine': 1, 'Creatinine': 1, 'Carnitine': 1, 'Xylan': 1, 'Chitin': 1}
}

# ==================== STYLE SETTINGS ====================
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.facecolor'] = DEEPSEA_SECONDARY['pearl']
plt.rcParams['figure.facecolor'] = 'white'


# ==================== HELPER FUNCTIONS ====================

def get_organism_color(organism_name):
    """Get color for an organism based on its phylum."""
    organism_key = None
    for key, display_name in ORGANISM_DISPLAY.items():
        if display_name == organism_name:
            organism_key = key
            break

    if organism_key and organism_key in ORGANISM_PHYLUM:
        phylum = ORGANISM_PHYLUM[organism_key]
        return PHYLUM_COLORS.get(phylum, DEEPSEA_SECONDARY['stone'])
    return DEEPSEA_SECONDARY['stone']


def categorize_quality(bitscore, coverage):
    """Categorize reaction quality based on bitscore and coverage"""
    if bitscore == 0:
        return 'none'
    elif bitscore >= 200 and coverage >= 70:
        return 'high'
    elif bitscore >= 100 or (bitscore >= 50 and coverage >= 50):
        return 'moderate'
    else:
        return 'weak'


def categorize_complex(subunit_qualities):
    """Categorize complex completeness based on subunit qualities"""
    if not subunit_qualities:
        return 'incomplete'
    high_count = sum(1 for q in subunit_qualities if q == 'high')
    moderate_count = sum(1 for q in subunit_qualities if q == 'moderate')
    none_count = sum(1 for q in subunit_qualities if q == 'none')
    total = len(subunit_qualities)

    if high_count == total:
        return 'complete'
    elif high_count > 0 and moderate_count > 0 and none_count == 0:
        return 'functional'
    elif moderate_count == total:
        return 'putative'
    else:
        return 'incomplete'


def extract_pathway_subsystems(hierarchy_string):
    """Extract subsystem/hierarchy levels from pathway metadata"""
    if not hierarchy_string or hierarchy_string == 'Unknown':
        return 'Unknown', []

    hierarchy_parts = [part.strip('|').strip() for part in hierarchy_string.split(',') if part.strip()]

    major_categories = {
        'Biosynthesis': ['Biosynthesis'],
        'Degradation': ['Degradation', 'DEGRADATION'],
        'Energy-Metabolism': ['Energy-Metabolism', 'ENERGY-METABOLISM'],
        'Detoxification': ['Detoxification', 'DETOXIFICATION'],
        'Metabolic-Clusters': ['Metabolic-Clusters', 'METABOLIC-CLUSTERS'],
        'Generation-of-Precursor-Metabolites-and-Energy': ['Generation-of-Precursor-Metabolites-and-Energy'],
        'Superpathways': ['Superpathways', 'SUPERPATHWAYS']
    }

    primary_subsystem = 'Other'
    for category, keywords in major_categories.items():
        if any(keyword in hierarchy_parts for keyword in keywords):
            primary_subsystem = category
            break

    return primary_subsystem, hierarchy_parts


# ==================== DATA LOADING FUNCTIONS ====================

def load_pathway_metadata(pathway_file):
    """ Load pathway metadata from gapseq pathway file"""
    pathway_metadata = {}
    try:
        with open(pathway_file, 'r') as f:
            lines = f.readlines()

        if len(lines) == 0:
            return pathway_metadata

        headers = lines[0].strip().split('\t')
        if 'id' not in headers or 'name' not in headers:
            return pathway_metadata

        id_idx = headers.index('id') if 'id' in headers else None
        name_idx = headers.index('name') if 'name' in headers else None
        hierarchy_idx = headers.index('hierarchy') if 'hierarchy' in headers else None
        rea_id_idx = headers.index('reaId') if 'reaId' in headers else None
        rea_nr_idx = headers.index('reaNr') if 'reaNr' in headers else None

        for line_num, line in enumerate(lines[1:], start=2):
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) > max(id_idx or 0, name_idx or 0):
                    pathway_id = parts[id_idx].strip() if id_idx is not None and id_idx < len(parts) else ''
                    pathway_name = parts[name_idx].strip() if name_idx is not None and name_idx < len(parts) else ''
                    hierarchy = parts[hierarchy_idx].strip() if hierarchy_idx is not None and hierarchy_idx < len(
                        parts) else ''
                    reactions = parts[rea_id_idx].strip() if rea_id_idx is not None and rea_id_idx < len(parts) else ''
                    rea_nr = parts[rea_nr_idx].strip() if rea_nr_idx is not None and rea_nr_idx < len(parts) else ''

                    pathway_id = pathway_id.strip('|').strip()
                    pathway_name = pathway_name.strip('|').strip()

                    try:
                        total_reactions = int(rea_nr) if rea_nr and rea_nr != '' else 0
                    except:
                        reaction_list = [r.strip() for r in reactions.split(',') if r.strip()]
                        total_reactions = len(reaction_list)

                    if pathway_id:
                        pathway_metadata[pathway_id] = {
                            'name': pathway_name,
                            'hierarchy': hierarchy,
                            'reactions': reactions.split(',') if reactions else [],
                            'total_reactions': total_reactions
                        }

    except Exception as e:
        pass

    return pathway_metadata


def load_organism_reactions(reactions_dir):
    """Load reactions for each organism and categorize them"""
    reaction_files = glob.glob(os.path.join(reactions_dir, "*-all-Reactions.tbl"))

    genome_to_organism = {v: k for k, v in ORGANISM_GENOME_MAPPING.items()}

    organism_reactions = {}
    reaction_annotations = {}

    for file in reaction_files:
        filename = os.path.basename(file)
        genome_key = filename.replace('-all-Reactions.tbl', '')
        organism = genome_to_organism.get(genome_key, genome_key)
        if organism not in ORGANISM_DISPLAY:
            continue
        organism_display = ORGANISM_DISPLAY[organism]

        try:
            high_quality_reactions = set()
            putative_reactions = set()
            regular_reactions = {}
            complex_reactions = {}

            with open(file, 'r') as f:
                lines = f.readlines()

            header_idx = None
            headers = None
            for i, line in enumerate(lines):
                if line.startswith('#'):
                    continue
                if line.startswith('rxn'):
                    headers = line.strip().split('\t')
                    header_idx = i
                    break

            if not headers:
                continue

            bitscore_idx = headers.index('bitscore') if 'bitscore' in headers else None
            qcovs_idx = headers.index('qcovs') if 'qcovs' in headers else None
            name_idx = headers.index('name') if 'name' in headers else None
            ec_idx = headers.index('ec') if 'ec' in headers else None
            pathway_idx = headers.index('pathway') if 'pathway' in headers else None

            subunit_col_idx = None
            for test_line in lines[header_idx + 1:min(header_idx + 10, len(lines))]:
                if test_line.startswith('#'):
                    continue
                parts = test_line.strip().split('\t')
                for idx in range(len(parts)):
                    if 'Subunit' in parts[idx]:
                        subunit_col_idx = idx
                        break
                if subunit_col_idx:
                    break

            for line in lines[header_idx + 1:]:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                rxn_id = parts[0]

                try:
                    bitscore = float(parts[bitscore_idx]) if bitscore_idx and parts[bitscore_idx] != 'NA' else 0
                    qcovs = float(parts[qcovs_idx]) if qcovs_idx and parts[qcovs_idx] != 'NA' else 0
                except:
                    bitscore = qcovs = 0

                reaction_name = parts[name_idx] if name_idx and name_idx < len(parts) else 'Unknown'
                ec_number = parts[ec_idx] if ec_idx and ec_idx < len(parts) else 'Unknown'
                pathways = parts[pathway_idx] if pathway_idx and pathway_idx < len(parts) else 'Unknown'

                is_complex = False
                subunit_info = None
                if subunit_col_idx and subunit_col_idx < len(parts):
                    subunit_info = parts[subunit_col_idx].strip()
                    if subunit_info and 'Subunit' in subunit_info:
                        is_complex = True

                if is_complex:
                    if rxn_id not in complex_reactions:
                        complex_reactions[rxn_id] = {
                            'subunits': {},
                            'name': reaction_name,
                            'ec': ec_number,
                            'pathways': pathways
                        }
                    complex_reactions[rxn_id]['subunits'][subunit_info] = {
                        'bitscore': bitscore,
                        'qcovs': qcovs,
                        'quality': categorize_quality(bitscore, qcovs)
                    }
                else:
                    if rxn_id not in regular_reactions or bitscore > regular_reactions[rxn_id]['bitscore']:
                        regular_reactions[rxn_id] = {
                            'bitscore': bitscore,
                            'qcovs': qcovs,
                            'quality': categorize_quality(bitscore, qcovs),
                            'name': reaction_name,
                            'ec': ec_number,
                            'pathways': pathways
                        }

            for rxn_id, metrics in regular_reactions.items():
                reaction_annotations[rxn_id] = {
                    'name': metrics['name'],
                    'ec': metrics['ec'],
                    'pathways': metrics['pathways']
                }
                quality = metrics['quality']
                if quality == 'high':
                    high_quality_reactions.add(rxn_id)
                    putative_reactions.add(rxn_id)
                elif quality == 'moderate':
                    putative_reactions.add(rxn_id)

            for complex_id, complex_data in complex_reactions.items():
                reaction_annotations[complex_id] = {
                    'name': complex_data['name'],
                    'ec': complex_data['ec'],
                    'pathways': complex_data['pathways']
                }
                subunit_qualities = [s['quality'] for s in complex_data['subunits'].values()]
                category = categorize_complex(subunit_qualities)
                if category == 'complete':
                    high_quality_reactions.add(complex_id)
                    putative_reactions.add(complex_id)
                elif category == 'functional':
                    putative_reactions.add(complex_id)
                elif category == 'putative':
                    putative_reactions.add(complex_id)

            organism_reactions[organism_display] = {
                'high_quality': high_quality_reactions,
                'high_plus_putative': putative_reactions
            }

        except Exception as e:
            pass

    return organism_reactions, reaction_annotations


# ==================== ANALYSIS FUNCTIONS ====================

def perform_pathway_enrichment(target_reactions, background_reactions, reaction_annotations, pathway_metadata=None):
    """ Perform pathway enrichment analysis using Fisher's exact test"""
    pathway_counts_target = defaultdict(int)
    pathway_counts_background = defaultdict(int)

    for reaction in target_reactions:
        if reaction in reaction_annotations:
            pathways = reaction_annotations[reaction]['pathways']
            if pathways and pathways != 'Unknown' and pathways != 'NA':
                pathway_list = [p.strip() for p in pathways.split('|') if p.strip()]
                for pathway in pathway_list:
                    if pathway:
                        pathway_counts_target[pathway] += 1

    for reaction in background_reactions:
        if reaction in reaction_annotations:
            pathways = reaction_annotations[reaction]['pathways']
            if pathways and pathways != 'Unknown' and pathways != 'NA':
                pathway_list = [p.strip() for p in pathways.split('|') if p.strip()]
                for pathway in pathway_list:
                    if pathway:
                        pathway_counts_background[pathway] += 1

    enrichment_results = []
    for pathway in pathway_counts_background.keys():
        a = pathway_counts_target.get(pathway, 0)
        b = len(target_reactions) - a
        c = pathway_counts_background[pathway] - a
        d = len(background_reactions) - len(target_reactions) - c

        if a == 0 or c == 0:
            continue

        pathway_name = pathway
        total_pathway_reactions = 'Unknown'
        found_in_dataset = pathway_counts_background[pathway]
        hierarchy_info = 'Unknown'
        subsystem = 'Unknown'

        pathway_key = pathway.strip('|')
        if pathway_metadata and pathway_key in pathway_metadata:
            pathway_name = pathway_metadata[pathway_key]['name']
            total_pathway_reactions = pathway_metadata[pathway_key]['total_reactions']
            hierarchy_info = pathway_metadata[pathway_key].get('hierarchy', 'Unknown')
            subsystem, _ = extract_pathway_subsystems(hierarchy_info)

        try:
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')

            target_freq = a / len(target_reactions) if len(target_reactions) > 0 else 0
            background_freq = pathway_counts_background[pathway] / len(background_reactions) if len(
                background_reactions) > 0 else 0
            fold_enrichment = target_freq / background_freq if background_freq > 0 else float('inf')

            enrichment_results.append({
                'Pathway_ID': pathway,
                'Pathway_Name': pathway_name,
                'Subsystem': subsystem,
                'Hierarchy': hierarchy_info,
                'Target_Count': a,
                'Target_Total': len(target_reactions),
                'Background_Count': pathway_counts_background[pathway],
                'Background_Total': len(background_reactions),
                'Total_Pathway_Reactions': total_pathway_reactions,
                'Found_In_Dataset': found_in_dataset,
                'Target_Frequency': target_freq,
                'Background_Frequency': background_freq,
                'Fold_Enrichment': fold_enrichment,
                'P_value': p_value,
                'Odds_Ratio': odds_ratio
            })
        except:
            continue

    if not enrichment_results:
        return pd.DataFrame()

    enrichment_df = pd.DataFrame(enrichment_results)
    enrichment_df = enrichment_df.sort_values('P_value')

    if len(enrichment_df) > 0:
        try:
            enrichment_df['FDR_P_value'] = false_discovery_control(enrichment_df['P_value'])
        except:
            enrichment_df['FDR_P_value'] = enrichment_df['P_value']

    return enrichment_df


def perform_missing_reaction_enrichment(organism_reactions, reaction_annotations, category, output_dir):
    """Perform pathway enrichment analysis on reactions missing from Marinacidobacterium_sp"""

    organisms = list(organism_reactions.keys())
    marina_reactions = set(organism_reactions['Marinacidobacterium sp.'][category])
    other_organisms = [org for org in organisms if org != 'Marinacidobacterium sp.']

    if len(other_organisms) == 0:
        return None

    shared_by_others = set(organism_reactions[other_organisms[0]][category])
    for organism in other_organisms[1:]:
        shared_by_others = shared_by_others.intersection(
            set(organism_reactions[organism][category])
        )

    all_other_reactions = set()
    for organism in other_organisms:
        all_other_reactions.update(organism_reactions[organism][category])

    missing_from_marina = shared_by_others - marina_reactions

    if len(missing_from_marina) == 0:
        return None

    pathway_metadata = load_pathway_metadata(PATHWAY_METADATA_FILE)

    enrichment_df = perform_pathway_enrichment(
        target_reactions=missing_from_marina,
        background_reactions=shared_by_others,
        reaction_annotations=reaction_annotations,
        pathway_metadata=pathway_metadata
    )

    if len(enrichment_df) > 0:
        csv_file = os.path.join(output_dir, f'missing_from_marina_pathway_enrichment_{category}.csv')
        enrichment_df.to_csv(csv_file, index=False)

        create_missing_reaction_dot_plot(enrichment_df, output_dir, category)
        create_pathway_missing_percentage_plot(
            missing_from_marina, shared_by_others, reaction_annotations,
            pathway_metadata, output_dir, category
        )

        missing_reactions_data = []
        for reaction in sorted(missing_from_marina):
            if reaction in reaction_annotations:
                ann = reaction_annotations[reaction]
                missing_reactions_data.append({
                    'Reaction_ID': reaction,
                    'Name': ann['name'],
                    'EC': ann['ec'],
                    'Pathways': ann['pathways']
                })
            else:
                missing_reactions_data.append({
                    'Reaction_ID': reaction,
                    'Name': 'Unknown',
                    'EC': 'Unknown',
                    'Pathways': 'Unknown'
                })

        missing_reactions_df = pd.DataFrame(missing_reactions_data)
        reactions_csv_file = os.path.join(output_dir, f'missing_from_marina_reactions_{category}.csv')
        missing_reactions_df.to_csv(reactions_csv_file, index=False)

        return enrichment_df
    else:
        return None


# ==================== VISUALIZATION FUNCTIONS ====================

def create_missing_reaction_dot_plot(enrichment_df, output_dir, category, top_n=30):
    """Create pathway enrichment dot plot for missing reactions"""
    if len(enrichment_df) == 0:
        return False

    sig_data = enrichment_df[
        (enrichment_df['P_value'] < 0.05) &
        (enrichment_df['Fold_Enrichment'] > 1.2)
        ].head(top_n)

    if len(sig_data) == 0:
        return False

    sig_data = sig_data.sort_values('Fold_Enrichment', ascending=True)

    fig, ax = plt.subplots(figsize=(12, 10))
    fig.patch.set_facecolor('white')
    ax.set_facecolor(DEEPSEA_SECONDARY['pearl'])

    colors = []
    for p_val in sig_data['P_value']:
        if p_val < 0.001:
            colors.append(DEEPSEA_PRIMARY['abyss'])
        elif p_val < 0.01:
            colors.append(DEEPSEA_PRIMARY['midnight'])
        elif p_val < 0.05:
            colors.append(DEEPSEA_PRIMARY['marine'])
        else:
            colors.append(DEEPSEA_PRIMARY['lagoon'])

    y_positions = np.arange(len(sig_data))
    dot_sizes = sig_data['Target_Count'] * 30

    scatter = ax.scatter(
        sig_data['FDR_P_value'] if 'FDR_P_value' in sig_data.columns else sig_data['P_value'],
        y_positions,
        c=colors,
        s=dot_sizes,
        alpha=0.8,
        edgecolors=DEEPSEA_PRIMARY['abyss'],
        linewidth=0.8,
        zorder=5
    )

    for i, (_, row) in enumerate(sig_data.iterrows()):
        label = f"{int(row['Target_Count'])}/{int(row['Background_Count'])}"
        ax.annotate(
            label,
            xy=(row.get('FDR_P_value', row['P_value']), i),
            xytext=(8, 0),
            textcoords='offset points',
            ha='left',
            va='center',
            fontsize=9,
            fontweight='bold',
            color=DEEPSEA_PRIMARY['abyss']
        )

    pathway_names = []
    for name in sig_data['Pathway_Name']:
        if len(name) > 60:
            pathway_names.append(name[:57] + '...')
        else:
            pathway_names.append(name)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(pathway_names, fontsize=10, fontweight='bold')

    ax.set_xscale('log')
    ax.set_xlabel('FDR-corrected P-value', fontsize=12, fontweight='bold', color=DEEPSEA_PRIMARY['abyss'])
    ax.set_ylabel('Metacyc Pathways', fontsize=12, fontweight='bold', color=DEEPSEA_PRIMARY['abyss'])

    title_text = f'Missing Reaction Pathways from Marinacidobacterium sp.\n({category.replace("_", " ").title()})'
    ax.set_title(title_text, fontsize=13, fontweight='bold', pad=20, color=DEEPSEA_PRIMARY['abyss'])

    ax.axvline(x=0.05, color=DEEPSEA_PRIMARY['abyss'], linestyle='--', alpha=0.6, linewidth=2, zorder=1)
    ax.axvline(x=0.01, color=DEEPSEA_SECONDARY['depth'], linestyle=':', alpha=0.5, linewidth=1.5, zorder=1)
    ax.axvline(x=0.001, color=DEEPSEA_SECONDARY['stone'], linestyle=':', alpha=0.5, linewidth=1.5, zorder=1)

    ax.grid(True, alpha=0.3, axis='x', color=DEEPSEA_SECONDARY['stone'])
    ax.set_axisbelow(True)

    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=DEEPSEA_PRIMARY['abyss'],
                   markersize=8, alpha=0.8, markeredgecolor=DEEPSEA_PRIMARY['abyss'],
                   label='P < 0.001', linewidth=0),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=DEEPSEA_PRIMARY['midnight'],
                   markersize=8, alpha=0.8, markeredgecolor=DEEPSEA_PRIMARY['abyss'],
                   label='0.001 ≤ P < 0.01', linewidth=0),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=DEEPSEA_PRIMARY['marine'],
                   markersize=8, alpha=0.8, markeredgecolor=DEEPSEA_PRIMARY['abyss'],
                   label='0.01 ≤ P < 0.05', linewidth=0)
    ]

    legend = ax.legend(handles=legend_elements, title='P-value',
                       loc='lower right', frameon=True, fancybox=True, shadow=True,
                       title_fontsize=11, fontsize=10)
    legend.get_frame().set_facecolor(DEEPSEA_SECONDARY['pearl'])
    legend.get_frame().set_edgecolor(DEEPSEA_SECONDARY['stone'])

    p_values = sig_data.get('FDR_P_value', sig_data['P_value'])
    ax.set_xlim([p_values.max() * 10, p_values.min() / 10])

    plt.tight_layout()

    safe_category = category.replace(' ', '_').replace('+', 'plus')
    output_file = os.path.join(output_dir, f'missing_from_marina_pathway_enrichment_ocean_{safe_category}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    output_file_pdf = os.path.join(output_dir, f'missing_from_marina_pathway_enrichment_ocean_{safe_category}.pdf')
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()

    return True


def create_pathway_missing_percentage_plot(missing_from_marina, shared_by_others, reaction_annotations,
                                           pathway_metadata, output_dir, category, min_reactions=2):
    """Create a plot showing percentage of each pathway that is missing from Marina"""

    pathway_total_counts = defaultdict(int)
    pathway_missing_counts = defaultdict(int)

    for reaction in shared_by_others:
        if reaction in reaction_annotations:
            pathways = reaction_annotations[reaction]['pathways']
            if pathways and pathways != 'Unknown' and pathways != 'NA':
                pathway_list = [p.strip() for p in pathways.split('|') if p.strip()]
                for pathway in pathway_list:
                    if pathway:
                        pathway_total_counts[pathway] += 1

    for reaction in missing_from_marina:
        if reaction in reaction_annotations:
            pathways = reaction_annotations[reaction]['pathways']
            if pathways and pathways != 'Unknown' and pathways != 'NA':
                pathway_list = [p.strip() for p in pathways.split('|') if p.strip()]
                for pathway in pathway_list:
                    if pathway:
                        pathway_missing_counts[pathway] += 1

    pathway_data = []
    for pathway, total_count in pathway_total_counts.items():
        if total_count >= min_reactions:
            missing_count = pathway_missing_counts.get(pathway, 0)
            percent_missing = (missing_count / total_count * 100) if total_count > 0 else 0

            pathway_name = pathway
            subsystem = 'Unknown'
            total_pathway_reactions = 'Unknown'

            pathway_key = pathway.strip('|')
            if pathway_metadata and pathway_key in pathway_metadata:
                pathway_name = pathway_metadata[pathway_key]['name']
                hierarchy_info = pathway_metadata[pathway_key].get('hierarchy', 'Unknown')
                subsystem, _ = extract_pathway_subsystems(hierarchy_info)
                total_pathway_reactions = pathway_metadata[pathway_key]['total_reactions']

            pathway_data.append({
                'Pathway_ID': pathway,
                'Pathway_Name': pathway_name,
                'Subsystem': subsystem,
                'Total_In_Dataset': total_count,
                'Missing_From_Marina': missing_count,
                'Percent_Missing': percent_missing,
                'Total_Pathway_Reactions': total_pathway_reactions
            })

    if not pathway_data:
        return False

    pathway_df = pd.DataFrame(pathway_data)
    pathway_df = pathway_df.sort_values('Percent_Missing', ascending=True)
    pathway_df = pathway_df[pathway_df['Percent_Missing'] > 0]
    pathway_df = pathway_df.tail(20)

    if len(pathway_df) == 0:
        return False

    fig, ax = plt.subplots(figsize=(14, 10))
    fig.patch.set_facecolor('white')

    colors = []
    for subsystem in pathway_df['Subsystem']:
        colors.append(SUBSYSTEM_COLORS.get(subsystem, SUBSYSTEM_COLORS['Unknown']))

    y_positions = np.arange(len(pathway_df))
    bars = ax.barh(y_positions, pathway_df['Percent_Missing'],
                   color=colors, alpha=0.8,
                   edgecolor=DEEPSEA_PRIMARY['abyss'], linewidth=0.5)

    for i, (_, row) in enumerate(pathway_df.iterrows()):
        label = f"{int(row['Missing_From_Marina'])}/{int(row['Total_In_Dataset'])}"
        ax.annotate(label, xy=(row['Percent_Missing'] / 2, i),
                    ha='center', va='center', fontsize=9, fontweight='bold',
                    color='white')

    pathway_names = []
    for name in pathway_df['Pathway_Name']:
        if len(name) > 50:
            pathway_names.append(name[:47] + '...')
        else:
            pathway_names.append(name)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(pathway_names, fontsize=10, fontweight='bold')
    ax.set_xlabel('Percentage of Pathway Missing from Marinacidobacterium sp (%)',
                  fontsize=12, fontweight='bold', color=DEEPSEA_PRIMARY['abyss'])
    ax.set_ylabel('Pathways', fontsize=12, fontweight='bold', color=DEEPSEA_PRIMARY['abyss'])

    title_text = f'Pathway Completeness Analysis\nMissing from Marinacidobacterium sp'
    ax.set_title(title_text, fontsize=13, fontweight='bold', pad=20, color=DEEPSEA_PRIMARY['abyss'])

    ax.grid(True, alpha=0.3, axis='x', color=DEEPSEA_SECONDARY['stone'])
    ax.set_axisbelow(True)
    ax.set_facecolor(DEEPSEA_SECONDARY['pearl'])
    ax.set_xlim(0, 100)

    unique_subsystems = pathway_df['Subsystem'].unique()
    legend_elements = []
    for subsystem in sorted(unique_subsystems):
        legend_elements.append(
            plt.Rectangle((0, 0), 1, 1,
                          facecolor=SUBSYSTEM_COLORS.get(subsystem, SUBSYSTEM_COLORS['Unknown']),
                          alpha=0.8, edgecolor=DEEPSEA_PRIMARY['abyss'],
                          label=subsystem.replace('-', ' '))
        )

    legend = ax.legend(handles=legend_elements, title='Subsystem Category',
                       loc='lower right', frameon=True, fancybox=True, shadow=True,
                       title_fontsize=11, fontsize=10)
    legend.get_frame().set_facecolor(DEEPSEA_SECONDARY['pearl'])
    legend.get_frame().set_edgecolor(DEEPSEA_SECONDARY['stone'])

    plt.tight_layout()

    safe_category = category.replace(' ', '_').replace('+', 'plus')
    output_file = os.path.join(output_dir, f'pathway_missing_percentage_ocean_{safe_category}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    output_file_pdf = os.path.join(output_dir, f'pathway_missing_percentage_ocean_{safe_category}.pdf')
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()

    csv_file = os.path.join(output_dir, f'pathway_missing_percentage_{safe_category}.csv')
    pathway_df.to_csv(csv_file, index=False)

    return True


def create_organism_upset_plot(organism_reactions, reaction_annotations, category, output_dir,
                               min_intersection_size=20):
    """Create UpSet plot matching the knockout"""
    organisms = list(organism_reactions.keys())

    organism_reaction_sets = {}
    for organism, data in organism_reactions.items():
        reactions = data[category]
        if reactions:
            organism_reaction_sets[organism] = reactions

    if not organism_reaction_sets:
        return False

    upset_data = from_contents(organism_reaction_sets)

    fig = plt.figure(figsize=(10, 6))
    fig.patch.set_facecolor('white')

    upset = UpSet(upset_data,
                  subset_size='count',
                  intersection_plot_elements=5,
                  show_counts=True,
                  sort_by='cardinality',
                  min_subset_size=min_intersection_size,
                  facecolor=DEEPSEA_PRIMARY['ocean'])

    upset.plot(fig=fig)

    try:
        for ax in fig.get_axes():
            if hasattr(ax, 'get_yticklabels'):
                labels = ax.get_yticklabels()
                if labels and len(labels) > 0:
                    first_text = labels[0].get_text()
                    if first_text in organism_reaction_sets.keys():
                        for label in labels:
                            organism_name = label.get_text()
                            color = get_organism_color(organism_name)
                            label.set_color(color)
                            label.set_fontweight('bold')
                            label.set_fontsize(10)
    except:
        pass

    title_map = {
        'high_quality': 'High Quality Reactions Shared Between Organisms',
        'high_plus_putative': 'High Quality + Putative Reactions Shared Between Organisms'
    }
    plt.suptitle(f'{title_map.get(category, category)}\nShowing intersections with ≥{min_intersection_size} reactions',
                 fontsize=14, fontweight='bold', color=DEEPSEA_PRIMARY['abyss'])

    all_reactions = set()
    for reactions in organism_reaction_sets.values():
        all_reactions.update(reactions)

    fig.text(0.5, 0.95, f'Total reactions: {len(all_reactions)} | Organisms: {len(organisms)}',
             ha='center', fontsize=11, style='italic', color=DEEPSEA_SECONDARY['depth'])

    legend_elements = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor=PHYLUM_COLORS['Pseudomonadota'],
               markersize=10, label='Pseudomonadota'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=PHYLUM_COLORS['Acidobacteriota'],
               markersize=10, label='Acidobacteriota'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=PHYLUM_COLORS['Bacteroidota'],
               markersize=10, label='Bacteroidota'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=PHYLUM_COLORS['Actinomycetota'],
               markersize=10, label='Actinomycetota')
    ]

    legend = fig.legend(handles=legend_elements, loc='upper right', title='Phylum',
                        frameon=True, fancybox=True, shadow=True)
    legend.get_frame().set_facecolor(DEEPSEA_SECONDARY['pearl'])
    legend.get_frame().set_edgecolor(DEEPSEA_SECONDARY['stone'])

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    category_safe = category.replace(' ', '_').replace('+', 'plus')
    output_file = os.path.join(output_dir, f'organism_upset_{category_safe}_ocean.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    output_file_pdf = os.path.join(output_dir, f'organism_upset_{category_safe}_ocean.pdf')
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()

    return True


# ==================== MAIN FUNCTION ====================
def main():
    """Main execution function"""

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    organism_reactions, reaction_annotations = load_organism_reactions(REACTIONS_DIR)

    if not organism_reactions:
        return False

    # Create UpSet plot
    create_organism_upset_plot(
        organism_reactions,
        reaction_annotations,
        'high_quality',
        OUTPUT_DIR,
        min_intersection_size=20
    )

    # Perform missing reaction enrichment analysis
    perform_missing_reaction_enrichment(
        organism_reactions,
        reaction_annotations,
        'high_quality',
        OUTPUT_DIR
    )

    return True


# ==================== EXECUTION ====================

if __name__ == "__main__":
    main()