#!/usr/bin/env python3
"""
Visualizes metabolic pathways and outputs reaction data as CSV files
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
import os
import glob
import warnings

warnings.filterwarnings('ignore')

# ==================== FILE PATHS ====================
REACTIONS_DIR = r".\09_gapseq_pathways"
TRANSPORTERS_DIR = r".\gapseq_transporters"
PATHWAYS_FILE = r".\meta_pwy.tbl"
OUTPUT_DIR = r".\pathway_matrices"

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
    'Pseudomonadota': '#a9def9',
    'Acidobacteriota': '#d0f4de',
    'Bacteroidota': '#e4c1f9',
    'Actinomycetota': '#fcb1a6'
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

PATHWAY_GROUPS = {
    'Taurine': ['PWY-5982', 'PWY-1263', 'PWY-1264', 'TAURINEDEG-PWY', 'PWY0-981', 'PWY-1281', 'PWY-6718'],
    'Xylan': ['PWY-6789', 'PWY-6717', 'XYLCAT-PWY', 'PWY-5516', 'PWY-6760', 'PWY-7294', 'PWY-8020',
              'PWY-8330', 'LARABITOLUTIL-PWY'],
    'Chitin': ['PWY-7118', 'PWY-6855', 'PWY-6902', 'PWY-6906', 'PWY-7822', 'PWY-6517', 'GLUAMCAT-PWY'],
    'Creatinine': ['CRNFORCAT-PWY', 'PWY-4722', 'PWY-4741'],
    'Carnitine': ['CARNMET-PWY', 'PWY-3641', 'PWY-3602', 'PWY-8307', 'PWY-3661']
}

SPECIFIC_REACTIONS = {
    'Chitin': ['3.2.1.132-RXN', '3.2.1.165-RXN']
}

COMPOUND_TRANSPORTERS = {
    'Taurine': ['taurine'],
    'Xylan': ['xylose', 'xylan'],
    'Chitin': ['chitin', 'chitobiose', 'n-acetylglucosamine', 'glcnac', 'glucosamine'],
    'Creatinine': ['creatinine', 'creatine'],
    'Carnitine': ['carnitine']
}

# ==================== STYLE SETTINGS ====================
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.size'] = 8


# ==================== DATA LOADING FUNCTIONS ====================
def load_genome_reactions_with_complex_tracking(reactions_dir):
    """Load genome reaction data tracking all subunits for complexes"""
    genome_reactions = {}
    reaction_files = glob.glob(os.path.join(reactions_dir, "*-all-Reactions.tbl"))

    for file in reaction_files:
        filename = os.path.basename(file)
        genome_key = filename.replace('-all-Reactions.tbl', '')
        try:
            reactions = {}
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

            if headers:
                bitscore_idx = headers.index('bitscore') if 'bitscore' in headers else None
                evalue_idx = headers.index('evalue') if 'evalue' in headers else None
                pident_idx = headers.index('pident') if 'pident' in headers else None
                qcovs_idx = headers.index('qcovs') if 'qcovs' in headers else None

                subunit_col_idx = None
                for test_line in lines[header_idx + 1:min(header_idx + 10, len(lines))]:
                    if test_line.startswith('#'):
                        continue
                    parts = test_line.strip().split('\t')
                    for idx in range(len(parts)):
                        if 'Subunit' in parts[idx] or (idx > 14 and parts[idx] in ['NA', '0', '1']):
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
                    metrics = {}

                    try:
                        if bitscore_idx and bitscore_idx < len(parts):
                            value = parts[bitscore_idx]
                            metrics['bitscore'] = float(value) if value != 'NA' else 0
                        else:
                            metrics['bitscore'] = 0

                        if evalue_idx and evalue_idx < len(parts):
                            value = parts[evalue_idx]
                            metrics['evalue'] = float(value) if value != 'NA' else 1.0
                        else:
                            metrics['evalue'] = 1.0

                        if pident_idx and pident_idx < len(parts):
                            value = parts[pident_idx]
                            metrics['pident'] = float(value) if value != 'NA' else 0
                        else:
                            metrics['pident'] = 0

                        if qcovs_idx and qcovs_idx < len(parts):
                            value = parts[qcovs_idx]
                            metrics['qcovs'] = float(value) if value != 'NA' else 0
                        else:
                            metrics['qcovs'] = 0
                    except:
                        continue

                    subunit_info = 'NA'
                    if subunit_col_idx and subunit_col_idx < len(parts):
                        subunit_info = parts[subunit_col_idx].strip()

                    if rxn_id not in reactions:
                        reactions[rxn_id] = {
                            'subunits': {},
                            'is_complex': False
                        }

                    if subunit_info and subunit_info != 'NA' and 'Subunit' in subunit_info:
                        reactions[rxn_id]['is_complex'] = True
                        reactions[rxn_id]['subunits'][subunit_info] = metrics
                    else:
                        if 'single' not in reactions[rxn_id]['subunits'] or metrics['bitscore'] > reactions[rxn_id][
                            'subunits'].get('single', {}).get('bitscore', 0):
                            reactions[rxn_id]['subunits']['single'] = metrics

            if reactions:
                genome_reactions[genome_key] = reactions

        except Exception as e:
            pass

    return genome_reactions


def load_genome_transporters(transporters_dir):
    """Load genome transporter data with metrics"""
    genome_transporters = {}
    transporter_files = glob.glob(os.path.join(transporters_dir, "*-Transporter.tbl"))

    for file in transporter_files:
        filename = os.path.basename(file)
        genome_key = filename.replace('-Transporter.tbl', '')

        try:
            transporters = {}
            with open(file, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()

            header_idx = None
            headers = None
            for i, line in enumerate(lines):
                if line.startswith('#'):
                    continue
                if line.startswith('id\t'):
                    headers = line.strip().split('\t')
                    header_idx = i
                    break

            if headers:
                id_idx = headers.index('id') if 'id' in headers else 0
                tc_idx = headers.index('tc') if 'tc' in headers else 1
                sub_idx = headers.index('sub') if 'sub' in headers else 2
                rea_idx = headers.index('rea') if 'rea' in headers else 4
                bitscore_idx = headers.index('bitscore') if 'bitscore' in headers else None
                qcovs_idx = headers.index('qcovs') if 'qcovs' in headers else None
                pident_idx = headers.index('pident') if 'pident' in headers else None

                for line in lines[header_idx + 1:]:
                    if line.startswith('#') or not line.strip():
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) <= max([id_idx, tc_idx, sub_idx, rea_idx]):
                        continue

                    transporter_id = parts[id_idx]
                    tc_class = parts[tc_idx] if tc_idx < len(parts) else ''
                    substrate = parts[sub_idx].lower() if sub_idx < len(parts) else ''
                    reactions = parts[rea_idx].split(',') if rea_idx < len(parts) else []

                    metrics = {
                        'substrate': substrate,
                        'tc_class': tc_class,
                        'reactions': reactions,
                        'bitscore': float(parts[bitscore_idx]) if bitscore_idx and bitscore_idx < len(parts) and parts[
                            bitscore_idx] != 'NA' else 0,
                        'qcovs': float(parts[qcovs_idx]) if qcovs_idx and qcovs_idx < len(parts) and parts[
                            qcovs_idx] != 'NA' else 0,
                        'pident': float(parts[pident_idx]) if pident_idx and pident_idx < len(parts) and parts[
                            pident_idx] != 'NA' else 0
                    }

                    simple_id = f"TC_{tc_class}" if tc_class else transporter_id.split('|')[-1]

                    if simple_id in transporters:
                        if metrics['bitscore'] > transporters[simple_id]['bitscore']:
                            transporters[simple_id] = metrics
                    else:
                        transporters[simple_id] = metrics

            if transporters:
                genome_transporters[genome_key] = transporters

        except Exception as e:
            pass

    return genome_transporters


# ==================== HELPER FUNCTIONS ====================
def filter_transporters_for_compound(transporters, compound):
    """Filter transporters relevant to a specific compound based on substrate"""
    if compound not in COMPOUND_TRANSPORTERS:
        return {}

    relevant_keywords = COMPOUND_TRANSPORTERS[compound]
    filtered = {}

    for transporter_id, metrics in transporters.items():
        substrate = metrics.get('substrate', '').lower()
        for keyword in relevant_keywords:
            if keyword.lower() in substrate:
                filtered[transporter_id] = metrics
                break

    return filtered


def get_best_transporter_per_organism(all_genome_transporters, compound, compound_organisms, min_bitscore=50):
    """Get only the best transporter (highest bitscore) for each organism for a given compound"""
    best_transporters = {}
    organism_best_transporter = {}

    for organism in compound_organisms:
        genome_key = ORGANISM_GENOME_MAPPING.get(organism, '')
        if genome_key in all_genome_transporters:
            filtered_transporters = filter_transporters_for_compound(
                all_genome_transporters[genome_key], compound
            )

            best_trans_id = None
            best_score = min_bitscore

            for trans_id, metrics in filtered_transporters.items():
                if metrics['bitscore'] > best_score:
                    best_score = metrics['bitscore']
                    best_trans_id = trans_id

            if best_trans_id:
                organism_best_transporter[organism] = best_trans_id
                if best_trans_id not in best_transporters:
                    best_transporters[best_trans_id] = filtered_transporters[best_trans_id]

    return best_transporters, organism_best_transporter


def get_ocean_themed_color(bitscore):
    """Get color based on bitscore threshold"""
    if bitscore == 0:
        return DEEPSEA_SECONDARY['pearl']
    elif bitscore < 100:
        intensity = bitscore / 100
        r = 0.5 + 0.3 * intensity
        g = 0.1 + 0.2 * intensity
        b = 0.3 + 0.2 * intensity
        return (r, g, b)
    elif bitscore < 200:
        intensity = (bitscore - 100) / 100
        r = 0.8 - 0.2 * intensity
        g = 0.7 + 0.2 * intensity
        b = 0.3 + 0.4 * intensity
        return (r, g, b)
    else:
        intensity = min((bitscore - 200) / 300, 1.0)
        r = 0.1 + 0.1 * intensity
        g = 0.6 + 0.3 * intensity
        b = 0.4 + 0.3 * intensity
        return (r, g, b)


def draw_complex_cell(ax, x, y, subunit_data, cell_size=0.8):
    """Draw a cell split by subunits for multi-subunit complexes"""
    subunit_list = list(subunit_data.keys())
    n_subunits = len(subunit_list)

    if n_subunits == 0:
        return

    if n_subunits == 1:
        metrics = subunit_data[subunit_list[0]]
        bitscore = metrics.get('bitscore', 0)
        coverage = metrics.get('qcovs', 0)

        if bitscore > 0:
            color = get_ocean_themed_color(bitscore)
            radius = cell_size * 0.4 * (0.5 + 0.5 * min(coverage / 100, 1))
            circle = patches.Circle((x, y), radius=radius,
                                    facecolor=color, edgecolor='black', linewidth=0.5)
            ax.add_patch(circle)

    elif n_subunits == 2:
        for i, subunit in enumerate(subunit_list):
            metrics = subunit_data[subunit]
            bitscore = metrics.get('bitscore', 0)

            if bitscore > 0:
                color = get_ocean_themed_color(bitscore)
            else:
                color = DEEPSEA_SECONDARY['fog']

            if i == 0:
                wedge = patches.Wedge((x, y), cell_size * 0.4, 90, 270,
                                      facecolor=color, edgecolor='black', linewidth=0.5)
            else:
                wedge = patches.Wedge((x, y), cell_size * 0.4, -90, 90,
                                      facecolor=color, edgecolor='black', linewidth=0.5)
            ax.add_patch(wedge)

    elif n_subunits == 3:
        for i, subunit in enumerate(subunit_list):
            metrics = subunit_data[subunit]
            bitscore = metrics.get('bitscore', 0)

            if bitscore > 0:
                color = get_ocean_themed_color(bitscore)
            else:
                color = DEEPSEA_SECONDARY['fog']

            start_angle = -90 + i * 120
            wedge = patches.Wedge((x, y), cell_size * 0.4, start_angle, start_angle + 120,
                                  facecolor=color, edgecolor='black', linewidth=0.5)
            ax.add_patch(wedge)

    elif n_subunits == 4:
        for i, subunit in enumerate(subunit_list):
            metrics = subunit_data[subunit]
            bitscore = metrics.get('bitscore', 0)

            if bitscore > 0:
                color = get_ocean_themed_color(bitscore)
            else:
                color = DEEPSEA_SECONDARY['fog']

            start_angle = -45 + i * 90
            wedge = patches.Wedge((x, y), cell_size * 0.4, start_angle, start_angle + 90,
                                  facecolor=color, edgecolor='black', linewidth=0.5)
            ax.add_patch(wedge)

    else:
        angle_per_subunit = 360 / n_subunits
        for i, subunit in enumerate(subunit_list):
            metrics = subunit_data[subunit]
            bitscore = metrics.get('bitscore', 0)

            if bitscore > 0:
                color = get_ocean_themed_color(bitscore)
            else:
                color = DEEPSEA_SECONDARY['fog']

            start_angle = -90 + i * angle_per_subunit
            wedge = patches.Wedge((x, y), cell_size * 0.4, start_angle, start_angle + angle_per_subunit,
                                  facecolor=color, edgecolor='black', linewidth=0.5)
            ax.add_patch(wedge)

    missing_subunits = sum(1 for metrics in subunit_data.values() if metrics.get('bitscore', 0) == 0)
    if missing_subunits > 0 and n_subunits > 1:
        ax.text(x, y, '!', fontsize=12, ha='center', va='center',
                color='red', fontweight='bold', zorder=10)


def create_reaction_dataframes(compound, all_reactions, compound_organisms, all_genome_reactions,
                               all_genome_transporters, organism_best_transporter, pathway_data,
                               reaction_to_pathway, output_dir):
    """Create and save DataFrames with reaction data for each compound"""

    bitscore_data = {}
    detailed_data = []
    metadata_data = []

    for organism in compound_organisms:
        genome_key = ORGANISM_GENOME_MAPPING.get(organism, '')
        organism_name = organism.replace('_', ' ')
        phylum = ORGANISM_PHYLUM.get(organism, 'Unknown')

        organism_row = {
            'Organism': organism_name,
            'Phylum': phylum,
            'Growth_on_compound': GROWTH_DATA[organism].get(compound, 0)
        }

        for reaction in all_reactions:
            bitscore = 0
            evalue = 1.0
            pident = 0
            qcovs = 0
            is_complex = False
            n_subunits = 0
            n_complete_subunits = 0

            if reaction.startswith('TC_'):
                if organism in organism_best_transporter and organism_best_transporter[organism] == reaction:
                    if genome_key in all_genome_transporters:
                        filtered_transporters = filter_transporters_for_compound(
                            all_genome_transporters[genome_key], compound
                        )
                        if reaction in filtered_transporters:
                            metrics = filtered_transporters[reaction]
                            bitscore = metrics.get('bitscore', 0)
                            evalue = metrics.get('evalue', 1.0)
                            pident = metrics.get('pident', 0)
                            qcovs = metrics.get('qcovs', 0)
            else:
                if genome_key in all_genome_reactions and reaction in all_genome_reactions[genome_key]:
                    rxn_data = all_genome_reactions[genome_key][reaction]
                    is_complex = rxn_data.get('is_complex', False)

                    max_bitscore = 0
                    for subunit_metrics in rxn_data['subunits'].values():
                        if subunit_metrics.get('bitscore', 0) > max_bitscore:
                            max_bitscore = subunit_metrics.get('bitscore', 0)
                            evalue = subunit_metrics.get('evalue', 1.0)
                            pident = subunit_metrics.get('pident', 0)
                            qcovs = subunit_metrics.get('qcovs', 0)

                    bitscore = max_bitscore
                    n_subunits = len(rxn_data['subunits'])
                    n_complete_subunits = sum(1 for m in rxn_data['subunits'].values()
                                              if m.get('bitscore', 0) > 0)

            organism_row[reaction] = bitscore

            detailed_record = {
                'Organism': organism_name,
                'Phylum': phylum,
                'Reaction': reaction,
                'BitScore': bitscore,
                'E_value': evalue,
                'Percent_Identity': pident,
                'Query_Coverage': qcovs,
                'Is_Complex': is_complex,
                'N_Subunits': n_subunits,
                'N_Complete_Subunits': n_complete_subunits,
                'Complex_Completeness': f"{n_complete_subunits}/{n_subunits}" if n_subunits > 0 else "NA",
                'Reaction_Present': 1 if bitscore > 0 else 0
            }
            detailed_data.append(detailed_record)

        bitscore_data[organism_name] = organism_row

    for reaction in all_reactions:
        pathways = reaction_to_pathway.get(reaction, [])
        pathway_names = []
        for pwy_id in pathways:
            if pwy_id in pathway_data:
                pathway_names.append(pathway_data[pwy_id]['name'])
            elif pwy_id == 'SPECIFIC':
                pathway_names.append('Compound-Specific Reaction')
            elif pwy_id == 'TRANSPORTER':
                pathway_names.append('Transporter')
            else:
                pathway_names.append(pwy_id)

        metadata_record = {
            'Reaction': reaction,
            'Pathway_IDs': ', '.join(pathways),
            'Pathway_Names': ', '.join(pathway_names),
            'Reaction_Type': 'Transporter' if reaction.startswith('TC_') else 'Enzymatic'
        }
        metadata_data.append(metadata_record)

    bitscore_df = pd.DataFrame.from_dict(bitscore_data, orient='index')
    detailed_df = pd.DataFrame(detailed_data)
    metadata_df = pd.DataFrame(metadata_data)

    summary_data = []
    for organism in compound_organisms:
        organism_name = organism.replace('_', ' ')
        phylum = ORGANISM_PHYLUM.get(organism, 'Unknown')

        org_data = bitscore_df.loc[organism_name]

        n_strong = sum(1 for r in all_reactions if org_data.get(r, 0) >= 200)
        n_moderate = sum(1 for r in all_reactions if 100 <= org_data.get(r, 0) < 200)
        n_weak = sum(1 for r in all_reactions if 0 < org_data.get(r, 0) < 100)
        n_absent = sum(1 for r in all_reactions if org_data.get(r, 0) == 0)

        summary_record = {
            'Organism': organism_name,
            'Phylum': phylum,
            'Growth_on_Compound': GROWTH_DATA[organism].get(compound, 0),
            'Total_Reactions': len(all_reactions),
            'Reactions_Present': n_strong + n_moderate + n_weak,
            'Strong_Reactions': n_strong,
            'Moderate_Reactions': n_moderate,
            'Weak_Reactions': n_weak,
            'Absent_Reactions': n_absent,
            'Percent_Coverage': round(100 * (n_strong + n_moderate + n_weak) / len(all_reactions), 2)
        }
        summary_data.append(summary_record)

    summary_df = pd.DataFrame(summary_data)

    create_latex_tables(compound, bitscore_df, detailed_df, metadata_df, summary_df, output_dir)

    summary_df = summary_df.sort_values(['Growth_on_Compound', 'Percent_Coverage'], ascending=[False, False])
    bitscore_df = bitscore_df.sort_values('Growth_on_compound', ascending=False)

    bitscore_file = os.path.join(output_dir, f"{compound.lower()}_reaction_bitscores.csv")
    bitscore_df.to_csv(bitscore_file)

    detailed_file = os.path.join(output_dir, f"{compound.lower()}_reaction_detailed.csv")
    detailed_df.to_csv(detailed_file, index=False)

    metadata_file = os.path.join(output_dir, f"{compound.lower()}_reaction_metadata.csv")
    metadata_df.to_csv(metadata_file, index=False)

    summary_file = os.path.join(output_dir, f"{compound.lower()}_organism_summary.csv")
    summary_df.to_csv(summary_file, index=False)

    excel_file = os.path.join(output_dir, f"{compound.lower()}_complete_analysis.xlsx")
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
        bitscore_df.to_excel(writer, sheet_name='BitScores')
        detailed_df.to_excel(writer, sheet_name='Detailed', index=False)
        metadata_df.to_excel(writer, sheet_name='Reaction_Metadata', index=False)

    return bitscore_df, detailed_df, metadata_df, summary_df


def create_latex_tables(compound, bitscore_df, detailed_df, metadata_df, summary_df, output_dir):
    """Generate simple LaTeX table"""

    def format_org_name(name):
        parts = name.split('_')
        if len(parts) >= 2:
            return f"{parts[0]}\\\\{' '.join(parts[1:])}"
        return name.replace('_', ' ')

    latex_file = os.path.join(output_dir, f"{compound.lower()}_simple_table.tex")

    with open(latex_file, 'w') as f:
        f.write(f"% Simple {compound} reaction data table\n")
        f.write("\\begin{table}[h]\n")
        f.write("\\centering\n")
        f.write("\\caption{" + compound + " reaction data}\n")
        f.write("\\begin{tabular}{llrrrr}\n")
        f.write("\\hline\n")
        f.write("Organism & Reaction & BitScore & Coverage & Identity & E-value\\\\\n")
        f.write("\\hline\n")

        reaction_cols = [c for c in bitscore_df.columns if c not in ['Phylum', 'Growth_on_compound']]

        for organism in bitscore_df.index:
            org_parts = organism.split()
            if len(org_parts) >= 2:
                org_name = f"{org_parts[0]}\\\\{' '.join(org_parts[1:])}"
            else:
                org_name = organism

            reactions_with_hits = []
            for rxn in reaction_cols:
                org_rxn_data = detailed_df[(detailed_df['Organism'] == organism) &
                                           (detailed_df['Reaction'] == rxn)]

                if not org_rxn_data.empty:
                    row_data = org_rxn_data.iloc[0]
                    if row_data['BitScore'] > 0:
                        reactions_with_hits.append({
                            'reaction': rxn,
                            'bitscore': row_data['BitScore'],
                            'coverage': row_data['Query_Coverage'],
                            'identity': row_data['Percent_Identity'],
                            'evalue': row_data['E_value']
                        })

            reactions_with_hits.sort(key=lambda x: x['bitscore'], reverse=True)

            for i, rxn_data in enumerate(reactions_with_hits):
                if i == 0:
                    f.write(org_name)
                else:
                    f.write("")

                if rxn_data['evalue'] == 0:
                    evalue_str = "0"
                else:
                    evalue_str = f"{rxn_data['evalue']:.2e}"

                f.write(f" & {rxn_data['reaction']} & {rxn_data['bitscore']:.0f} & ")
                f.write(f"{rxn_data['coverage']:.1f} & {rxn_data['identity']:.1f} & ")
                f.write(f"{evalue_str}\\\\\n")

            if reactions_with_hits:
                f.write("\\hline\n")

        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

    return latex_file, None


def create_tiled_heatmap_with_pathways(compound, pathways, all_genome_reactions, all_genome_transporters,
                                       pathway_data, output_dir):
    """Create tiled heatmap"""

    compound_organisms = [org for org, growth in GROWTH_DATA.items() if growth.get(compound, 0) == 1]

    if not compound_organisms:
        return

    pathway_reactions_map = {}
    reaction_to_pathway = {}
    for pwy in pathways:
        if pwy in pathway_data:
            pathway_reactions_map[pwy] = pathway_data[pwy]['reactions']
            for reaction in pathway_data[pwy]['reactions']:
                if reaction not in reaction_to_pathway:
                    reaction_to_pathway[reaction] = []
                reaction_to_pathway[reaction].append(pwy)

    all_reactions = []
    pathway_boundaries = []

    for pwy in pathways:
        if pwy in pathway_reactions_map:
            pwy_reactions = pathway_reactions_map[pwy]
            for rxn in pwy_reactions:
                if rxn not in all_reactions:
                    all_reactions.append(rxn)
            if pwy_reactions:
                current_pos = len(all_reactions) - 1
                pathway_boundaries.append((current_pos, pwy))

    if compound in SPECIFIC_REACTIONS:
        specific_rxns = SPECIFIC_REACTIONS[compound]

        if all_reactions:
            pathway_boundaries.append((len(all_reactions) - 1, 'PATHWAY'))

        for rxn in specific_rxns:
            if rxn not in all_reactions:
                all_reactions.append(rxn)
                if rxn not in reaction_to_pathway:
                    reaction_to_pathway[rxn] = []
                reaction_to_pathway[rxn].append('SPECIFIC')

        if specific_rxns:
            pathway_boundaries.append((len(all_reactions) - 1, 'SPECIFIC'))

    best_transporters, organism_best_transporter = get_best_transporter_per_organism(
        all_genome_transporters, compound, compound_organisms
    )

    compound_transporter_ids = list(best_transporters.keys())

    if compound_transporter_ids:
        if all_reactions:
            pathway_boundaries.append((len(all_reactions) - 1, 'REACTIONS'))

        for transporter_id in sorted(compound_transporter_ids):
            all_reactions.append(transporter_id)
            reaction_to_pathway[transporter_id] = ['TRANSPORTER']

        pathway_boundaries.append((len(all_reactions) - 1, 'TRANSPORTER'))

    if not all_reactions:
        return

    create_reaction_dataframes(
        compound, all_reactions, compound_organisms, all_genome_reactions,
        all_genome_transporters, organism_best_transporter, pathway_data,
        reaction_to_pathway, output_dir
    )

    n_organisms = len(compound_organisms)
    n_reactions = len(all_reactions)

    fig_width = max(14, n_reactions * 0.8)
    fig_height = max(10, n_organisms * 0.6 + 2)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 12], hspace=0.02)

    ax_pathway = fig.add_subplot(gs[0])
    ax_pathway.set_facecolor('white')

    ax = fig.add_subplot(gs[1])
    ax.set_facecolor(DEEPSEA_SECONDARY['pearl'])

    tile_color1 = DEEPSEA_SECONDARY['pearl']
    tile_color2 = '#f0f1f2'

    for i in range(n_organisms):
        for j in range(n_reactions):
            if (i + j) % 2 == 0:
                rect = patches.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                         facecolor=tile_color2,
                                         edgecolor='none', alpha=0.3, zorder=0)
                ax.add_patch(rect)

    pathway_colors = [
        DEEPSEA_PRIMARY['abyss'],
        DEEPSEA_PRIMARY['ocean'],
        DEEPSEA_PRIMARY['twilight'],
        DEEPSEA_PRIMARY['marine'],
        DEEPSEA_PRIMARY['reef'],
        DEEPSEA_PRIMARY['lagoon'],
        DEEPSEA_SECONDARY['urchin'],
        DEEPSEA_SECONDARY['current'],
        DEEPSEA_SECONDARY['starfish'],
        DEEPSEA_SECONDARY['jellyfish']
    ]

    pathway_color_map = {pwy: pathway_colors[i % len(pathway_colors)]
                         for i, pwy in enumerate(pathways)}
    pathway_color_map['SPECIFIC'] = '#ff6b6b'
    pathway_color_map['TRANSPORTER'] = '#00d9ff'

    bar_height = 1.0
    for i, reaction in enumerate(all_reactions):
        if reaction in reaction_to_pathway:
            pathways_for_reaction = reaction_to_pathway[reaction]
            n_pathways = len(pathways_for_reaction)
            segment_height = bar_height / n_pathways

            for j, pwy in enumerate(pathways_for_reaction):
                color = pathway_color_map.get(pwy, '#666666')
                y_pos = j * segment_height
                ax_pathway.bar(i, segment_height, bottom=y_pos, color=color,
                               edgecolor='none', linewidth=0, alpha=0.8)

    ax_pathway.set_xlim(-0.5, n_reactions - 0.5)
    ax_pathway.set_ylim(0, bar_height)
    ax_pathway.set_xticks([])
    ax_pathway.set_yticks([])
    ax_pathway.set_title('PATHWAY ASSIGNMENT', fontsize=12, fontweight='bold',
                         color=DEEPSEA_PRIMARY['midnight'], pad=5)

    for spine in ax_pathway.spines.values():
        spine.set_visible(False)

    for boundary_pos, pwy in pathway_boundaries[:-1]:
        ax.axvline(x=boundary_pos + 0.5, color=DEEPSEA_SECONDARY['stone'],
                   linestyle='--', linewidth=1, alpha=0.5)
        ax_pathway.axvline(x=boundary_pos + 0.5, color=DEEPSEA_SECONDARY['stone'],
                           linestyle='--', linewidth=1, alpha=0.5)

    for i, organism in enumerate(compound_organisms):
        genome_key = ORGANISM_GENOME_MAPPING.get(organism, '')

        for j, reaction in enumerate(all_reactions):
            if reaction.startswith('TC_'):
                if organism in organism_best_transporter and organism_best_transporter[organism] == reaction:
                    if genome_key in all_genome_transporters:
                        filtered_transporters = filter_transporters_for_compound(
                            all_genome_transporters[genome_key], compound
                        )
                        if reaction in filtered_transporters:
                            metrics = filtered_transporters[reaction]
                            if metrics.get('bitscore', 0) > 0:
                                draw_complex_cell(ax, j, i, {'single': metrics})
            else:
                if genome_key in all_genome_reactions and reaction in all_genome_reactions[genome_key]:
                    rxn_data = all_genome_reactions[genome_key][reaction]

                    any_hit = any(m.get('bitscore', 0) > 0 for m in rxn_data['subunits'].values())

                    if any_hit:
                        draw_complex_cell(ax, j, i, rxn_data['subunits'])

    ax.set_xlim(-0.5, n_reactions - 0.5)
    ax.set_ylim(-0.5, n_organisms - 0.5)
    ax.set_xticks(range(n_reactions))
    ax.set_yticks(range(n_organisms))

    organism_labels = [org.replace('_', ' ') for org in compound_organisms]
    ax.set_yticklabels(organism_labels, fontsize=12, fontweight='bold')

    for i, organism in enumerate(compound_organisms):
        phylum = ORGANISM_PHYLUM.get(organism, 'Unknown')
        color = PHYLUM_COLORS.get(phylum, '#666666')
        ax.get_yticklabels()[i].set_color(color)

    ax.set_xticklabels(all_reactions, rotation=45, ha='right', fontsize=12, fontweight='bold')

    ax.grid(False)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.tick_params(length=0)

    ax.set_xlabel('REACTIONS', fontsize=12, fontweight='bold', color=DEEPSEA_PRIMARY['midnight'])
    ax.set_ylabel('ORGANISMS', fontsize=12, fontweight='bold', color=DEEPSEA_PRIMARY['midnight'])

    title_text = f'{compound} Metabolic Pathways'
    fig.suptitle(title_text, fontsize=16, fontweight='bold', color=DEEPSEA_PRIMARY['abyss'], y=0.98)

    pathway_legend_elements = []
    for pwy in pathways[:10]:
        if pwy in pathway_data:
            name = pathway_data[pwy]['name']
            color = pathway_color_map.get(pwy, '#666666')
            pathway_legend_elements.append(
                patches.Patch(facecolor=color, edgecolor='black', linewidth=0.5, label=name)
            )

    if compound in SPECIFIC_REACTIONS:
        pathway_legend_elements.append(
            patches.Patch(facecolor=pathway_color_map['SPECIFIC'], edgecolor='black',
                          linewidth=0.5, label=f'{", ".join(SPECIFIC_REACTIONS[compound])}')
        )

    if compound_transporter_ids:
        pathway_legend_elements.append(
            patches.Patch(facecolor=pathway_color_map['TRANSPORTER'], edgecolor='black',
                          linewidth=0.5, label=f'Transporters ({len(compound_transporter_ids)} types)')
        )

    pathway_legend = ax.legend(handles=pathway_legend_elements,
                               loc='upper left', bbox_to_anchor=(1.02, 1),
                               title='PATHWAYS', fontsize=12, ncol=2,
                               frameon=True, fancybox=False, shadow=True,
                               edgecolor=DEEPSEA_PRIMARY['ocean'],
                               facecolor=DEEPSEA_SECONDARY['shell'])
    pathway_legend.get_title().set_fontweight('bold')
    pathway_legend.get_title().set_fontsize(12)

    ax.add_artist(pathway_legend)
    color_legend_elements = [
        patches.Patch(facecolor=get_ocean_themed_color(50), edgecolor='black',
                      linewidth=0.5, label='0-100 (Weak)'),
        patches.Patch(facecolor=get_ocean_themed_color(150), edgecolor='black',
                      linewidth=0.5, label='100-200 (Moderate)'),
        patches.Patch(facecolor=get_ocean_themed_color(250), edgecolor='black',
                      linewidth=0.5, label='200+ (Strong)'),
    ]

    color_legend = ax.legend(handles=color_legend_elements,
                             loc='upper left', bbox_to_anchor=(1.02, 0.6),
                             title='BITSCORE STRENGTH', fontsize=12,
                             frameon=True, fancybox=False, shadow=True,
                             edgecolor=DEEPSEA_PRIMARY['ocean'],
                             facecolor=DEEPSEA_SECONDARY['shell'])
    color_legend.get_title().set_fontweight('bold')
    color_legend.get_title().set_fontsize(12)

    ax.add_artist(color_legend)
    complex_legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                   markersize=10, label='Single/Complete'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='none',
                   markeredgecolor='black', markersize=10, label='Split = Complex',
                   markeredgewidth=1.5),
        plt.Line2D([0], [0], marker='$!$', color='red', markersize=8,
                   linestyle='', label='Incomplete'),
        patches.Patch(facecolor=DEEPSEA_SECONDARY['fog'], edgecolor='black',
                      linewidth=0.5, label='Missing subunit'),
    ]

    complex_legend = ax.legend(handles=complex_legend_elements,
                               loc='upper left', bbox_to_anchor=(1.02, 0.35),
                               title='COMPLEX STATUS', fontsize=12,
                               frameon=True, fancybox=False, shadow=True,
                               edgecolor=DEEPSEA_PRIMARY['ocean'],
                               facecolor=DEEPSEA_SECONDARY['shell'])
    complex_legend.get_title().set_fontweight('bold')
    complex_legend.get_title().set_fontsize(12)

    ax.add_artist(complex_legend)
    phylum_legend_elements = [
        patches.Patch(facecolor=PHYLUM_COLORS['Pseudomonadota'], edgecolor='black',
                      linewidth=0.5, label='Pseudomonadota'),
        patches.Patch(facecolor=PHYLUM_COLORS['Acidobacteriota'], edgecolor='black',
                      linewidth=0.5, label='Acidobacteriota'),
        patches.Patch(facecolor=PHYLUM_COLORS['Bacteroidota'], edgecolor='black',
                      linewidth=0.5, label='Bacteroidota'),
        patches.Patch(facecolor=PHYLUM_COLORS['Actinomycetota'], edgecolor='black',
                      linewidth=0.5, label='Actinomycetota'),
    ]

    phylum_legend = ax.legend(handles=phylum_legend_elements,
                              loc='upper left', bbox_to_anchor=(1.02, 0.1),
                              title='PHYLUM', fontsize=12,
                              frameon=True, fancybox=False, shadow=True,
                              edgecolor=DEEPSEA_PRIMARY['ocean'],
                              facecolor=DEEPSEA_SECONDARY['shell'])
    phylum_legend.get_title().set_fontweight('bold')
    phylum_legend.get_title().set_fontsize(12)

    output_file = os.path.join(output_dir, f"pathway_{compound.lower()}_matrix_complex_viz.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.3)

    plt.show()
    plt.close()

    return True


# ==================== MAIN FUNCTION ====================

def main():
    """Main execution function"""

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load pathway metadata
    try:
        pathways_df = pd.read_csv(PATHWAYS_FILE, sep='\t', on_bad_lines='skip', engine='python')
        pathways_df['id'] = pathways_df['id'].str.strip('|')

        pathway_data = {}
        for _, row in pathways_df.iterrows():
            pwy_id = row['id']
            pwy_name = row['name']
            reactions = []
            if pd.notna(row['reaId']):
                reactions = [r.strip() for r in str(row['reaId']).split(',')]
            pathway_data[pwy_id] = {
                'name': pwy_name,
                'reactions': reactions
            }
    except Exception as e:
        pathway_data = {}

    # Load genome reactions
    all_genome_reactions = load_genome_reactions_with_complex_tracking(REACTIONS_DIR)

    if not all_genome_reactions:
        return

    # Load genome transporters
    all_genome_transporters = load_genome_transporters(TRANSPORTERS_DIR)

    if not all_genome_transporters:
        all_genome_transporters = {}

    # Process each compound
    total_plots_created = 0
    for compound, pathways in PATHWAY_GROUPS.items():
        try:
            result = create_tiled_heatmap_with_pathways(
                compound, pathways, all_genome_reactions,
                all_genome_transporters, pathway_data, OUTPUT_DIR
            )
            if result:
                total_plots_created += 1
        except Exception as e:
            pass

    return total_plots_created



if __name__ == "__main__":
    main()