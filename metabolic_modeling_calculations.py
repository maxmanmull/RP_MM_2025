#!/usr/bin/env python3
"""
Taurine Metabolic Modeling/Data Extraction and Analysis
"""
import os
import glob
import pandas as pd
import numpy as np
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
from collections import defaultdict
import warnings

warnings.filterwarnings('ignore')

# ==================== FILE PATHS ====================
MODELS_DIR = r".\gapseq_models"
DRAFT_DIR = r".\11_gapseq_draft"
SEED_PWY_FILE = r".\seed_pwy.tbl"  #get from gapseq git repo

# ==================== OUTPUT FILES ====================
OUTPUT_GROWTH_RATES = 'growth_rates.csv'
OUTPUT_EXCHANGE_FLUXES = 'exchange_fluxes.csv'
OUTPUT_COMPREHENSIVE = 'reaction_analysis.csv'

# ==================== CONSTANTS ====================
FLUX_THRESHOLD = 1e-6
FVA_FRACTION = 1.0
KO_SENSITIVITY_THRESHOLD = 0.05
CONSTRAINEDNESS_THRESHOLDS = {
    'highly_constrained': 0.1,  # Less than 10% variability
    'moderately_constrained': 0.5  # Less than 50% variability
}

# Compound name mappings for standardization
COMPOUND_NAME_MAPPINGS = {
    'taurine': 'Taurine', 'taur': 'Taurine', 'tau': 'Taurine',
}

# Patterns to remove from compound names
PATTERNS_TO_REMOVE = [
    '-e0', '_e0', '[e0]', '[e]', '_e', '-e',
    'EX_', 'M_', 'EX ', '__'
]


# ==================== ANALYZER CLASS ====================

class TaurineModelAnalyzer:
    def __init__(self, models_dir, draft_dir, seed_pwy_file=None):
        self.models_dir = models_dir
        self.draft_dir = draft_dir
        self.models = {}
        self.draft_models = {}
        self.results = {}
        self.exchange_fluxes = {}
        self.gapfill_stats = {}
        self.gapfill_reactions = {}
        self.reaction_info = {}
        self.all_fluxes = {}

        if seed_pwy_file:
            self.load_seed_pathways(seed_pwy_file)

    def load_models(self):
        """Load all Taurine SBML models"""
        taurine_files = glob.glob(os.path.join(self.models_dir, "*Taurine.xml"))

        for filepath in sorted(taurine_files):
            filename = os.path.basename(filepath)
            org_id = filename.replace('_Taurine.xml', '').replace('Taurine.xml', '')
            org_id = org_id.replace('_genome', '').replace('genome', '')

            try:
                model = read_sbml_model(filepath)
                model.solver = 'glpk'

                biomass_rxn = None
                for rxn in model.reactions:
                    if 'biomass' in rxn.id.lower() or rxn.objective_coefficient != 0:
                        biomass_rxn = rxn.id
                        break

                self.models[org_id] = {
                    'model': model,
                    'biomass_rxn': biomass_rxn,
                    'n_reactions': len(model.reactions),
                    'n_metabolites': len(model.metabolites)
                }

            except Exception as e:
                pass

    def load_draft_models(self):
        """Load draft models for gap-filling analysis"""
        for org_id in self.models.keys():
            patterns = [
                f"{org_id}-draft.xml",
                f"{org_id}_genome-draft.xml",
                f"{org_id}genome-draft.xml"
            ]

            for pattern in patterns:
                draft_path = os.path.join(self.draft_dir, pattern)
                if os.path.exists(draft_path):
                    try:
                        draft = read_sbml_model(draft_path)
                        draft.solver = 'glpk'
                        self.draft_models[org_id] = draft
                        break
                    except:
                        pass

    def run_fba_analysis(self):
        """Run FBA/pFBA analysis"""
        for org_id, model_data in self.models.items():
            model = model_data['model']
            biomass_rxn = model_data['biomass_rxn']

            try:
                fba_solution = model.optimize()

                if fba_solution.status == 'optimal':
                    growth_rate = fba_solution.objective_value
                    solution_to_use = fba_solution
                    method = "FBA"

                    try:
                        pfba_solution = pfba(model)
                        if pfba_solution.status == 'optimal':
                            if biomass_rxn:
                                pfba_growth = pfba_solution.fluxes[biomass_rxn]
                                solution_to_use = pfba_solution
                                method = "pFBA"
                                growth_rate = pfba_growth
                    except Exception as e:
                        pass

                    self.results[org_id] = {
                        'growth_rate': growth_rate,
                        'method': method,
                        'status': solution_to_use.status,
                        'fluxes': solution_to_use.fluxes.to_dict(),
                        'biomass_rxn': biomass_rxn
                    }

                    self.all_fluxes[org_id] = solution_to_use.fluxes.to_dict()
                    self.extract_exchange_fluxes(org_id, model, solution_to_use.fluxes)

                else:
                    self.results[org_id] = {
                        'growth_rate': 0,
                        'method': 'Failed',
                        'status': fba_solution.status,
                        'fluxes': {},
                        'biomass_rxn': biomass_rxn
                    }

            except Exception as e:
                self.results[org_id] = {
                    'growth_rate': 0,
                    'method': 'Error',
                    'status': 'error',
                    'fluxes': {},
                    'biomass_rxn': biomass_rxn
                }

    def extract_exchange_fluxes(self, org_id, model, fluxes):
        """Extract and clean exchange reaction fluxes"""
        exchanges = {}

        for rxn in model.reactions:
            if rxn.boundary:
                flux = fluxes[rxn.id]

                if abs(flux) > FLUX_THRESHOLD:
                    if rxn.metabolites:
                        met = list(rxn.metabolites.keys())[0]
                        compound = self.clean_compound_name(met.name if met.name else met.id)
                    else:
                        compound = self.clean_compound_name(rxn.name if rxn.name else rxn.id)

                    if compound in exchanges:
                        exchanges[compound] += flux
                    else:
                        exchanges[compound] = flux

        self.exchange_fluxes[org_id] = exchanges

    def clean_compound_name(self, name):
        """Clean and standardize compound names"""
        if not name:
            return "Unknown"

        name = str(name)

        for pattern in PATTERNS_TO_REMOVE:
            name = name.replace(pattern, ' ')

        name = ' '.join(name.replace('_', ' ').split())

        name_lower = name.lower().strip()
        for old, new in COMPOUND_NAME_MAPPINGS.items():
            if name_lower == old or name_lower == f"{old} e0":
                return new

        if name and name[0].islower():
            name = name.capitalize()

        return name

    def analyze_gapfilling_with_fva(self):
        """Analyze gap-filled reactions using FVA with detailed constrainedness metrics"""
        from cobra.flux_analysis import flux_variability_analysis

        for org_id in self.models.keys():
            if org_id not in self.draft_models or org_id not in self.results:
                continue

            model = self.models[org_id]['model']
            draft = self.draft_models[org_id]
            fluxes = self.results[org_id]['fluxes']
            baseline_growth = self.results[org_id]['growth_rate']

            if baseline_growth <= 0:
                continue

            model_rxns = set(r.id for r in model.reactions)
            draft_rxns = set(r.id for r in draft.reactions)
            gapfilled = model_rxns - draft_rxns

            active_gapfilled = set()
            inactive_gapfilled = set()

            for rxn_id in gapfilled:
                if rxn_id in fluxes and abs(fluxes[rxn_id]) > FLUX_THRESHOLD:
                    active_gapfilled.add(rxn_id)
                else:
                    inactive_gapfilled.add(rxn_id)

            if not active_gapfilled:
                continue

            fva_100 = flux_variability_analysis(
                model,
                reaction_list=list(active_gapfilled),
                fraction_of_optimum=FVA_FRACTION,
                processes=1
            )

            highly_constrained = set()
            moderately_constrained = set()
            flexible = set()
            fva_details = {}

            for rxn_id in active_gapfilled:
                observed_flux = fluxes.get(rxn_id, 0)

                min_100 = fva_100.loc[rxn_id, 'minimum']
                max_100 = fva_100.loc[rxn_id, 'maximum']
                range_100 = max_100 - min_100
                relative_var = range_100 / (abs(observed_flux) + 1e-9)

                if relative_var < CONSTRAINEDNESS_THRESHOLDS['highly_constrained']:
                    highly_constrained.add(rxn_id)
                elif relative_var < CONSTRAINEDNESS_THRESHOLDS['moderately_constrained']:
                    moderately_constrained.add(rxn_id)
                else:
                    flexible.add(rxn_id)

                fva_details[rxn_id] = {
                    'observed_flux': observed_flux,
                    'fva_100_min': min_100,
                    'fva_100_max': max_100,
                    'fva_100_range': range_100,
                    'relative_variability': relative_var,
                    'highly_constrained': rxn_id in highly_constrained,
                    'moderately_constrained': rxn_id in moderately_constrained,
                    'flexible': rxn_id in flexible,
                    'is_gapfilled': True,
                    'is_active': True
                }

                if rxn_id not in self.reaction_info:
                    try:
                        rxn = model.reactions.get_by_id(rxn_id)
                        self.reaction_info[rxn_id] = {
                            'name': rxn.name if rxn.name else rxn_id,
                            'subsystem': rxn.subsystem if rxn.subsystem else "Unknown",
                            'reaction': rxn.reaction if hasattr(rxn, 'reaction') else "",
                            'gene_rule': rxn.gene_reaction_rule if hasattr(rxn, 'gene_reaction_rule') else "",
                            'lower_bound': rxn.lower_bound,
                            'upper_bound': rxn.upper_bound
                        }
                    except:
                        self.reaction_info[rxn_id] = {
                            'name': rxn_id,
                            'subsystem': "Unknown",
                            'reaction': "",
                            'gene_rule': "",
                            'lower_bound': 0,
                            'upper_bound': 0
                        }

            ko_sensitive = set()
            ko_details = {}

            for rxn_id in highly_constrained:
                try:
                    with model:
                        model.reactions.get_by_id(rxn_id).knock_out()
                        ko_solution = model.optimize()

                        if ko_solution.status == 'optimal':
                            ko_growth = ko_solution.objective_value
                            reduction = (baseline_growth - ko_growth) / baseline_growth if baseline_growth > 0 else 1

                            ko_details[rxn_id] = {
                                'ko_growth': ko_growth,
                                'growth_reduction': reduction,
                                'is_ko_sensitive': reduction > KO_SENSITIVITY_THRESHOLD
                            }

                            if reduction > KO_SENSITIVITY_THRESHOLD:
                                ko_sensitive.add(rxn_id)
                        else:
                            ko_details[rxn_id] = {
                                'ko_growth': 0,
                                'growth_reduction': 1.0,
                                'is_ko_sensitive': True
                            }
                            ko_sensitive.add(rxn_id)
                except:
                    ko_details[rxn_id] = {
                        'ko_growth': None,
                        'growth_reduction': None,
                        'is_ko_sensitive': None
                    }

            for rxn_id, ko_info in ko_details.items():
                if rxn_id in fva_details:
                    fva_details[rxn_id].update({
                        'ko_tested': True,
                        'ko_growth': ko_info['ko_growth'],
                        'growth_reduction': ko_info['growth_reduction'],
                        'ko_sensitive': ko_info['is_ko_sensitive']
                    })

            for rxn_id in fva_details:
                if rxn_id not in ko_details:
                    fva_details[rxn_id].update({
                        'ko_tested': False,
                        'ko_growth': None,
                        'growth_reduction': None,
                        'ko_sensitive': False
                    })

            self.gapfill_reactions[org_id] = {
                'all': gapfilled,
                'active': active_gapfilled,
                'inactive': inactive_gapfilled,
                'highly_constrained': highly_constrained,
                'moderately_constrained': moderately_constrained,
                'flexible': flexible,
                'ko_sensitive': ko_sensitive,
                'fva_details': fva_details
            }

    def create_comprehensive_reaction_analysis(self):
        """Create a reaction analysis showing all unique reactions and their classification"""
        all_reactions = set()

        for org_id, model_data in self.models.items():
            model = model_data['model']
            for rxn in model.reactions:
                all_reactions.add(rxn.id)

                if rxn.id not in self.reaction_info:
                    self.reaction_info[rxn.id] = {
                        'name': rxn.name if rxn.name else rxn.id,
                        'subsystem': rxn.subsystem if rxn.subsystem else "Unknown",
                        'reaction': rxn.reaction if hasattr(rxn, 'reaction') else "",
                        'gene_rule': rxn.gene_reaction_rule if hasattr(rxn, 'gene_reaction_rule') else "",
                        'lower_bound': rxn.lower_bound,
                        'upper_bound': rxn.upper_bound
                    }

        reaction_analysis = []

        for rxn_id in sorted(all_reactions):
            rxn_info = self.reaction_info.get(rxn_id, {
                'name': rxn_id,
                'subsystem': 'Unknown',
                'reaction': '',
                'gene_rule': ''
            })

            for org_id in sorted(self.models.keys()):
                model = self.models[org_id]['model']
                reaction_exists = any(r.id == rxn_id for r in model.reactions)

                if not reaction_exists:
                    continue

                is_gapfilled = False
                is_gapfilled_active = False
                is_gapfilled_active_highly_constrained = False
                is_gapfilled_active_highly_constrained_ko_sensitive = False
                observed_flux = 0
                growth_reduction = None

                if org_id in self.gapfill_reactions:
                    gf_data = self.gapfill_reactions[org_id]

                    if rxn_id in gf_data.get('all', set()):
                        is_gapfilled = True

                        if rxn_id in gf_data.get('active', set()):
                            is_gapfilled_active = True

                            if org_id in self.all_fluxes:
                                observed_flux = self.all_fluxes[org_id].get(rxn_id, 0)

                            if rxn_id in gf_data.get('highly_constrained', set()):
                                is_gapfilled_active_highly_constrained = True

                                if rxn_id in gf_data.get('ko_sensitive', set()):
                                    is_gapfilled_active_highly_constrained_ko_sensitive = True

                                    if 'fva_details' in gf_data and rxn_id in gf_data['fva_details']:
                                        growth_reduction = gf_data['fva_details'][rxn_id].get('growth_reduction')

                reaction_analysis.append({
                    'Reaction_ID': rxn_id,
                    'Organism': org_id,
                    'Reaction_Name': rxn_info['name'],
                    'Subsystem': rxn_info['subsystem'],
                    'Gene_Rule': rxn_info.get('gene_rule', ''),
                    'Is_Gapfilled': is_gapfilled,
                    'Is_Gapfilled_Active': is_gapfilled_active,
                    'Is_Gapfilled_Active_Highly_Constrained': is_gapfilled_active_highly_constrained,
                    'Is_Gapfilled_Active_Highly_Constrained_KO_Sensitive': is_gapfilled_active_highly_constrained_ko_sensitive,
                    'Observed_Flux': observed_flux,
                    'Growth_Reduction_on_KO': growth_reduction
                })

        df_comprehensive = pd.DataFrame(reaction_analysis)
        df_comprehensive = df_comprehensive.sort_values(['Organism', 'Reaction_ID'])

        return df_comprehensive

    def load_seed_pathways(self, seed_pwy_file):
        """Load SEED pathway mappings from seed_pwy.tbl file"""
        self.seed_pathways = {}
        self.reaction_to_pathways = defaultdict(list)
        self.pathway_info = {}

        try:
            seed_df = pd.read_csv(seed_pwy_file, sep='\t', dtype=str)
            return True
        except FileNotFoundError:
            return False
        except Exception as e:
            return False

    def save_streamlined_results(self):
        """Save results """

        # 1. Growth rates
        growth_data = []
        for org, result in self.results.items():
            model_data = self.models.get(org, {})
            growth_data.append({
                'Organism': org,
                'Growth_Rate_hr-1': result['growth_rate'],
                'Method': result['method'],
                'Status': result['status'],
                'Biomass_Reaction': result['biomass_rxn'],
                'Total_Reactions': model_data.get('n_reactions', 0),
                'Total_Metabolites': model_data.get('n_metabolites', 0)
            })

        if growth_data:
            df_growth = pd.DataFrame(growth_data)
            df_growth = df_growth.sort_values('Growth_Rate_hr-1', ascending=False)
            df_growth.to_csv(OUTPUT_GROWTH_RATES, index=False)

        # 2. Exchange fluxes
        flux_data = []
        for org, exchanges in self.exchange_fluxes.items():
            for compound, flux in exchanges.items():
                flux_data.append({
                    'Organism': org,
                    'Compound': compound,
                    'Flux_mmol/gDW/hr': flux,
                    'Direction': 'Uptake' if flux < 0 else 'Secretion',
                    'Abs_Flux': abs(flux)
                })

        if flux_data:
            df_flux = pd.DataFrame(flux_data)
            df_flux = df_flux.sort_values(['Organism', 'Abs_Flux'], ascending=[True, False])
            df_flux.to_csv(OUTPUT_EXCHANGE_FLUXES, index=False)

        # 3. reaction analysis
        df_comprehensive = self.create_comprehensive_reaction_analysis()
        if df_comprehensive is not None:
            df_comprehensive.to_csv(OUTPUT_COMPREHENSIVE, index=False)


# ==================== MAIN FUNCTION ====================

def main():
    """Main analysis pipeline"""

    analyzer = TaurineModelAnalyzer(MODELS_DIR, DRAFT_DIR)

    analyzer.load_seed_pathways(SEED_PWY_FILE)
    analyzer.load_models()

    if analyzer.models:
        analyzer.load_draft_models()
        analyzer.run_fba_analysis()

        if analyzer.exchange_fluxes:
            analyzer.analyze_gapfilling_with_fva()
            analyzer.save_streamlined_results()
            return True
        else:
            return False
    else:
        return False


# ==================== EXECUTION ====================

if __name__ == "__main__":
    main()