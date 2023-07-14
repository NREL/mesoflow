#########################################################################
#  Write mesoflow simulation files from a pre-defined kinetic scheme
#  Requires two csv files (see example directory for expected format):
#   - `species.csv` containing species names and properties
#   - `reactions.csv` containing reaction definitions and stoichiometry
#  Writes the following mesoflow files:
#   - species.cpp & species.H containing species definitions and reactions
#   - transport.H containing transport properties
#   - thermo.H  containing thermodynamic properties
#   - userfuncs.cpp containing user-defined variables
#
#   Written by Meagan Crowley, NREL 2023
#########################################################################

import rxnscheme_writer
import os
import pandas as pd

# Read species.csv and reactions.csv into dataframes
project_dir='./gas_phase'

if not os.path.isdir(project_dir):
    os.makedirs(project_dir)

reactions_df = pd.read_csv(f'{project_dir}/reactions.csv')
species_df = pd.read_csv(f'{project_dir}/species.csv')
soldiff = '1.e-10' # solid diffusion coefficient, SI units
fludiff = '4.27e-5' # fluid diffusion coefficient, SI units

# Write species.cpp
rxnscheme_writer.write_speciescpp(species_df, project_dir)
rxnscheme_writer.write_rxn_arrays(reactions_df, project_dir)

# Write transport.H
rxnscheme_writer.write_transport(species_df, soldiff, fludiff, project_dir)

# write species.H
rxnscheme_writer.write_speciesH(species_df, reactions_df, project_dir)

# write thermo.H
rxnscheme_writer.write_thermo(species_df, project_dir)

# write userfuncs.cpp
rxnscheme_writer.write_userfuncs(species_df, project_dir)
