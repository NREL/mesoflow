# rxnscheme_writer

This tool is used to convert tabulated reaction mechanisms to mesoflow simulation files.

# To run:
- Run with python3
- Required libraries: pandas
- Two .csv files are required:
    - `reactions.csv` containing a matrix of reactions in your mechanism with stoichiometric coefficients. The convention assumes negative coefficients for reactants and positive for products
    - `species.csv` containing information about each species in the mechanism. All units are assumed SI by convention. For solid species, advect_flag should be set to 0 and gamma (heat capacity ratio) should be set to NaN. Gas-phase species should have advect_flag set to 1 and gamma can be user-specified.
    - See examples directory for the expected formatting
- Run the write-mesoflow-files.py with python, altering the "project_dir" path, and solid/gas diffusion coefficients to your use case.

# Output
write-mesoflow-files.py will write species.cpp, species.H, and transport.H with the information from your csv files. These can be copied into an existing mesoflow simulation directory.
