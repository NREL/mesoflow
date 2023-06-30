# rxnscheme_writer

This tool is used to convert tabulated reaction mechanisms to mesoflow simulation files.

# To run:
- Run with python3
- Required libraries: pandas
- Run write-mesoflow-files.py with python, altering the "project_dir" path and solid/gas diffusion coefficients to your use case.
- Two .csv files are required:
    - `reactions.csv` containing a matrix of reactions in your mechanism with stoichiometric coefficients. The convention assumes negative coefficients for reactants and positive for products
    - `species.csv` containing information about each species in the mechanism. All units are assumed SI by convention. For solid species, advect_flag should be set to 0 and gamma (heat capacity ratio) should be set to NaN. Gas-phase species should have advect_flag set to 1 and gamma can be user-specified.
  
    **Important requirements for species.csv:**
    - Gas species should be defined first
    - The first species is assumed to be the primary vapor
    - See example directory for the expected formatting and sample simulation files


# Output
write-mesoflow-files.py will write the following files with the information from your supplied csv files:
  - species.cpp, species.H, transport.H, thermo.H, and userfuncs.cpp
 These can be copied into an existing mesoflow simulation directory. Compile a new executable with these new files and run your simulation.

# Disclaimer
These have been tested using the simulation files supplied in the example directory, but may not be compatible with other use cases. Please submit an issue on GitHub or email meagan.crowley@nrel.gov if you run into issues.
