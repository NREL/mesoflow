# Zero-D Solver


## How to run:
1. Copy species.cpp and species.H files from an existing mechanism into this directory. 
2. Change H2_ID and S1_ID in driver.cpp to match the primary vapor and active site species in your mechanism. Specify initial concentrations in mol/m^3
3. Change finaltime in driver.cpp if necessary, will run to 1.0 s by default
4. Compile & run executable
5. Plot concentrations when run has finished using plot-quants.py, changing the species names to match your mechanism. The order of columns is the order of your species as defined in species.H.

## Troubleshooting
- *Important*: Delete or rename quants.dat.0 in between runs or else it will have the data from the next run appended to it.
- If you see nans in quants.dat.0, reduce the timestep (dt in driver.cpp) or check your mechanism for errors

