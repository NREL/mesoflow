# Mesoflow
## A mesoscale modeling tool for heterogenous gas-solid reacting flows

Mesoflow is a continuum scale simulation tool developed specifically for modeling transport and chemistry at the mesoscale. 
Our solver utilizes Cartesian block-structured adaptive mesh refinement to resolve complex surface morphologies 
(of catalysts/biomass particles among others) directly obtained from X-ray tomography data. 
An immersed boundary based formulation enables rapid representation of complex geometries prevalent in most mesoporous interfaces. 
The solver is developed on top of open-source performance portable library, AMReX, providing parallel execution capabilities on 
current and upcoming high-performance-computing (HPC) architectures. 

### Mesoscale simulations zeolite catalyst deactivation during catalytic upgrading of biomass pyrolysis vapors
see Sitaraman et al.,AIChE 2020 (https://www.nrel.gov/docs/fy21osti/78376.pdf) and
Ciesielski et al.,Energy & Fuels 35 (18), 14382-14400	6	2021 (https://pubs.acs.org/doi/full/10.1021/acs.energyfuels.1c02163)
![mesoflow_pic](https://user-images.githubusercontent.com/7399475/185490891-67a7d1c9-bfd4-4691-b735-2ff12ef89b12.png)

### Mesoscale simulation of flow through a porous biomass particle
see Crowley et al.,Frontiers in Energy Research 10 (https://www.frontiersin.org/articles/10.3389/fenrg.2022.850630/full)
![mesoscale_permeability](https://user-images.githubusercontent.com/7399475/185494788-904e5b21-62bc-4b65-af10-b9ff08440a43.png)




# Build instructions

* gcc and an MPI library (openMPI/MPICH) for CPU builds. cuda-11.0 is also required for GPU builds
* This solver depends on the AMReX library - clone from https://github.com/AMReX-Codes/amrex
* Set AMREX_HOME to the path of amrex (e.g. $export AMREX_HOME=/path/to/amrex
* go to any of the test cases in tests or model folder (e.g. cd models/blockCatalyst1d)
* build executable using the GNUMakefile - do $make for CPU build or do $make USE_CUDA=TRUE for GPU build

# Run instructions

* By default MPI is enabled in all builds, you can turn it off by doing $make USE_MPI=FALSE
* For parallel execution do $mpirun -n <procs> mesoflow3d.gnu.MPI.ex inputs
* For serial builds do $./mesoflow3d.gnu.ex inputs
* For GPU execution make sure the number of ranks match the number of GPUs on the machine. 
  For example, if you have 2 GPUs on a node, do $mpirun -n 2 mesoflow3d.gnu.MPI.CUDA.ex inputs

