AMREX_HOME ?= ../../../../..
MFLO_HOME = /scratch/mcrowley2/mesoflow

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM       = 3

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE 
USE_CUDA   = FALSE

USE_SEC_ORDER_FLUX=TRUE
USE_IMM_BOUNDARY=TRUE

Bpack   := ./Make.package 
Blocs   := . 

include ${MFLO_HOME}/tests/Make.mflo
