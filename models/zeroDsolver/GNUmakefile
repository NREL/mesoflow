AMREX_HOME ?= ../../../../..

EBASE = zerod

PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE

DIM       = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE 
USE_CUDA   = FALSE

USE_SEC_ORDER_FLUX=TRUE
USE_IMM_BOUNDARY=TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Base/Make.package
include ./Make.package 

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
