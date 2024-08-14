PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 3

COMP       = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = FALSE

Bpack   := ./Make.package
Blocs   := .

include ../Make.mflo

ifeq ($(USE_SUNDIALS),TRUE)
# NOTE: SUNDIALS_ROOT must point to the directory where sundials is installed
# A good check is to see if $(SUNDIALS_ROOT)/lib has a bunch of libsundials_ files
# To run with sundials, enabled, please compile with USE_SUNDIALS = TRUE
SUNDIALS_ROOT ?= $(TOP)/../../../../../sundials/instdir
SUNDIALS_LIB_DIR ?= $(SUNDIALS_ROOT)/lib64

USE_CVODE_LIBS ?= TRUE
USE_ARKODE_LIBS ?= TRUE

DEFINES += -DAMREX_USE_SUNDIALS
INCLUDE_LOCATIONS += $(SUNDIALS_ROOT)/include
LIBRARY_LOCATIONS += $(SUNDIALS_LIB_DIR)

LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_cvode -lsundials_arkode
LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_core -lsundials_nvecserial
LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_nvecmanyvector
endif
