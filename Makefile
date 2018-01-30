SUBSYSNAME := phasta
MODULENAME := SimpleAdapt
ISEXECUTABLE := 1
OUTNAME := SimpleAdapt
DEVROOT := /projects/tools/Simmetrix.develop

#Specify compilers 
CC := mpicc
CXX := mpicxx
FC := gfortran

SIM := 1
DEFS := $(DEFS) -DSIM

MODELER := parasolid
NODEP := 1
NOSHARED := 1
MESHSIM := /projects/tools/SimmetrixTest/12.0-171109dev
PARALLEL := openmpi

ARCHOS := x86_64_linux
ifeq ($(ARCHOS),x86_64_linux)
   SIMARCH := x64_rhel7_gcc48
   ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread
   DEFS := $(DEFS) -DB64
endif
dirs := .

ifeq ($(SIM),1)
 ifeq ($(PARALLEL),mpich1)
    SIMPAR := -mpich
 endif
 ifeq ($(PARALLEL),mvapich-1.2-gcc41)
    SIMPAR := -mpich
 endif
 ifeq ($(PARALLEL),openmpi)
   SIMPAR := -openmpi
 endif
 ifeq ($(PARALLEL),sgimpi)
  SIMPAR := -mpisgi
 endif
endif

DEPS := $(DEPS) \
        phasta/phShape \
        phasta/phUtil/LU \
        phasta/phastaIO
LIBS := $(LIBS) phShape LU phastaIO

ifeq ($(SIM),1)
  ifeq ($(MODELER),parasolid)
    PARASOLID = $(MESHSIM)/lib/$(SIMARCH)/psKrnl
    PSKRNL = -lpskernel
    MESHSIMMODELVERS = -lSimParasolid280
  else
    MODELLIB = -lSim$(MODELER)
  endif

  INCLUDES := $(INCLUDES) -I$(MESHSIM)/include
  DEFS := $(DEFS) -DSIM -DSIM_$(shell echo $(MODELER)| tr '[a-z]' '[A-Z]')

  MESHSIMFLAGS := -L$(MESHSIM)/lib/$(SIMARCH)
ifeq ($(MODELER),parasolid)
  MESHSIMFLAGS := $(MESHSIMFLAGS)  -L$(PARASOLID)
endif
  MESHSIMFLAGS := $(MESHSIMFLAGS) \
                -lSimPartitionedMesh-mpi -lSimAdvMeshing -lSimMeshing \
                -lSimPartitionWrapper-$(PARALLEL) -lSimPartitionedMesh-mpi $(MESHSIMMODELVERS) -lSimMeshTools -lSimModel\
                $(PSKRNL)
#               -lSimPartitionWrapper-$(PARALLEL)14 -lSimPartitionedMesh-mpi $(MESHSIMMODELVERS) -lSimMeshTools -lSimModel\
                $(PSKRNL) 
endif

ifeq ($(SIM),1)
LDFLAGS := $(LDFLAGS) $(MESHSIMFLAGS) $(ETCFLAGS) -lnsl -lm -lpthread
endif

# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common
