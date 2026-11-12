#===============================================================================
# Makefile for TDGP3D -- the 3D trans-dimensional Gaussian-process (TDGP) Bayesian
# MT inversion (tdgp3d_*.f90), with the V6 3D MT forward code as the forward
# operator.
#
# The V6 forward sources are compiled from MT3D_DIR (default ../CODES_MPIV6, i.e.
# this folder next to CODES_MPIV6) exactly as they are; object and .mod files go
# in this folder and nothing is written to MT3D_DIR.
#
# Usage:  make CLUSTER=olivia <target>          (CLUSTER=betzy: the same settings)
#         make CLUSTER=olivia MT3D_DIR=/path/to/CODES_MPIV6 <target>
#
#   tdgp3d_test_gp   standalone 3D GP test, no MPI, no forward code:  ./tdgp3d_test_gp
#   tdgp3d_single    one chain, no parallel tempering:  mpirun -np <nproc> ./tdgp3d_single
#   tdgp3d_mpi       parallel-tempered sampler:  mpirun -np <1+N_chains*nproc_per_team> ./tdgp3d_mpi
#   clean            removes this folder's objects, .mod files and the executables
#
# Use plain make (no -j): the object lists are in module-dependency order.
#
# The run folder needs tdgp.settings (see tdgp.settings.example), MT3D_V6.dat (the
# forward mesh, V6 format) and the ModEM data file named by dataFile.
#===============================================================================

MT3D_DIR ?= ../CODES_MPIV6

CLUSTER_LC := $(shell echo $(CLUSTER) | tr A-Z a-z)

# Olivia and Betzy (Sigma2), Intel toolchain (module load intel/2023b), as the V6
# Makefile. MKL provides LAPACK (for the GP) and the FFTW3 Fortran interface that
# V6's FFT_2D calls.
ifneq (,$(filter $(CLUSTER_LC),olivia betzy))
   FC     = mpiifort
   FFLAGS = -O2 -fpp -qopenmp
   LIBS   = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
else ifneq (clean,$(MAKECMDGOALS))
   $(error CLUSTER not given or not known. Usage: make CLUSTER=olivia <target> (or CLUSTER=betzy). Targets: tdgp3d_test_gp tdgp3d_single tdgp3d_mpi clean)
endif

default: tdgp3d_mpi

#-------------------------------------------------------------------------------
# V6 forward objects (from $(MT3D_DIR)): the modules and the per-frequency solve
# only -- not the V6 driver, mesh_master or the output routines.
#-------------------------------------------------------------------------------
V6_OBJS = MODULE.o MODULE_MESHMASTER.o interp_util.o MODEL.o READ_PARA.o SZ_ALLOCATE.o \
          WAVE_NUMBER.o GAUSSA.o WAVE_MIGRATION.o GS_FFT_WAVE_NUMBER.o GMATR1.o csem_x_1d.o \
          GMATR2.o GaussPoint.o WaveNumberSequence.o FFT_2D.o GaussFFT_2D.o \
          fxfyfz0_x.o MT3DA.o fxfyfz_x.o fxfyfz0_y.o fxfyfz_y.o

#-------------------------------------------------------------------------------
# TDGP3D objects. tdgp_data_type.f90 and tdgp_mpi_topology.f90 are copies of the 2D
# code's files (mare2dem_source_Bayesian).
#-------------------------------------------------------------------------------
TDGP3D_OBJS = tdgp_data_type.o tdgp3d_gp.o tdgp3d_modem.o tdgp3d_com.o tdgp_mpi_topology.o \
              tdgp3d_lib_interface.o tdgp3d_rjmcmc.o tdgp3d_pt.o

tdgp3d_test_gp: tdgp_data_type.o tdgp3d_gp.o
	$(FC) $(FFLAGS) $^ $(LIBS) -o $@

tdgp3d_single: $(V6_OBJS) $(TDGP3D_OBJS)
	$(FC) $(FFLAGS) $^ $(LIBS) -o $@

tdgp3d_mpi: $(V6_OBJS) $(TDGP3D_OBJS) tdgp3d_main.o
	$(FC) $(FFLAGS) $^ $(LIBS) -o $@

clean:
	rm -f *.o *.mod *__genmod.f90 tdgp3d_mpi

#-------------------------------------------------------------------------------
# Compiling: this folder's sources first, then the V6 sources in $(MT3D_DIR).
#-------------------------------------------------------------------------------
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: $(MT3D_DIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

.PHONY: default clean
