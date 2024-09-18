#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                                                       :::
#                          LTRANS Makefile                              :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#==========================================================================
#                         USER-DEFINED OPTIONS                            =
#==========================================================================
OMP :=     # turn on to enable openmp (turn it off to have a sequential executable)
Debug :=    # turn on to activate debug flags
#
#    Turn one of the following on to link the libraries and flags specific to a given cluster:
IFORT_UMCES :=     
GFORTRAN_UMCES :=    
GFORTRAN_USGS :=   
PGI_USGS   :=
Fluxus := 
INTEL_g100 := on      
#==========================================================================
#================   setting some default options, the =====================
#================   user shouldnt need to change them =====================
#==========================================================================
FOPENMP :=   #  
HDF5_LIBDIR= # by default, not needed
EXTRA_LIBS_AND_FLAGS= # by default, not needed

ifdef OMP
  FOPENMP := on 
endif
HIDE_UNUSED_VARIABLES := on # set-up only for fluxus 

#==========================================================================
# =======   from now on, modifications should be done only      ===========
# =======   to insert the configuration for a new cluster or    ===========
# =======   to to change the libraries linked by a given cluser ===========
#==========================================================================
ifdef IFORT_UMCES 
  FC = ifort
  Intel := on
  Gnu := 
  NETCDF_INCDIR = /usr/local/include
  NETCDF_LIBDIR = /usr/local/lib
  FFLAGS = -vec-report0 -fp-model precise -mcmodel=medium -I$(NETCDF_INCDIR)
  Debug=
endif

ifdef GFORTRAN_UMCES
  Intel := 
  Gnu := on
  FC = gfortran
  NETCDF_INCDIR = /usr/local/include
  NETCDF_LIBDIR = /usr/local/lib
  FFLAGS = -march=k8 -ffast-math -fno-cx-limited-range -O3 -funroll-loops --param max-unroll-times=4 -ffree-line-length-none -I$(NETCDF_INCDIR)
  Debug=
endif

ifdef GFORTRAN_USGS
  FC = gfortran
  Intel := 
  Gnu := on
  NETCDF_INCDIR = /usr/local/include
  NETCDF_LIBDIR = /usr/local/lib
  FFLAGS = -Wall  -Warray-bounds -fbacktrace -ffree-line-length-none  -I$(NETCDF_INCDIR)
  Debug=
endif

ifdef PGI_USGS
  FC = pgf90
  Intel := 
  Gnu := on
  NETCDF_INCDIR = /share/apps/netcdf/include
  NETCDF_LIBDIR = /share/apps/netcdf/lib
  HDF5_LIBDIR=/share/apps/hdf5/lib
  EXTRA_LIBS_AND_FLAGS= -lz -lm -L/share/apps/szip/lib -lsz
  FFLAGS := -g -I$(NETCDF_INCDIR)
  Debug=
endif

ifdef Fluxus
  Intel := 
  Gnu := on
  ifdef HIDE_UNUSED_VARIABLES
     F_HIDEUNUSED = -Wno-unused-parameter  -Wno-unused-dummy-argument  -Wno-unused
  else
     F_HIDEUNUSED = 
  endif
  FC = mpif90
  F90C = mpif90
  CC = mpicc
  NETCDF_INCDIR =  /usr/include -I/usr/lib64/gfortran/modules/ 
  NETCDF_LIBDIR = /usr/lib -L/usr/lib64  -lnetcdff -lnetcdf
  ifdef Debug 
    FFLAGS =  -O0 -pg -g  -fbounds-check -fimplicit-none -ffpe-trap=invalid,overflow -pedantic -ffixed-line-length-132 -cpp  $(F_HIDEUNUSED)   -I$(NETCDF_INCDIR)  
  else
    FFLAGS =   -O3  -Wall  -Warray-bounds -fbacktrace -ffree-line-length-none -cpp $(F_HIDEUNUSED) -I$(NETCDF_INCDIR) 
  endif
endif

ifdef INTEL_g100
  Intel := on
  Gnu := 
  DEBUGFLAGS= -O3 -xcore-avx2
  ifdef Intel
   FC = ifort
   F90C = ifort
   CC = icc
   ifdef Debug
     DEBUGFLAGS= -g -check all -traceback -debug all -check uninit -check bounds -check pointer -trace -qopt-report=5 -qopt-report-phase=all -fp-stack-check -fpe0 -ftrapuv -CB -assume byterecl
   endif
   ifdef DebugOptimized
     DEBUGFLAGS= -g -traceback -trace -fpe0 -ftrapuv   -ggdb -g3  -check all -debug all -assume byterecl
   endif
   ifdef DebugO0
     DEBUGFLAGS= -g -trace -traceback -fpe0 -fp-stack-check -ftrapuv -CB -qopt-report=5 -qopt-report-phase=all  -assume byterecl
   endif
  else
   FC = gfortran
   F90C = gfortran
   CC = gcc
   ifdef Debug
     DEBUGFLAGS= -g -fsanitize=address -fbacktrace 
   endif
  endif
    NETCDF_INCDIR = ${NETCDF_FORTRAN_HOME}/include -I ${NETCDF_C_HOME}/include
    NETCDF_LIBDIR = ${NETCDF_FORTRAN_HOME}/lib     -L ${NETCDF_C_HOME}/lib
    HDF5_LIBDIR = ${HDF5_HOME}/lib  -L ${CURL_HOME}/lib 
  #FFLAGS = -O3 -fno-math-errno -unroll=3 -cpp  -extend-source -assume byterecl -I$(NETCDF_INCDIR)
  FFLAGS =  -D_FILE_OFFSET_BITS=64 -O3 -funroll-loops -xcore-avx2 -axCORE-AVX2,AVX -fpp -assume byterecl -I$(NETCDF_INCDIR)
  ifdef Debug 
    FFLAGS = $(DEBUGFLAGS) -O0 -fno-omit-frame-pointer -cpp -ggdb -g3 -I$(NETCDF_INCDIR) 
  endif
endif



#------------------------------------------------
#    Set NetCDF Library Locations.
#    If NetCDF was compiled with HDF5, set:
#        HDF5 := on
#    Otherwise, leave blank:
#        HDF5 :=
#------------------------------------------------

HDF5 := on 

#==========================================================================
# End of user-defined options. Nothing should be changed below this point =
#==========================================================================

ifdef FOPENMP
  ifdef Intel
    FOPENMP_FLAG = -qopenmp
  else
    FOPENMP_FLAG = -fopenmp
  endif
  LTRANS_NAME=LTRANS_Zlev_OMP
else
  FOPENMP_FLAG =
  LTRANS_NAME=LTRANS_Zlev_SEQ
endif
LTRANSEXE=${LTRANS_NAME}.exe
ifdef Debug
  LTRANSEXE=${LTRANS_NAME}_dbg.exe
endif

COMPILINGMESSAGE1= "Executable is ../$(LTRANSEXE)"

ifdef FOPENMP
  COMPILINGMESSAGE2= "- Compilation included OPENMP flag $(FOPENMP_FLAG). To compile instead in sequential deactivate the OMP flag in the makefile"
else
  COMPILINGMESSAGE2= "- The executable is for sequential use only. To compile an OpenMP version activate the OMP flag in the makefile"
endif

ifdef Debug
  COMPILINGMESSAGE3= "- WARNING !!! NON-OPTIMIZED EXECUTABLE compiled with Debug flags!!! "
else
  COMPILINGMESSAGE3= "- No debugging flags used for the compilation, the executable is optimised. To create instead a debugging executable please activate the Debug flag in the makefile"
endif


OBJS          = gridcell_module.o interpolation_module.o parameter_module.o \
				point_in_polygon_module.o random_module.o tension_module.o  \
				conversion_module.o oil_module_Zlev.o stokes_drift_module.o  \
				hydrodynamic_module.o norm_module.o                         \
				boundary_module.o hor_turb_module.o settlement_module.o     \
				ver_turb_module.o stranding_module.o behavior_module.o                 

ifdef HDF5
	LIBS      = -L$(NETCDF_LIBDIR) -lnetcdf -lnetcdff -L$(HDF5_LIBDIR) -lhdf5_hl -lhdf5 $(EXTRA_LIBS_AND_FLAGS) -lcurl
else
	LIBS      = -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf
endif

LTRANS : $(OBJS)
	@echo "  Compiling LTRANS_vZ.f90"
	@echo "$(F90C) $(FOPENMP_FLAG)  $(FFLAGS) -o ../$(LTRANSEXE) LTRANS.f90 $(OBJS) $(LIBS)"
	@$(F90C) $(FOPENMP_FLAG)  $(FFLAGS) -o ../$(LTRANSEXE) LTRANS.f90 $(OBJS) $(LIBS)
	#@echo "  Removing .o and .mod files"
	#@\rm *.o *.mod
	@echo "  "
	@echo "####################################################################################"
	@echo "  "
	@echo "  Compilation Successfully Completed"
	@echo "  "
	@echo "  $(COMPILINGMESSAGE1)"
	@echo "  $(COMPILINGMESSAGE2)"
	@echo "  $(COMPILINGMESSAGE3)"
	@echo "  "
	@echo "####################################################################################"
	@echo "  "

%.o: %.f90
	@echo "  Compiling $<"
	@$(FC) $(FOPENMP_FLAG)  $(FFLAGS) -c $<

clean:
	@echo "####################################################################################"
	\rm *.o *.mod *.optrpt


