project('LTRANS_Zlev',
        ['fortran'],
        default_options : [
                          'buildtype=release', # {plain, debug, debugoptimized, release, minsize, custom}
                          'warning_level=0',  # From 0 = none to 3 = highesti
                          ]
        )

_compiler = meson.get_compiler('fortran')
#######################################################################
#----------------- BEGINNING OF THE CUSTOMISABLE AREA -----------------
#######################################################################
message('compiler is '+ _compiler.get_id())
#------------------- eventually setup optimization arguments -------------
_bt = get_option('buildtype')
if _bt == 'debugoptimized' or _bt == 'release'
   if _compiler.get_id() == 'intel'
    add_global_arguments([  '-funroll-loops',
                            '-assume','byterecl',
                            '-xcore-avx2',
                            '-axCORE-AVX2,AVX'
                          ],language : ['fortran'] )
  elif _compiler.get_id() == 'gcc'
    add_global_arguments([  '-march=native',
                            '-ffast-math',
                            '-ftree-vectorize',
                            '-ftree-vectorizer-verbose=0',
                            '-funroll-loops',
                          ],language : ['fortran'] )
  endif
endif
if _bt == 'debug'
  if _compiler.get_id() == 'intel'
    add_global_arguments([ 
                          '-trace',
                          '-qopt-report=5',
                          '-qopt-report-phase=all',
                          '-diag-disable 10397', #disables info on qoptr files creation
                          '-traceback',
                          '-fp-stack-check',
                          '-fpe0',
                          '-ftrapuv',
                          '-CB',#'-fp-model source','-debug all',
                          ],language : ['fortran'] )
  else
    add_global_arguments([
                          '-fbacktrace','-fbounds-check',
                          '-Wall','-Wextra',
                          '-fimplicit-none',
                          '-fcheck=all',
                          '-Wtabs',
                          '-ffpe-trap=zero,invalid,overflow',
                          ],language : ['fortran'] )
  endif

endif

#######################################################################
#-------------------- END OF THE CUSTOMISABLE AREA --------------------
#######################################################################

#--------------- check eventual incompatibilities ---------------------------
if get_option('use_OMP')==true
   message('Compiling with OpenMP flags')
   #omp = dependency('openmp')
   if _compiler.get_id() == 'intel'
      add_global_arguments(['-qopenmp','-fpp','-auto'],language : ['fortran'] )
   elif _compiler.get_id() == 'gcc'
      add_global_arguments(['-fopenmp'],language : ['fortran'] )
   endif
endif

###############################################################################
message('------------------- preparing global settings --------------------')
###############################################################################

#if get_option('Verbose')
#  warning('Verbose option is activated, it induces verbose LTRANS output slowing down A LOT the LTRANS program')
#endif

#--------------- setting up pragma precompilation directives ----------------
_config_data = configuration_data()

_config_data.set('use_OMP',get_option('use_OMP'))

if _compiler.get_id() == 'intel'
    add_global_arguments('-fpp',language : ['fortran'] )
elif _compiler.get_id() == 'gcc'
  add_global_arguments('-cpp',language : ['fortran'] )
endif

###############################################################################
message('------------------- preparing LTRANS_Zlev build target ----------------')
###############################################################################
_LTRANS_dependencies= []
_LTRANS_include_directories=[]
_LTRANS_source_files=[]
_LTRANS_fortran_args = []

# create list of dependencies, containing for each one the array : [library_name, installation_path, compiler]
_lib_properties=[['netcdf','NETCDF_C_HOME_PATH',_compiler] ] 
_lib_properties +=[['netcdff','NETCDF_F_HOME_PATH',_compiler]] 
if get_option('netcdf_compiled_with_hdf5') == true
   _lib_properties +=[['hdf5_hl','HDF5_HOME_PATH',_compiler]]      
   _lib_properties +=[['hdf5','HDF5_HOME_PATH',_compiler]]      
   _lib_properties +=[['curl','CURL_HOME_PATH',_compiler]]      
endif
 

foreach _libproperty : _lib_properties
    _lib=_libproperty[0] 
    _opt=_libproperty[1]
    _dir=get_option(_opt)
    _comp=_libproperty[2] 
    if _dir != 'none'
       message('Using user-provided '+_lib+' directory '+_dir)
       _LTRANS_dependencies += [ _compiler.find_library(_lib,dirs:[ _dir+'/lib/', _dir+'/lib64/' ] )  ]
       _LTRANS_include_directories += [ include_directories(_dir+'/include/') ]
    else
       _LTRANS_dependencies += [dependency(_lib,
                                 not_found_message:_lib+' dependency couldnt be found using pkg-config,'
                                      +' meson cannot proceed because '
                                      +' both library and include paths are required to use this libray\.' 
                                      +' Please provide installation path through command line '
                                      +' -D'+_opt+'=PATH_TO_INSTALLATION_DIRECTORY'
                                      +' or add file .pc path to environmental variable PKG_CONFIG_PATH'
                                      +' (path might be $<LIBNAME>_HOME/lib/pkgconfig)',
                                 )]
    endif
endforeach

#------------------- add common source files ---------------------------
_LTRANS_source_files += [
                         'LTRANS_src/ver_turb_module.f90',
                         'LTRANS_src/tension_module.f90',
                         'LTRANS_src/settlement_module.f90',
                         'LTRANS_src/random_module.f90',
                         'LTRANS_src/point_in_polygon_module.f90',
                         'LTRANS_src/norm_module.f90',
                         'LTRANS_src/interpolation_module.f90',
                         'LTRANS_src/hor_turb_module.f90',
                         'LTRANS_src/gridcell_module.f90',
                         'LTRANS_src/conversion_module.f90',
                         'LTRANS_src/boundary_module.f90',
                         'LTRANS_src/parameter_module.f90',
                         'LTRANS_src/oil_module_Zlev.f90',
                         'LTRANS_src/LTRANS.f90',
                         'LTRANS_src/hydrodynamic_module.f90',
                         'LTRANS_src/behavior_module.f90',
                        ]


_LTRANS_exe=executable('LTRANS_Zlev',
           _LTRANS_source_files ,
           fortran_args : _LTRANS_fortran_args,
           include_directories : _LTRANS_include_directories,
           dependencies : _LTRANS_dependencies,
           )
subdir('SIM/tests')
