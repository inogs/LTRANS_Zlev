_prefix = meson.current_source_dir()

message('- - - - - [ meson configuring running tests ] - - - - - - - - - - ')

#_get_MITgcm_outputs = find_program('get_MITgcm_outputs.sh')
#test('--------------- get MITgcm outputs from zenodo -----------------------------------',
#          _get_MITgcm_outputs,
#          workdir:meson.current_source_dir()+'/..',
#          timeout:1800,
#          priority:100)

python_exe = find_program('python3', 'python')
params = files('check_python_libs.py')
r = run_command(python_exe, params,  check: true)
if r.returncode() != 0
  error('could not run python with the required python module libraries, next test will fail, please set up an adeguate python environment and try again.')
endif


_create_fake_input_files = find_program('create_fake_input_files.sh')
test('--------------- create fake input files (grid and fields)--------------------------',
       _create_fake_input_files, 
       args:[
             _prefix+'/'
            ],
       workdir:meson.current_build_dir(),
       timeout:300,
       priority:100)

_synch = find_program('synch.sh')
test('--------------- barrier insuring that the input files are present -----------------------',
    _synch,
    args:[],
    workdir:meson.current_build_dir(),
    priority:99, 
    is_parallel : false)

_priority=98

_run_and_check_success = find_program('run_and_check_success.sh')
_threshold='1E-05'
#_threshold='none'

_SorZlevels=['Slevels','Zlevels']
_SubCases=['Cur','Win','CurWin',
           'Cur_B998','CurWin_B998',
           'CurWin_StkfromWind_B998',
           'Win_StkfromWind',
           'CurWin_StkfromWind',
           'StkRead','StkRead_B998',
           'Cur_StkRead','Win_StkRead','CurWin_StkRead',
           'CurWin_StkRead_B998',]
foreach _leveltype : _SorZlevels
  foreach _subcase : _SubCases
    test('--------------- running LTRANS_'+_leveltype+'_'+_subcase+' case -----------------------------------------',
          _run_and_check_success,
          args:[
                _LTRANS_exe.full_path(),
                'LTRANS_'+_leveltype+'_'+_subcase+'.data',
                'solutions/LTRANS_'+_leveltype+'_'+_subcase+'-endfile.csv',
                'LTRANS_'+_leveltype+'_'+_subcase+'-endfile.csv',
                'Adjacentelements-'+_leveltype+'.data',
                _prefix+'/',
                _threshold,
               ],
          workdir:meson.current_build_dir(),
          priority:_priority,
          timeout:300,
          is_parallel : true)
  endforeach
endforeach
