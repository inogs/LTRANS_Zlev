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
test('--------------- running LTRANS with Zlevels case -----------------------------------------',
          _run_and_check_success,
          args:[
                _LTRANS_exe.full_path(),
                'LTRANS_fake_Zlevels.data',
                'fake_Zlevels-endfile.ref',
                'fake_Zlevels-endfile.csv',
                'Adjacentelements-Zlevels.data',
                _prefix+'/',
                _threshold,
               ],
          workdir:meson.current_build_dir(),
          priority:_priority,
          timeout:300,
          is_parallel : true)
test('--------------- running LTRANS with Slevels case -----------------------------------------',
          _run_and_check_success,
          args:[
                _LTRANS_exe.full_path(),
                'LTRANS_fake_Slevels.data',
                'fake_Slevels-endfile.ref',
                'fake_Slevels-endfile.csv',
                'Adjacentelements-Slevels.data',
                _prefix+'/',
                _threshold,
               ],
          workdir:meson.current_build_dir(),
          priority:_priority,
          timeout:300,
          is_parallel : true)
