_prefix = meson.current_source_dir()

message('- - - - - [ meson configuring running tests ] - - - - - - - - - - ')

_get_MITgcm_outputs = find_program('get_MITgcm_outputs.sh')
test('--------------- get MITgcm outputs from zenodo -----------------------------------',
          _get_MITgcm_outputs,
          workdir:meson.current_source_dir()+'/..',
          timeout:1800,
          priority:100)

_synch = find_program('synch.sh')
test('--------------- barrier insuring MITgcm_output files are downloaded --------------',
    _synch,
    args:[],
    workdir:meson.current_build_dir(),
    priority:99, 
    is_parallel : false)

_priority=98

_run_and_check_success = find_program('run_and_check_success.sh')
_threshold='1E-05'
#_threshold='none'
test('--------------- running LTRANS ---------------------------------------------------',
          _run_and_check_success,
          args:[
                _LTRANS_exe.full_path(),
                'LTRANS_boxes_NiNj64.data',
                'boxes_NiNj64-endfile.ref',
                'boxes_NiNj64-endfile.csv',
                _prefix+'/',
                _threshold,
               ],
          workdir:meson.current_build_dir(),
          priority:_priority,
          timeout:180,
          is_parallel : true)
