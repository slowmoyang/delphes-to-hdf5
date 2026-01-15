project(
  'delphes-to-hdf5',
  'cpp',
  version : '1.0.0',
  default_options : ['warning_level=3']
)

###############################################################################
# dependencies
###############################################################################

# root
root_config = find_program('root-config', required: true)
root_compile_args = run_command(root_config, '--cflags', check: true).stdout().strip().split(' ')
root_link_args = run_command(root_config, '--ldflags', '--libs', check: true).stdout().strip().split(' ')
root_incdir = run_command(root_config, '--incdir', check: true).stdout().strip().split(' ')
root_dep = declare_dependency(
    compile_args: root_compile_args,
    link_args: root_link_args + ['-lEG'],
    include_directories: root_incdir
)

# delphes
delphes_config = find_program('delphes-config.py', required: true)
delphes_compile_args = run_command(delphes_config, '--cflags', check: true).stdout().strip().split(' ')
delphes_link_args = run_command(delphes_config, '--ldflags', '--libs', check: true).stdout().strip().split(' ')
delphes_incdir = run_command(delphes_config, '--incdir', check: true).stdout().strip().split(' ')
delphes_dep = declare_dependency(
  compile_args: delphes_compile_args,
  link_args: delphes_link_args,
  include_directories: delphes_incdir,
)

# dependencies installed by wrap
tomlplusplus_dep = dependency(
  'tomlplusplus',
  required: true,
  static: true,
)

argparse_dep = dependency(
  'argparse',
  required: true,
)

indicators_dep = dependency(
  'indicators',
  required: true,
)

spdlog_dep = dependency(
  'spdlog',
  required: true,
)

highfive_dep = dependency(
  'HighFive',
  required: true,
)

# dependencies
dependencies = [
  root_dep,
  delphes_dep,
  argparse_dep,
  indicators_dep,
  tomlplusplus_dep,
  spdlog_dep,
  highfive_dep,
]

###############################################################################
# hi
###############################################################################

executable(
  'delphes-to-hdf5',
  './src/delphes-to-hdf5.cc',
  dependencies: dependencies,
  install : true,
)
#
# executable(
#   'test-highfive',
#   './src/test-highfive.cc',
#   dependencies: [highfive_dep],
#   install: true,
# )
