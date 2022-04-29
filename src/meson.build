project('bpp', ['cpp','c'],
        version: '4.4.1',
        default_options : [
          'optimization=3',     # -O3
          'warning_level=1',    # -Wall
          'buildtype=release',  # specify -Ddebug=true on commmand line for -g
          'b_ndebug=if-release'
        ],
        meson_version: '>= 0.56',
        license: 'AGPLv3')

cc = meson.get_compiler('c')

warn_flags = ['-Wsign-compare']
add_project_arguments(['-D_GNU_SOURCE'] + cc.get_supported_arguments(warn_flags), language: 'c')

avx_sources  = ['core_partials_avx.c' , 'core_likelihood_avx.c' ]
avx2_sources = ['core_partials_avx2.c', 'core_likelihood_avx2.c']

if host_machine.cpu_family().startswith('x86')
  if get_option('have_sse3')
    add_project_arguments(['-DHAVE_SSE3','-msse3'], language: 'c')
  endif

  if get_option('have_avx')
    add_project_arguments(['-DHAVE_AVX'], language: 'c')
  endif

  if get_option('have_avx2')
    add_project_arguments(['-DHAVE_AVX2'], language: 'c')
  endif

  if get_option('have_avx')
    libavx  = static_library('avx' , avx_sources,  c_args: ['-mavx'])
  else
    libavx = dependency('', required: false)
  endif

  if get_option('have_avx2')
    libavx2 = static_library('avx2', avx2_sources, c_args: ['-mavx2','-mfma'])
  else
    libavx2 = dependency('', required: false)
  endif
endif

sources = [
  'bpp.c',
  'rtree.c',
  'util.c',
  'arch.c',
  'phylip.c',
  'msa.c',
  'maps.c',
  'locus.c',
  'mapping.c',
  'compress.c',
  'hash.c',
  'list.c',
  'stree.c',
  'random.c',
  'gtree.c',
  'core_partials.c',
  'core_pmatrix.c',
  'core_likelihood.c',
  'output.c',
  'core_partials_sse.c',
  'dlist.c',
  'allfixed.c',
  'core_likelihood_sse.c',
  'prop_mixing.c',
  'method.c',
  'delimit.c',
  'prop_rj.c',
  'summary.c',
  'cfile.c',
  'hardware.c',
  'revolutionary.c',
  'diploid.c',
  'dump.c',
  'load.c',
  'summary11.c',
  'simulate.c',
  'cfile_sim.c',
  'gamma.c',
  'prop_gamma.c',
  'threads.c',
  'treeparse.c',
  'parsemap.c',
  'msci_gen.c',
  'constraint.c',
  'debug.c',
  'lswitch.c',
  'ming2.c'
]

# On some platforms there isn't a separate libm
libm = cc.find_library('m', required: false)

threads = dependency('threads')

bpp = executable('bpp',
                 sources,
                 include_directories: include_directories('.'),
                 dependencies: [libm, threads],
                 link_with: [libavx, libavx2],
                 install: true)

#### Summary

summary({'host': host_machine.system()
        }, section: 'Architecture')

summary({'optimization': get_option('optimization'),
         'debug': get_option('debug')
        },section: 'Configuration')


