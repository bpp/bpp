project('bpp', 'c',
        version: '4.4.1',
        default_options : [
          'optimization=3',     # -O3
          'warning_level=1',    # -Wall
          'buildtype=release',  # specify -Ddebug=true on commmand line for -g
          'b_ndebug=if-release'
        ],
        meson_version: '>= 0.56',
        license: 'AGPLv3')

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
  'dlist.c',
  'allfixed.c',
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

cc = meson.get_compiler('c')

warn_flags = ['-Wsign-compare']
add_project_arguments(['-D_GNU_SOURCE','-DSIMDCONFIG'] + cc.get_supported_arguments(warn_flags), language: 'c')

cdata = configuration_data()

simd = import('unstable-simd')
sse3_avx = simd.check('sse3_avx',
                      sse3: ['core_partials_sse.c','core_likelihood_sse.c'],
                      avx:  ['core_partials_avx.c','core_likelihood_avx.c'],
                      compiler: cc)

avx2 = simd.check('avx', avx2: ['core_partials_avx2.c','core_likelihood_avx2.c'], compiler: cc, c_args: ['-mfma'])
simdlibs = [sse3_avx[0], avx2[0]]
cdata.merge_from(avx2[1])
cdata.merge_from(sse3_avx[1])

configure_file(output : 'simdconfig.h', configuration : cdata)

# On some platforms there isn't a separate libm
libm = cc.find_library('m', required: false)

threads = dependency('threads')

bpp = executable('bpp',
                 sources,
                 include_directories: include_directories('.'),
                 dependencies: [libm, threads],
                 link_with: simdlibs,
                 install: true)

#### Summary

summary({'host': host_machine.system()
        }, section: 'Architecture')

summary({'optimization': get_option('optimization'),
         'debug': get_option('debug')
        },section: 'Configuration')


