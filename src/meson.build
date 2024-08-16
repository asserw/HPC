lib_sources = [
    'hello.cpp',
    'verlet.cpp',
    'xyz.cpp',
    'lj_direct_summation.cpp',
    'compute_kinetic_energy.cpp',
    'neighbors.cpp',
    'lj_direct_summation_neighbors.cpp',
    'berendsen_thermostat.cpp',
]

lib = static_library(
    'lib',  # Unique name
    lib_sources,
    dependencies : [eigen, mpi]
)

lib_incdirs = include_directories('.')
