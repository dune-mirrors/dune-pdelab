#ifndef DUNE_PDELAB_COMMON_CHECKS_HH
#define DUNE_PDELAB_COMMON_CHECKS_HH

#ifndef DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY
#ifndef NDEBUG
#define DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY 1
#else
#define DUNE_PDELAB_ENABLE_CHECK_ASSEMBLY 0
#endif
#endif

#endif // DUNE_PDELAB_COMMON_CHECKS_HH
