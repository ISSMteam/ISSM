#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#if defined(_HAVE_ADJOINTMPI_)
#include "externals/ampi_interface_realreverse.cpp"
#elif defined(_HAVE_MEDIPACK_)
#include "medi/medi.cpp"
#else
#error "Cannot compile without MeDiPack and AdjointMpi"
#endif
