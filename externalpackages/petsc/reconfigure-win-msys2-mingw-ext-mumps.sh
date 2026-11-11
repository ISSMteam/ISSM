#!/bin/bash
set -eu


## Constants
#
PETSC_ARCH_RECONF="arch-mswin-reconf"
PETSC_DIR="${ISSM_DIR}/externalpackages/petsc/src" # DO NOT CHANGE THIS
PREFIX="${ISSM_DIR}/externalpackages/petsc/install" # Set to location where external package should be installed (UNIX)

# Reconfigure with external MUMPS installation.
#
cd ${PETSC_DIR}
echo "Reconfiguring with external MUMPS installation..."
${PETSC_ARCH_RECONF}/lib/petsc/conf/reconfigure-${PETSC_ARCH_RECONF}.py \
    --with-shared-libraries=0 \
    --with-debugging=0 \
    --with-valgrind=0 \
    --with-x=0 \
    --with-ssl=0 \
    --with-cc=/ucrt64/bin/mpicc.exe \
    --with-cxx=/ucrt64/bin/mpicxx.exe \
    --with-fc=/ucrt64/bin/mpifort.exe \
    --with-mpiexec="${MPIEXEC_DIR}/mpiexec.exe" \
    --with-metis-lib="-L${PREFIX}/lib -lmetis" \
    --with-metis-include="${PREFIX}/include" \
    --with-parmetis-lib="-L${PREFIX}/lib -lparmetis" \
    --with-parmetis-include="${PREFIX}/include" \
    --with-scalapack-lib="-L${PREFIX}/lib -lscalapack" \
    --with-mumps-include="${MUMPS_ROOT}/include" \
    --with-mumps-lib="-L${MUMPS_ROOT}/lib -ldmumps -lsmumps -lcmumps -lzmumps -lmumps_common -lpord"

# Compile and install
make PETSC_DIR="${PETSC_DIR}" PETSC_ARCH="${PETSC_ARCH_RECONF}" all
make PETSC_DIR="${PETSC_DIR}" PETSC_ARCH="${PETSC_ARCH_RECONF}" install
