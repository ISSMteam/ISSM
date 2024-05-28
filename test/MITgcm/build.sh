#!/bin/bash
#This script compiles and links MITgcm

#recover hostname and model path:
hostname="$1"
modelpath="$2"

if [ -e ~/.bashrc ]; then
    source ~/.bashrc
fi

# Get MITgcm code, if needed
if [ ! -d "$modelpath/../MITgcm/install" ]; then
    cd $modelpath/../MITgcm
    source install.sh
    cd $modelpath
fi

# Create build directory, if needed
cd $modelpath
if [ ! -d "build" ]; then mkdir build; fi
cd build

#create MITgcm makefile for this run, if needed
if [ ! -f Makefile ]; then
	case $hostname in
		"pleiades")
			$modelpath/../MITgcm/install/tools/genmake2 -of $SLR_DIR/models/ice-ocean/configs/linux_amd64_gfortran+mpi_ice_nas -mo ../code -rd $modelpath/../MITgcm/install
			;;
		"babylon")
			$modelpath/../MITgcm/install/tools/genmake2 -of $modelpath/../MITgcm/install/tools/build_options/linux_amd64_ifort -mpi -mo $modelpath/../MITgcm/code -rd $modelpath/../MITgcm/install
			export LD_LIBRARY_PATH="$ISSM_DIR/externalpackages/petsc/install/lib:/dartfs-hpc/admin/opt/el7/intel/compilers_and_libraries_2019.3.199/linux/compiler/lib/intel64:$ISSM_DIR/externalpackages/triangle/install/lib"
			;;
		"amundsen")
			export LDADD="-L$ISSM_DIR/externalpackages/petsc/install/lib -lmpi -lmpifort"
			$modelpath/../MITgcm/install/tools/genmake2 -mpi -mo $modelpath/../MITgcm/code -rd $modelpath/../MITgcm/install
			;;
		*)
			$modelpath/../MITgcm/install/tools/genmake2 -mpi -mo $modelpath/../MITgcm/code -rd $modelpath/../MITgcm/install
			;;
	esac
fi

#create MITgcm code links for this run, if needed
if [ ! -f BUILD_INFO.h ]; then
    make depend
fi

#run make command
STR=`uname -v`
SUB='ARM64'
if [[ "$STR" == *"$SUB"* ]]; then
    arch -arm64 make -j &> Makefile.log
else
    make -j 4 &> Makefile.log   
fi
