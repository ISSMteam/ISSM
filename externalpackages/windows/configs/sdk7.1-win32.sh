# -----------------------------------------------------------------------------
# 
# The following variables can vary from system to system and standard
# installation is assumed throughout. If different directories were used when
# installing MSVC, or Win SDK.
#
# TODO: Bring out the host machine arhictecture specific stuff from INCLUDE,
# LIB and LIBPATH.
#
# -----------------------------------------------------------------------------

# The version of Visual Studio is 10.0. Newer versions should work as well.
export MSVC_DIR='C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\'

# SDK change from 7 to 8 involved changing the locations of important 
# libraries. If you wish to use 8.0 or 8.1 then you need to provide the 
# directory of 'Windows Kits' instead.
export WIN_SDK_DIR='C:\Program Files (x86)\Microsoft SDKs\Windows\v7.1A\'

# Unfotunately, standard installation directories of Matlab usually include 
# white spaces that are not handled correctly by the command interpreter. As
# such, the directory where one would find Matlab headers and libraries might
# need to be provided as environment variables.
export MATLAB_DIR_WIN=`cygpath -w ${MATLAB_DIR}`
export MATLAB_DIR_LIB="${MATLAB_DIR_WIN}\\extern\\lib\\microsoft"

# Information about the .NET framework is required to run the MSVC toolchain
export FrameworkDir='C:\Windows\Microsoft.NET\Framework\'
export FrameworkVersion=v4.0.30319

# LIB and LIBPATH seem redundant, but MSVC linker and compiler use different 
# variables for the same purpose.
export INCLUDE="${MSVC_DIR}include;${WIN_SDK_DIR}Include;${MATLAB_DIR_WIN}\\extern\\include;"
export LIB="${MSVC_DIR}lib;${WIN_SDK_DIR}Lib;${MATLAB_DIR_LIB}"
export LIBPATH="${FrameworkDir}${FrameworkVersion};${MATLAB_DIR_LIB};"
export LIBPATH="${LIBPATH}${MSVC_DIR}lib;${MSVC_DIR}bin;${WIN_SDK_DIR}Lib;"

export MSVC_DIR_UNIX=`cygpath -u "${MSVC_DIR}"`
export PATH="${MSVC_DIR_UNIX}/bin:$PATH"
