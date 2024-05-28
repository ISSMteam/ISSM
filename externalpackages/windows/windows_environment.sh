#This file sources all relevant scripts to setup the paths to windows compilers.

#Where are the configuration files for each compiler: 
config_dir="$ISSM_DIR/externalpackages/windows/configs"

#your choise of compiler: 
# 1: sdk 7.1 32 bits on Win7
# 2: sdk 7.1 64 bits on Win7
# 3: intel compiler on Win7
# 4: intel compiler on WinXP

#Determine OS version using uname: 
version=`uname -m | grep x86_64`
if [[ $version == "" ]];then
	compiler=1
else
	compiler=2
fi

#If you want to override and use intel compilers: 
#compiler=3;


#source corresponding environment variables: 

if [[ "$compiler" == "1" ]]; then 
	source $config_dir/sdk10.0-win32.sh
elif [[ "$compiler" == "2" ]]; then 
	source $config_dir/sdk10.0-win64.sh
elif [[ "$compiler" == "3" ]]; then 
	source $config_dir/intel-win7.sh
else 
	source $config_dir/intel-winXP.sh
fi

#finally, out of ISSM_DIR, we need to create an ISSM_DIR_WIN variable for Matlab to pick up on.
ISSM_DIR_WIN=`cygpath -m $ISSM_DIR`
export ISSM_DIR_WIN

#Now source for fortran environment: 
#source $ISSM_DIR/externalpackages/windows/fortran_environment.sh
