#!/bin/bash
#
# android_aux.sh serves as an auxiliary script for all installation
# scripts with the Android suffix.
#
# TODO: include M4 macros for generic options.
#

sdk_rev=20
ndk_rev=8
api_levels="android-14,android-15,android-16"
host_triplet="arm-linux-androideabi"
default_droid="android-4.1"

step=0;
j=1;

echo ""
echo "This install script utilizes 'android_aux.sh' to allow for options."
echo "For usage information enter: '--help'"
echo Number of arguments is: $#

for arg in $* 
do 
    if [[ "$arg" =~ --step=([0-9])* ]]; then
        step=${BASH_REMATCH[1]}; 
        echo "Setting step to: " $step
    elif [[ "$arg" == "--help" ]]; then
        echo ""
        echo "USAGE: $ install.sh [--step=#] [-j#]"
        echo ""
        echo "Where '#' is some integer."
        echo "To check the number of steps check the install script."
        echo ""
        exit 1;
    elif [[ "$arg" =~ -j=([1-9]+[0-9]*) ]]; then
        j=${BASH_REMATCH[1]}; 
        echo "Number of jobs set to: " $j
    else
        echo "Option not recognized"
        exit 1;
    fi
done
