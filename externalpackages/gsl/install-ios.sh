#!/bin/bash
set -eu

source $ANDROID_DIR/android_aux.sh
export DEVROOT="/Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer/"
export SDKROOT="$DEVROOT/SDKs/iPhoneOS6.0.sdk/"
export CC="$DEVROOT/usr/bin/gcc -arch armv7"
export CFLAGS="-O -arch armv7 -isysroot $SDKROOT"
export CXXFLAGS="-O -arch armv7 -isysroot $SDKROOT"
export CPP=$DEVROOT/usr/bin/cpp
export CXX=$DEVROOT/usr/bin/g++
export AR=$DEVROOT/usr/bin/ar
export AS=$DEVROOT/usr/bin/as
export LIBTOOL=$DEVROOT/usr/bin/libtool
export LDFLAGS="-lstdc++ -arch armv7 -isysroot $SDKROOT"
export STRIP=$DEVROOT/usr/bin/strip
export RANLIB=$DEVROOT/usr/bin/ranlib


if [[ $step == "1" || $step == "0" ]]; then

    #Some cleanup
    rm -rf src install gsl-1.15
    mkdir src install

    #Download from ISSM server
    $ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/gsl-1.' 'gsl-1.15.tar.gz'

    #Untar 
    tar -zxvf  gsl-1.15.tar.gz

    #Move gsl into src directory
    mv gsl-1.15/* src
    rm -rf gsl-1.15
fi

#Configure gsl
if [[ $step == "2" || $step == "0" ]]; then
    cd src

    patch Makefile.am < ./../Makefile.am.patch

    autoreconf -if

    ./configure \
		--build="i386-apple-darwin10.8.0" \
		--host="arm-apple-darwin10" \
	    --prefix="$ISSM_DIR/externalpackages/gsl/install"
fi

#Compile gsl
if [[ $step == "3" || $step == "0" ]]; then
	cd $ISSM_DIR/externalpackages/gsl/src

    if [ $# -eq 0 ]; then
	    make 
    else
	    make -j $j 
    fi

    make install
fi

