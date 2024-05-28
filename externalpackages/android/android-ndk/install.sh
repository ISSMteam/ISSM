#!/bin/bash

# This installs the Android NDK (Native Development Kit)
# which is needed for the compilation of C/C++ code into the 
# native architecture of the Android device.

# The android-ndk-r8-darwin-x86.tar.bz2 ndk.tar.bz2 file was downloaded 
# from the android developer website: 
source $ANDROID_DIR/android_aux.sh

if [[ $step == "1" || $step == "0" ]]; then
    # Cleanup the install
    rm -rf install
    rm -rf $ANDROID_DIR/arm-linux-android-${default_droid}

    # Download from ISSM server
    $ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/android-ndk-r'${ndk_rev}'-darwin-x86.tar.bz2' 'ndk.tar.bz2'

    # Untar
    tar -zxvf ndk.tar.bz2
    rm -rf ndk.tar ndk.tar.bz2

    # Move to install
    mv android-ndk-r${ndk_rev} install
fi

# Create Standalone Development Directory
# Note: API-14 corresponds to Android 4.0.
if [[ $step == "2" || $step == "0" ]]; then
    $ANDROID_DIR/android-ndk/install/build/tools/make-standalone-toolchain.sh --platform=android-14 --install-dir=$ANDROID_DIR/android-ndk/install/arm-linux-android-install
fi
