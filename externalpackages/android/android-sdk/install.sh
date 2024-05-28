#!/bin/bash

# This installs the Android SDK (Software Development Kit)
# which is needed for the compilation of the Java project. 

source $ANDROID_DIR/android_aux.sh

# Different steps here. 
#   0: do all
#   1: install sdk, ant and sdk tools
#   2: install an emulator.
#   3: test the emulator
#   4: cleanup

present_dir=`pwd`;
sd_card="issm-sdcard"

if [[ $step == "1" || $step == "0" ]]; then

	# Cleanup the install
	rm -rf install

	# Download from ISSM server
	$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/android-sdk_r'$sdk_rev'-macosx.zip' 'android-sdk_r'${sdk_rev}'-macosx.zip'

	# Install Android SDK and NDK.
	unzip -o android-sdk_r${sdk_rev}-macosx.zip

	# Move to install
	mv -f android-sdk-macosx install

	# Post_install configuration: 
	# We need specific settings for specific platforms, for the SDK to 
	# function properly

	# For now, we need to install:  
	# android sdk platform tools  
	# and a specific android api: API 16, API 15 and API 14
	# Note: API 16, API 15 and 14 correspond to Android 4.1, 4.0.3 and 4.0 respectively. 

	cd install/tools/ && source ./android update sdk -t platform-tool,${api_levels},system-image --no-ui

fi

if [[ $step == "2" || $step == "0" ]]; then

	# Once this is done, we need to install an emulator. Location will default to ~/.android/avd, 
	# which we will move to $ISSM_DIR/externalpackages/android-emulators.  
	# For now, it's called: Android-4.0.3

	# Here we delete the Android-4.0.3 device if it already exists.
	cd $present_dir/install/tools

    if [ -e $ANDROID_DIR/android-emulators/$default_droid ] 
    then
        echo "Deleting previously created device: $default_droid"
	    ./android delete avd -n $default_droid
    fi

	# Android will prompt the user to specify hardware emulation options. For now, default
	# default settings will suffice. Press 'enter' to take default settings or enter 'no'.

	./android create avd -f -n $default_droid -t 1 -p $ANDROID_DIR/android-emulators/$default_droid --abi armeabi-v7a
    echo "Creating an SD Card"
    ./mksdcard -l $sd_card 2G $ANDROID_DIR/android-emulators/$sd_card.img
fi

if [[ $step == "3" || $step == "0" ]]; then
    # Here we will start up our default emulator to test that it is working properly.
    # Once the device has booted we will use the Android Debug Bridge tool to gain
    # a terminal in our device.

	cd $present_dir/install/tools
	./emulator -avd $default_droid -sdcard $ANDROID_DIR/android-emulators/$sd_card.img &

    cd ../platform-tools
    ./adb wait-for-device shell
fi

if [[ $step == "4" || $step == "0" ]]; then
	rm -rf install
fi
