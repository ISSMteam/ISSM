#!/bin/bash
set -eu

#Some cleanup
rm -rf src install ffmpeg-1.1.2
mkdir src install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/ffmpeg-1.1.2.tar.gz' 'ffmpeg-1.1.2.tar.gz'

#Untar 
tar -zxvf  ffmpeg-1.1.2.tar.gz

#Move ffmpeg into src directory
mv ffmpeg-1.1.2/* src
rm -rf ffmpeg-1.1.2

#Configure ffmpeg
cd src

export CFLAGS=" -arch x86_64"

./configure --prefix="$ISSM_DIR/externalpackages/ffmpeg/install" --disable-yasm

#Compile ffmpeg
if [ $# -eq 0 ]; then
	make
else
	make -j $1
fi
make install 
