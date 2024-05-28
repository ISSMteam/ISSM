#!/bin/bash
#build windows archive of binaries.

#Some local script functions 
function today_date {
suffix=`date | awk '{printf("%s-%s-%s",$2,$3,$6);}'` 
echo $suffix;
}

#Create tar file, with today's date in the title;
today=`today_date`

cd $ISSM_DIR/bin

#Filter out .svn files
rm -rf list
ls *.mexw32 | grep -v "\.svn" > list;

tar zcvf ../issm-1.0-win-$today.tar.gz  `cat list`
rm -rf list

