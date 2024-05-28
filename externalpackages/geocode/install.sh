#/bin/bash
unzip geoCode.zip
sed -i -e 's/http:/https:/g' geoCode.m
