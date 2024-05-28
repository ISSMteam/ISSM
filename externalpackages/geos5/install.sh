#!/bin/bash
set -eu

#Some cleanup
rm -rf  GEOSagcm
#svn download
svn --username eric.larour@jpl.nasa.gov checkout http://geos5.org/svn/branches/Fortuna-2_5_p1 GEOSagcm
