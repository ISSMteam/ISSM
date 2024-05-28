Instructions for getting MITgcm:   http://mitgcm.org/ 


(1)
add to your .bashrc file
# CVS 
export CVSROOT=':pserver:cvsanon@mitgcm.org:/u/gcmpack'

source .bashrc

cvs login 
 ( enter the CVS password: "cvsanon" )

cvs co -P MITgcm

(
or if you don't want the verification packages
cvs co -P MITgcm_code
)


(2)
and if you don't like CVS  go to :

http://mitgcm.org/download/

and download the most recent checkpoint:
e.g.
     MITgcm_c62w.tar.gz





==========================================================
Instructions for generating and running a 1-CPU experiment
==========================================================

  cd MITgcm/verification/lab_sea
  cd build
  cp ../code/*.h ../code/packages.conf .
  ../../../tools/genmake2
  make depend
  make
  cd ../input
  ../build/mitgcmuv > output.txt

Use matlab script to look at the output
  cd ../../../verification/lab_sea/matlab
  matlab
  lookat_ice  (you might have to modify the script)



================================================================
Instructions for running the "weddell" 200x160x50 configuration
================================================================
face=6; ix=101:300; jx=290:449; kx=1:50;

1. Obtain copies of following directories:
 ftp://ecco2.jpl.nasa.gov/data1/weddell/code
 ftp://ecco2.jpl.nasa.gov/data1/weddell/run_template
 ftp://ecco2.jpl.nasa.gov/data1/data/era40/era40_ecmwf_blend
 ftp://ecco2.jpl.nasa.gov/data1/data/blend_forcing/cube59_GPCP

2. Get and compile code:
 cvs co MITgcm_code
 cd MITgcm
 mkdir bin exe
 cd bin
 ../tools/genmake2 -mods=../../code
 make depend
 make -j

3. Model execution:
 cd ../exe
 cp ../../run_template/* .
 cp ../bin/mitgcmuv .
 ./mitgcmuv >& output.txt &





