#!/bin/bash
set -eu

#Do these commands once:
# cvs -d:pserver:cvsanon@mitgcm.org:/u/gcmpack login
# ( enter the CVS password: "cvsanon" )

#Some cleanup
rm -rf install   bin exe

#add cvs repository
export CVSROOT=':pserver:cvsanon@mitgcm.org:/u/gcmpack'

echo loging into MITgcm cvs: provide password, which is cvsanon
cvs login

#Download code from server
cvs login
cvs co -P MITgcm

#move
mv MITgcm install


#compile code
cd install
mkdir bin exe
cd bin
../tools/genmake2 -mods=../../code
make depend
make -j 8



