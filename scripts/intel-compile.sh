#!/bin/bash
#Why don't we just type make? (shouldn't automake have taken care of this?)
#The problem is the /Fe option from the intel compiler, which we weren't able to 
#get automake to recognize. End result is that every file compiled is not named libISSM_a-name, 
#but just name.  This makes the creation of libISSM.a impossible, as none of its objects 
#can be found with the correct name. 
#As a fix, we rename the objects, and then link.

#First compile.
#make

#Then change the names
list=`ls *.obj | grep -v libISSM_a`
for i in `echo $list`
do
	mv $i libISSM_a-$i
done

#Now create the library out the .obj files
rm -rf libISSM.a
lib.exe /out:libISSM.a *.obj
