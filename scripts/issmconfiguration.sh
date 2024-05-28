#/bin/bash
#This script picks up whatever configuration files exists in trunk/configs, 
#and offers the user the choice to reconfigure the ISSM compilation using
#a given configuration file: 

#keep track of present directory: 
presendir=`pwd`

if test -d "$ISSM_DIR/configs" ; then
	cd $ISSM_DIR/configs
	LIST=`ls`
	
	if test -d "$JPL_SVN/usr/$USER/configs"; then
		cd $JPL_SVN/usr/$USER/configs 
		LIST2=`ls`
	fi
	
	#print choices
	COUNT=0;
	printf 'ISSM wide configurations\n'
	for STEP in $LIST
	do
		let COUNT=$COUNT+1
		printf '%3i: %s\n' $COUNT $STEP
	done
	printf 'Personal configuration\n'
	for STEP in $LIST2
	do
		let COUNT=$COUNT+1
		printf '%3i: %s\n' $COUNT $STEP
	done

	echo -n "Configuration choice: "
	read choice 

	#Now go backto the list and retrieve the name of the configuration file: 
	COUNT=0;
	for STEP in $LIST
	do
		let COUNT=$COUNT+1
		if [[ $COUNT == $choice ]]; then
			configurename=$STEP
		fi
	done
	for STEP in $LIST2
	do
		let COUNT=$COUNT+1
		if [[ $COUNT == $choice ]]; then
			configurename=$STEP
		fi
	done

	#Now go ahead and configure: 
	echo ""
	echo "Configuring ISSM with following configs: $configurename"
	echo ""

	cd $ISSM_DIR 

	#at this point, was a cleanup of the archive requested? 
	if [ "$1" == "clean" ]; 
	then 
		make uninstall && make distclean
	fi

	source ./scripts/automakererun.sh 
	if [ -f configs/$configurename ]; then 
		source configs/$configurename
	else 
		source $JPL_SVN/usr/$USER/configs/$configurename
	fi
	
	#we are done, go back to original directory: 
	cd $presendir
else
	echo "Configuration directory does not exist!"
	exit
fi

#alias aut='a=`pwd` && cd $ISSM_DIR && ./scripts/automakererun.sh && ./configs/config-macosx64-larour-nopetsc.sh'
