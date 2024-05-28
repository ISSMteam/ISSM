/*!\file:  DescriptorIndex: return type of qmu variable: indexed, scaled, nodal or regular
 * + figure out the descriptor root. 
 * Ex: scaled_Thickness_1 should return SCALEDENUM, fill root with Thickness, and initialize index 
 * to 1.
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "../Enum/Enum.h"
#include "../io/io.h"
#include "../Exceptions/exceptions.h"

int  DescriptorIndex(char* root, int* pindex,char* descriptor){ //We assume root has already been allocated, and we just have to copy into it.

	char * pch=NULL;

	/*retrieve first token, separated by underscore: */
	pch = strtok (descriptor,"_");
	if(!pch)_error_("descriptor " << descriptor << " is not correctly formatted!");

	if (strncmp(pch,"scaled",6)==0){
		/*we have a scaled variable. recover the root: */
		pch = strtok (NULL, "_");
		if(!pch)_error_("scaled descriptor " << descriptor << " is not correctly formatted!");
		memcpy(root,pch,(strlen(pch)+1)*sizeof(char));

		/*now recover  the index if it exists: */
		pch = strtok (NULL, "_");
		if(!pch){
			*pindex=-1;
		}
		else{
			sscanf(pch,"%i",pindex);
		}
		return ScaledEnum;
	}
	else if (strncmp(pch,"indexed",7)==0){
		/*we have an indexed variable. recover the root: */
		pch = strtok (NULL, "_");
		if(!pch)_error_("indexed descriptor " << descriptor << " is not correctly formatted!");
		memcpy(root,pch,(strlen(pch)+1)*sizeof(char));
		/*now recover  the index: */
		pch = strtok (NULL, "_");
		if(!pch)_error_("indexed descriptor " << descriptor << " is not correctly formatted!");
		sscanf(pch,"%i",pindex);
		return IndexedEnum;
	}
	else if (strncmp(pch,"nodal",5)==0){
		/*we have an indexed variable. recover the root: */
		pch = strtok (NULL, "_");
		if(!pch)_error_("nodal descriptor " << descriptor << " is not correctly formatted!");
		memcpy(root,pch,(strlen(pch)+1)*sizeof(char));
		/*now recover  the index: */
		pch = strtok (NULL, "_");
		if(!pch)_error_("nodal descriptor " << descriptor << " is not correctly formatted!");
		sscanf(pch,"%i",pindex);
		return NodalEnum;
	}
	else{
		/*We don't have _ in the name, this is a regular variable: */
		memcpy(root,pch,(strlen(pch)+1)*sizeof(char));
		*pindex=-1;
		return RegularEnum;
	}
}
