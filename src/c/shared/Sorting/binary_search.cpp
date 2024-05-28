/*!\file:  binary_search.cpp
 * \brief  binary search on an integer array, that is already sorted
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../Exceptions/exceptions.h"
#include <stdio.h>

int binary_search(int* poffset,int target,int* sorted_integers,int num_integers){ /*{{{*/

	/*output: */
	int offset=0;  //offset, if found
	int found=0;   //found=0 if target is not found, 1 otherwise.

	/*intermediary: */
	int *beg = NULL;
	int *end = NULL;
	int *mid = NULL;

	_assert_(sorted_integers);

	// point to beginning and end of the array
	beg=sorted_integers;
	end=sorted_integers+num_integers;
	mid=beg+(int)(num_integers/2);

	if (*beg==target){
		found=1;
		offset=0;
	}
	else if(*(end-1)==target){
		found=1;
		offset=num_integers-1;
	}
	else{
		while((beg <= end) && (*mid != target)){
			// is the target in lower or upper half?
			if (target < *mid) {
				end = mid - 1;  //new end
				mid = beg + (end-beg)/2;  //new middle
			}
			else {
				beg = mid + 1;  //new beginning
				mid = beg + (end-beg)/2;  //new middle
			}
		}

		//did we find the target?
		if (*mid == target) {
			found=1;
			offset=mid-sorted_integers;
		}
		else {
			found=0;
		}
	}

	/*Assign output pointers:*/
	*poffset=offset;

	/*Return result: */
	return found;
} /*}}}*/
