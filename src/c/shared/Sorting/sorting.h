/*!\file: sorting.h
 * \brief prototypes for sorting.h
 */ 

#ifndef _SORTING_H_
#define  _SORTING_H_

int binary_search(int* poffset,int target,int* sorted_integers,int num_integers);
template <typename doubletype> int binary_search(int* poffset,doubletype target,doubletype* list,int length){ /*{{{*/
	/*
	 *             l[0]  l[1]  l[2]        l[n]  l[n+1]   l[length-1]
	 *     <-------+-----+-----+-----+ ... +-----+........+-------------->
	 * offset: -1     0     1     2           n              length -1
	 *  
	 *  offset = -1        target < list[0]
	 *  offset = n         list[n] <= target < list[n+1]
	 *  offset = length-1  list[length-1] <= target
	 */

	/*output: */
	int offset = 0;
	int found  = 0;

	/*intermediary: */
	int n0 = 0;
	int n1 = int(length/2);
	int n2 = length-1;

	if(target<list[n0]){
		found  = 1;
		offset = -1;
	}
	else if(target>=list[n2]){
		found  = 1;
		offset = length-1;
	}
	else{
		for(;;){
			/*did we find the target?*/
			if(list[n1]<=target && list[n1+1]>target){
				found  = 1;
				offset = n1;
				break;
			}
			else if(target < list[n1]){
				n2 = n1;
				n1 = n0 + int((n2-n0)/2);
			}
			else{
				n0 = n1;
				n1 = n0 + int((n2-n0)/2);
			}
		}
	}

	/*Assign output pointer and return*/
	*poffset=offset;
	return found;
} /*}}}*/

#endif //ifndef _SORTING_H_
