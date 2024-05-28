#ifndef _HEAPSORT_H_
#define _HEAPSORT_H_

/*Sort a list of size n*/
template<class T> inline void  HeapSort(T *c,long n){ /*{{{*/

	/*Intermediaries*/
	int i,j,l,r;
	T   crit;

	/*return if size <=1*/
	if(n<=1) return;

	/*Initialize variables*/
	l=n/2+1; 
	r=n;
	c--; //the array must starts at 1 and not 0 

	/*Sorting algorithm*/
	for(;;){
		if(l<=1){
			crit   = c[r];
			c[r--] = c[1];
			if(r==1){
				c[1]=crit; 
				return;
			}
		}
		else{
			crit = c[--l]; 
		}
		j=l;
		for(;;){
			i = j;
			j = 2*j;
			if  (j>r){c[i]=crit;break;}
			if ((j<r) && (c[j] < c[j+1])) j++;//c[j+1]> c[j] -> take j+1 instead of j (larger value)
			if (crit < c[j]) c[i]=c[j];       //c[j]  > crit -> put this large value in i(<j)
			else{c[i]=crit;break;}            //c[j]  < crit -> put crit in i (<j)
		}
	}
}
/*}}}*/

/*Sort a list of size n and returns ordering*/
template<class T> inline void  HeapSort(int** porder,T* c,int n){ /*{{{*/

	/*Intermediaries*/
	int  i,j,l,r;
	T    crit;
	int  pos;
	int* order = NULL;

	/*return if size <=1*/
	if(n<=1) return;

	/*Initialize variables*/
	l=n/2+1; 
	r=n;
	c--; //the array must starts at 1 and not 0 
	order = new int[n];
	for(i=0;i<n;i++) order[i]=i;
	order--;

	/*Sorting algorithm*/
	for(;;){
		if(l<=1){
			crit  =c[r]; pos=order[r];
			c[r--]=c[1]; order[r+1]=order[1];
			if (r==1){
				c[1]=crit; order[1]=pos;
				order++;
				*porder=order;
				return;
			}
		}
		else  {crit=c[--l]; pos=order[l];}
		j=l;
		for(;;){
			i=j;
			j=2*j;
			if  (j>r) {c[i]=crit;order[i]=pos;break;}
			if ((j<r) && (c[j] < c[j+1]))j++;
			if (crit < c[j]){
				c[i]=c[j];
				order[i]=order[j];
			}
			else{
				c[i]=crit;order[i]=pos;
				break;
			}
		}
	}

	/*Make cppcheck happy*/
	*porder=order;
}/*}}}*/

#endif
