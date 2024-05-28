/*!\file:  SparseRow: 
 * \brief implementation of a sparse row, which can then be used to make a sparse matrix.
 */ 

#ifndef _SPARSE_ROW_H_
#define _SPARSE_ROW_H_

/*Headers:*/
#include "../toolkitsenums.h"
#include "../../shared/shared.h"

template <class doubletype> 
class SparseRow{

	public:

		int         M; //real size
		int         ncols; //number of non-zeros 
		int*        indices;
		doubletype* values;

		/*SparseRow constructors, destructors*/
		SparseRow(){ /*{{{*/
			 M=0;
			 ncols=0;
			 indices=NULL;
			 values=NULL;
		} /*}}}*/
		SparseRow(int in_M){/*{{{*/

			M=in_M;
			ncols=0;
			indices=NULL;
			values=NULL;

		} /*}}}*/
		~SparseRow(){/*{{{*/
			if(ncols){
				xDelete<int>(indices);
				xDelete<doubletype>(values);
			}
		} /*}}}*/

		/*SparseRow specific routines*/
		void Echo(){ /*{{{*/
			int i;

			for(i=0;i<ncols;i++){
				_printf_("(" << indices[i] << "," << values[i] << ") ");
			}
		} /*}}}*/
		void SetValues(int numvalues,int* cols,doubletype* vals, int* mods){ /*{{{*/

			int count;
			int i,j;

			if(!M)_error_("unknow dimension for this sparse row!");

			/*Deallocate if already allocated: */
			if(ncols){
				xDelete<int>(indices);
				xDelete<doubletype>(values);
			}

			/*check numvalues: */
			if(!numvalues)return;

			/*Go through cols and resolve duplicates: */
			for(i=0;i<numvalues;i++){
				for(j=i+1;j<numvalues;j++){
					if (cols[j]==cols[i]){
						if (mods[j]==ADD_VAL){
							vals[i]+=vals[j];
						}
						else vals[i]=vals[j];
						cols[j]=-1;
					}
					/*Ensure that this value will not be used anymore: */
				}
			}

			/*Now go through cols once more, and retrieve only what we need: */
			ncols=0;
			for(i=0;i<numvalues;i++)if(cols[i]>=0)ncols++;

			/*Allocate and fill: */
			indices=xNew<int>(ncols); _assert_(indices);
			values=xNewZeroInit<doubletype>(ncols);  _assert_(values);

			count=0;
			for(i=0;i<numvalues;i++){
				if(cols[i]>=0){
					indices[count]=cols[i];
					values[count]=vals[i];
					count++;
				}
			}

			if(count!=ncols)_error_("counter problem during set values operations");
		} /*}}}*/
		doubletype Norm(NormMode mode){ /*{{{*/

			int i;
			doubletype norm=0.;

			switch(mode){
				case NORM_INF:
					for(i=0;i<ncols;i++){
						norm+=fabs(values[i]);
					}
					return norm;
					break; 
				case NORM_FROB:
					for(i=0;i<ncols;i++){
						norm+=values[i]*values[i];
					}
					return norm;
					break; 

				default:
					_error_("unknown norm !");
					break;
			}
		}
		/*}}}*/
		doubletype Mult(doubletype* X){ /*{{{*/

			int i;
			doubletype mult=0;

			for(i=0;i<ncols;i++){
				mult+=values[i]*X[indices[i]];
			}

			return mult;
		}
		/*}}}*/
		int Nnz(void){ /*{{{*/

			return ncols;
		}
		/*}}}*/
		void SetIrnJcnA(int* irn_loc,int* jcn_loc,doubletype* a_loc,int i_index,int count){/*{{{*/
			int i;

			for(i=0;i<ncols;i++){
				irn_loc[count+i]=i_index;
				jcn_loc[count+i]=indices[i]+1; //fortran indexing
				a_loc[count+i]=values[i];
			}
		}
		/*}}}*/
};
#endif //#ifndef _SPARSE_ROW_H_
