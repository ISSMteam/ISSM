/*!\file Contour.h
 * \brief: header file for Contour object
 */

#ifndef _CONTOUR_H_
#define _CONTOUR_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "../datastructures/datastructures.h"
/*}}}*/

template <class doubletype>
class Contour: public Object{

	public: 

		int         id;
		int         nods;     //number of vertices in the contour
		doubletype *x;
		doubletype *y;
		bool        closed;   //is this contour closed?

		/*Contour constructors, destructors :*/
		Contour(){/*{{{*/
			this->id     = 0;
			this->nods   = 0;
			this->x      = NULL;
			this->y      = NULL;
			this->closed = false;
		}
		/*}}}*/
		Contour(int pid,int pnods, doubletype* px, doubletype* py,bool pclosed){/*{{{*/

			this->id     = pid;
			this->nods   = pnods;
			this->closed = pclosed;
			if(nods){
				this->x=xNew<doubletype>(nods);
				xMemCpy<doubletype>(this->x,px,nods);
				this->y=xNew<doubletype>(nods);
				xMemCpy<doubletype>(this->y,py,nods);
			}
		}
		/*}}}*/
		~Contour(){/*{{{*/
			xDelete<doubletype>(this->x);
			xDelete<doubletype>(this->y);
		}
		/*}}}*/

		/*Object virtual function resolutoin: */
		Object* copy() {/*{{{*/

			Contour* contour = new Contour(this->id,this->nods,this->x,this->y,this->closed);

			return (Object*) contour;
		}
		/*}}}*/
		void DeepEcho(void){/*{{{*/
			this->Echo();
		}
		/*}}}*/
		void Echo(void){/*{{{*/
			_printf_(" Contour: " << id << "\n");
			_printf_("    nods: " << nods << "\n");
			_printf_("  closed: " << (closed?"true":"false") << "\n");
			if(nods){
				_printf_("   x , y:\n");
				for(int i=0;i<nods;i++){
					_printf_(i << ": " << x[i] << " | " << y[i] << "\n");
				}
			}
		}
		/*}}}*/
		int Id(void){/*{{{*/
			return id;
		}
		/*}}}*/
		void Marshall(MarshallHandle* marshallhandle){/*{{{*/
			_error_("not implemented yet!"); 
		} 
		/*}}}*/
		int ObjectEnum(void){/*{{{*/
			return ContourEnum;
		}
		/*}}}*/

};

#endif  /* _CONTOUR_H_ */
