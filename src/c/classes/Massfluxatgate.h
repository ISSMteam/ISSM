/*!\file Massfluxatgate.h
 * \brief: header file for Massfluxatgate object
 */

#ifndef _MASSFLUXATGATE_H_
#define _MASSFLUXATGATE_H_

/*Headers:*/
/*{{{*/
#include "./Definition.h"
#include "../datastructures/datastructures.h"
#include "./Elements/Element.h"
#include "./Elements/Elements.h"
#include "./FemModel.h"
/*}}}*/

template <class doubletype> 
class Massfluxatgate: public Object, public Definition{

	public: 

		int         definitionenum;
		char*       name;
		int         numsegments;
		doubletype *x1;
		doubletype *y1;
		doubletype *x2;
		doubletype *y2;
		int*        elements;

		/*Massfluxatgate constructors, destructors :*/
		Massfluxatgate(){/*{{{*/
			this->definitionenum        = -1;
			this->name        = NULL;
			this->numsegments = 0;
			this->elements    = NULL;
			this->x1				= NULL;
			this->x2				= NULL;
			this->y1				= NULL;
			this->y2				= NULL;
		}
		/*}}}*/
		Massfluxatgate(char* in_name, int in_definitionenum, int in_numsegments, doubletype* in_segments) {/*{{{*/

			int i;

			this->definitionenum=in_definitionenum;

			this->name   = xNew<char>(strlen(in_name)+1);
			xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

			this->numsegments=in_numsegments;

			if(this->numsegments){
				this->x1=xNew<doubletype>(this->numsegments);
				this->x2=xNew<doubletype>(this->numsegments);
				this->y1=xNew<doubletype>(this->numsegments);
				this->y2=xNew<doubletype>(this->numsegments);
				this->elements=xNew<int>(this->numsegments);

				for(i=0;i<this->numsegments;i++){
					this->x1[i]=in_segments[5*i+0];
					this->y1[i]=in_segments[5*i+1];
					this->x2[i]=in_segments[5*i+2];
					this->y2[i]=in_segments[5*i+3];
					this->elements[i]=reCast<int,doubletype>(in_segments[5*i+4]);
				}
			}
		}
		/*}}}*/
		Massfluxatgate(char* in_name, int in_definitionenum, int in_numsegments, doubletype* in_x1, doubletype* in_y1, doubletype* in_x2, doubletype* in_y2,int* in_elements){/*{{{*/

			this->definitionenum=in_definitionenum;
			this->name   = xNew<char>(strlen(in_name)+1);
			xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

			this->numsegments=in_numsegments;

			if(this->numsegments){
				this->x1=xNew<doubletype>(this->numsegments); xMemCpy<doubletype>(this->x1,in_x1,this->numsegments);
				this->y1=xNew<doubletype>(this->numsegments); xMemCpy<doubletype>(this->y1,in_y1,this->numsegments);
				this->x2=xNew<doubletype>(this->numsegments); xMemCpy<doubletype>(this->x2,in_x2,this->numsegments);
				this->y2=xNew<doubletype>(this->numsegments); xMemCpy<doubletype>(this->y2,in_y2,this->numsegments);
				this->elements=xNew<int>(this->numsegments); xMemCpy<int>(this->elements,in_elements,this->numsegments);

			}
		}
		/*}}}*/
		~Massfluxatgate(){/*{{{*/
			if(this->numsegments){
				xDelete<doubletype>(this->x1);
				xDelete<doubletype>(this->y1);
				xDelete<doubletype>(this->x2);
				xDelete<doubletype>(this->y2);
				xDelete<int>(this->elements);
			}
			xDelete<char>(this->name);
		}
		/*}}}*/

		/*Object virtual function resolutoin: */
		Object* copy() {/*{{{*/
			return new Massfluxatgate(this->name,this->definitionenum,this->numsegments,this->x1,this->y1,this->x2,this->y2,this->elements); 
		}
		/*}}}*/
		void DeepEcho(void){/*{{{*/
			this->Echo();
		}
		/*}}}*/
		void Echo(void){/*{{{*/
			_printf_(" Massfluxatgate: " << name << " " << this->definitionenum << "\n");
			_printf_("    numsegments: " << numsegments << "\n");
			if(numsegments){
				_printf_("   element: x1, y1, x2, y2:\n");
				for(int i=0;i<numsegments;i++){
					_printf_(elements[i] << " " << x1[i] << " " << y1[i] << " " << x2[i] << " " << y2[i] << "\n");
				}
			}
		}
		/*}}}*/
		int Id(void){/*{{{*/
			return -1;
		}
		/*}}}*/
		int ObjectEnum(void){/*{{{*/
			return MassfluxatgateEnum;
		}
		/*}}}*/
		void Marshall(MarshallHandle* marshallhandle){/*{{{*/

			int object_enum = MassfluxatgateEnum;
			marshallhandle->call(object_enum);

			marshallhandle->call(this->definitionenum);
			marshallhandle->call(this->name);
			marshallhandle->call(this->numsegments);
			marshallhandle->call(this->x1,this->numsegments);
			marshallhandle->call(this->x2,this->numsegments);
			marshallhandle->call(this->y1,this->numsegments);
			marshallhandle->call(this->y2,this->numsegments);
			marshallhandle->call(this->elements,this->numsegments);
		} 
		/*}}}*/
		/*Definition virtual function resolutoin: */
		int DefinitionEnum(){/*{{{*/

			return this->definitionenum;
		}
		/*}}}*/
		char* Name(){/*{{{*/

			char* name2=xNew<char>(strlen(this->name)+1);
			xMemCpy(name2,this->name,strlen(this->name)+1);

			return name2;
		}
		/*}}}*/
		 IssmDouble Response(FemModel* femmodel){/*{{{*/

			int          i,j;
			Element     *element       = NULL;
			IssmDouble mass_flux     = 0;
			IssmDouble all_mass_flux = 0;

			/*Go through segments, and then elements, and figure out which elements belong to a segment. 
			 * When we find one, use the element to compute the mass flux on the segment: */
			for(i=0;i<numsegments;i++){
				for(j=0;j<femmodel->elements->Size();j++){
					element=(Element*)femmodel->elements->GetObjectByOffset(j);
					if (element->Id()==this->elements[i]){
						/*We found the element which owns this segment, use it to compute the mass flux: */
						mass_flux+=element->MassFlux(x1[i],y1[i],x2[i],y2[i],elements[i]);
						break;
					}
				}
			}

			ISSM_MPI_Allreduce ( (void*)&mass_flux,(void*)&all_mass_flux,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
			mass_flux=all_mass_flux;
			return mass_flux;
		}
			/*}}}*/
};

#endif  /* _MASSFLUXATGATE_H_ */
