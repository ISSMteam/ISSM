/*!\file Masscon.h
 * \brief: header file for Masscon object
 */

#ifndef _MASSCON_H_
#define _MASSCON_H_

/*Headers:*/
/*{{{*/
#include "./Definition.h"
#include "../datastructures/datastructures.h"
#include "./Elements/Element.h"
#include "./Elements/Elements.h"
#include "./FemModel.h"
#include "../classes/Params/Parameters.h"
/*}}}*/

class Masscon: public Object, public Definition{

	public: 

		int         definitionenum;
		char*       name;
		IssmDouble* levelset;
		int         M;

		/*Masscon constructors, destructors :*/
		Masscon(){/*{{{*/

			this->definitionenum = -1;
			this->name = NULL;
			this->levelset=NULL;
			this->M=0;

		}
		/*}}}*/
		Masscon(char* in_name, int in_definitionenum, IssmDouble* levelsetin, int Min){ /*{{{*/

			this->definitionenum=in_definitionenum;
			this->name   = xNew<char>(strlen(in_name)+1);
			xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

			this->levelset   = xNew<IssmDouble>(Min);
			xMemCpy<IssmDouble>(this->levelset, levelsetin, Min);

			this->M=Min;

		}
		/*}}}*/
		~Masscon(){/*{{{*/
			if(this->name)xDelete(this->name); 
			if(this->levelset)xDelete(this->levelset);
		}
		/*}}}*/
		/*Object virtual function resolutoin: */
		Object* copy() {/*{{{*/
			Masscon* mf = new Masscon(this->name,this->definitionenum,this->levelset,this->M);
			return (Object*) mf;
		}
		/*}}}*/
		void DeepEcho(void){/*{{{*/
			this->Echo();
		}
		/*}}}*/
		void Echo(void){/*{{{*/
			_printf_(" Masscon: " << this->name << " " << this->definitionenum << "\n");
			_printf_("    levelset: " << this->levelset << "\n");
			_printf_("    M: " << this->M << "\n");
		}
		/*}}}*/
		int Id(void){/*{{{*/
			return -1;
		}
		/*}}}*/
		void Marshall(MarshallHandle* marshallhandle){/*{{{*/
			_error_("not implemented yet!"); 
		} 
		/*}}}*/
		int ObjectEnum(void){/*{{{*/
			return MassconEnum;
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

			 int i;
			 IssmDouble mass_t=0.;
			 IssmDouble all_mass_t=0.;

			 for(i=0;i<femmodel->elements->Size();i++){
				 Element* element=(Element*)femmodel->elements->GetObjectByOffset(i);
				 mass_t+=element->Masscon(this->levelset);
			 }

			 ISSM_MPI_Allreduce ( (void*)&mass_t,(void*)&all_mass_t,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
			 mass_t=all_mass_t;

			 return mass_t;
		 }
			/*}}}*/
};

#endif  /* _MASSCON_H_ */
