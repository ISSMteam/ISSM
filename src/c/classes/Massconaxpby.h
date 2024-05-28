/*!\file Massconaxpby.h
 * \brief: header file for Massconaxpby object
 */

#ifndef _MASSCON_AXPBY_H_
#define _MASSCON_AXPBY_H_

/*Headers:*/
/*{{{*/
#include "./Definition.h"
#include "../datastructures/datastructures.h"
#include "./Elements/Element.h"
#include "./Elements/Elements.h"
#include "./FemModel.h"
#include "../classes/Params/Parameters.h"
int OutputDefinitionsResponsex(IssmDouble* pres, FemModel* femmodel,const char* output_string);
/*}}}*/
class Massconaxpby: public Object, public Definition{

	public: 

		int         definitionenum;
		char*       name;
		char*       namex;
		char*       namey;
		IssmDouble  alpha;
		IssmDouble  beta;

		/*Massconaxpby constructors, destructors :*/
		Massconaxpby(){/*{{{*/

			this->definitionenum = -1;
			this->name = NULL;
			this->namex = NULL;
			this->namey = NULL;
			this->alpha=UNDEF;
			this->beta=UNDEF;

		}
		/*}}}*/
		Massconaxpby(char* in_name,int in_definitionenum, char* in_namex, char* in_namey, IssmDouble in_alpha,IssmDouble in_beta){ /*{{{*/

			this->definitionenum = in_definitionenum;
			this->name   = xNew<char>(strlen(in_name)+1);
			xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

			this->namex   = xNew<char>(strlen(in_namex)+1);
			xMemCpy<char>(this->namex,in_namex,strlen(in_namex)+1);

			this->namey   = xNew<char>(strlen(in_namey)+1);
			xMemCpy<char>(this->namey,in_namey,strlen(in_namey)+1);

			this->alpha=in_alpha;
			this->beta=in_beta;

		}
		/*}}}*/
		~Massconaxpby(){/*{{{*/
			if(this->name)xDelete(this->name); 
			if(this->namex)xDelete(this->namex); 
			if(this->namey)xDelete(this->namey); 
		}
		/*}}}*/
		/*Object virtual function resolutoin: */
		Object* copy() {/*{{{*/
			Massconaxpby* mf = new Massconaxpby(this->name,this->definitionenum,this->namex,this->namey, this->alpha, this->beta);
			return (Object*) mf;
		}
		/*}}}*/
		void DeepEcho(void){/*{{{*/
			this->Echo();
		}
		/*}}}*/
		void Echo(void){/*{{{*/
			_printf_(" Massconaxpby: " << this->name << " " << this->definitionenum << "\n");
			_printf_("    namex: " << this->namex << "\n");
			_printf_("    namey: " << this->namey << "\n");
			_printf_("    alpha: " << this->alpha << "\n");
			_printf_("    beta: " << this->beta << "\n");
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
			return MassconaxpbyEnum;
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

			 int ierr;
			 IssmDouble xresponse,yresponse;

			 /*Get response from both masscons: */
			 ierr=OutputDefinitionsResponsex(&xresponse,femmodel,this->namex);
			 if(ierr) _error_("not found");
			 ierr=OutputDefinitionsResponsex(&yresponse,femmodel,this->namey);
			 if(ierr) _error_("not found");

			 return this->alpha*xresponse+this->beta*yresponse;
		 }
			/*}}}*/
};

#endif  /* _MASSCON_H_ */
