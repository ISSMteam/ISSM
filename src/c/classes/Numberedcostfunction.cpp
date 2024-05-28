/*!\file Numberedcostfunction.cpp
 * \brief: implementation for the Numberedcostfunction object
 */
/*Include files: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Headers:*/
//#include "./Definition.h"
//#include "../datastructures/datastructures.h"
#include "./classes.h"
#include "./Elements/Element.h"
#include "./Elements/Elements.h"
#include "./FemModel.h"
#include "./ExternalResults/ExternalResult.h"
#include "./ExternalResults/Results.h"
#include "../modules/SurfaceAbsVelMisfitx/SurfaceAbsVelMisfitx.h"
#include "../modules/SurfaceRelVelMisfitx/SurfaceRelVelMisfitx.h"
#include "../modules/SurfaceLogVelMisfitx/SurfaceLogVelMisfitx.h"
#include "../modules/SurfaceLogVxVyMisfitx/SurfaceLogVxVyMisfitx.h"
#include "../modules/ThicknessAbsMisfitx/ThicknessAbsMisfitx.h"
#include "../modules/ThicknessAlongGradientx/ThicknessAlongGradientx.h"
#include "../modules/ThicknessAcrossGradientx/ThicknessAcrossGradientx.h"
#include "../modules/RheologyBbarAbsGradientx/RheologyBbarAbsGradientx.h"
#include "../modules/DragCoefficientAbsGradientx/DragCoefficientAbsGradientx.h"

/*}}}*/

/*Numberedcostfunction constructors, destructors :*/
Numberedcostfunction::Numberedcostfunction(){/*{{{*/

	this->definitionenum = -1;
	this->name = NULL;
	this->number_cost_functions = -1;
	this->cost_functions_list = NULL;

}
/*}}}*/
Numberedcostfunction::Numberedcostfunction(char* in_name, int in_definitionenum,int number_cost_functions_in,int* cost_functions_list_in){/*{{{*/
	_assert_(number_cost_functions_in>0); 
	_assert_(cost_functions_list_in); 

	this->definitionenum=in_definitionenum;
	this->name   = xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->number_cost_functions = number_cost_functions_in;
	this->cost_functions_list = xNew<int>(number_cost_functions_in);

	for(int i=0;i<number_cost_functions_in;i++){
		this->cost_functions_list[i] = cost_functions_list_in[i];
	}
}
/*}}}*/
Numberedcostfunction::~Numberedcostfunction(){/*{{{*/
	xDelete<int>(this->cost_functions_list);
	if(this->name)xDelete(this->name);
}
/*}}}*/

/*Object virtual function resolutoin: */
Object* Numberedcostfunction::copy() {/*{{{*/
	Numberedcostfunction* out = new Numberedcostfunction(this->name,this->definitionenum,this->number_cost_functions,this->cost_functions_list);
	return (Object*)out;
}
/*}}}*/
void Numberedcostfunction::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Numberedcostfunction::Echo(void){/*{{{*/
	_printf_(" Numberedcostfunction: " << this->name << " " << this->definitionenum << "\n");
	_printf_("    number_cost_functions: "<<this->number_cost_functions<<"\n");
	_printf_("    ");
	for(int i=0;i<this->number_cost_functions;i++){
		_printf_(this->cost_functions_list[i]<< "  ");
	}
	_printf_("\n");
}
/*}}}*/
int Numberedcostfunction::Id(void){/*{{{*/
	return -1;
}
/*}}}*/
void Numberedcostfunction::Marshall(MarshallHandle* marshallhandle){/*{{{*/
	_error_("not implemented yet!"); 
} 
/*}}}*/
int Numberedcostfunction::ObjectEnum(void){/*{{{*/
	return NumberedcostfunctionEnum;
}
/*}}}*/

/*Definition virtual function resolutoin: */
int Numberedcostfunction::DefinitionEnum(){/*{{{*/
	return this->definitionenum;
}
/*}}}*/
char* Numberedcostfunction::Name(){/*{{{*/

	char* name2=xNew<char>(strlen(this->name)+1);
	xMemCpy(name2,this->name,strlen(this->name)+1);

	return name2;
}
/*}}}*/
IssmDouble Numberedcostfunction::Response(FemModel* femmodel){/*{{{*/

	 _assert_(number_cost_functions>0 && number_cost_functions<1e3); 

	 /*output:*/
	 IssmDouble value;
	 IssmDouble value_sum = 0.;

		/*Scalar control output*/
	 for(int i=0;i<this->number_cost_functions;i++){
		 switch(this->cost_functions_list[i]){
			 case SurfaceAbsVelMisfitEnum:
				 SurfaceAbsVelMisfitx(&value,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters);        
				 break;
			 case SurfaceRelVelMisfitEnum:            
				 SurfaceRelVelMisfitx(&value, femmodel->elements,femmodel->nodes,femmodel-> vertices,femmodel-> loads,femmodel-> materials,femmodel->parameters); 
				 break;
			 case SurfaceLogVelMisfitEnum:            
				 SurfaceLogVelMisfitx(&value,femmodel-> elements,femmodel->nodes,femmodel-> vertices,femmodel-> loads,femmodel-> materials,femmodel->parameters); 
				 break;
			 case SurfaceLogVxVyMisfitEnum:           
				 SurfaceLogVxVyMisfitx(&value,femmodel-> elements,femmodel->nodes,femmodel-> vertices,femmodel-> loads,femmodel-> materials,femmodel->parameters); 
				 break;
			 case ThicknessAbsMisfitEnum:             
				 ThicknessAbsMisfitx(&value,femmodel-> elements,femmodel->nodes,femmodel-> vertices,femmodel-> loads,femmodel-> materials,femmodel-> parameters); 
				 break;
			 case ThicknessAbsGradientEnum:             
				 femmodel->ThicknessAbsGradientx(&value);
				 break;
			 case ThicknessAlongGradientEnum:         
				 ThicknessAlongGradientx(&value,femmodel-> elements,femmodel->nodes,femmodel-> vertices,femmodel-> loads,femmodel-> materials,femmodel-> parameters); 
				 break;
			 case ThicknessAcrossGradientEnum:        
				 ThicknessAcrossGradientx(&value,femmodel-> elements,femmodel->nodes,femmodel-> vertices,femmodel-> loads,femmodel-> materials,femmodel-> parameters); 
				 break;
			 case RheologyBbarAbsGradientEnum:        
				 RheologyBbarAbsGradientx(&value,femmodel-> elements,femmodel->nodes,femmodel-> vertices,femmodel-> loads,femmodel-> materials,femmodel-> parameters); 
				 break;
			 case DragCoefficientAbsGradientEnum:     
				 DragCoefficientAbsGradientx(&value,femmodel-> elements,femmodel->nodes,femmodel-> vertices,femmodel-> loads,femmodel-> materials,femmodel-> parameters); 
				 break;
			 default:
				 _error_(EnumToStringx(this->cost_functions_list[i])<<" not supported");
		 }
		 _printf0_("#"<<i+1<<": "<<value<<" ");
		 value_sum += value;
 }
	 _printf0_("\n");

	 /*done:*/
	return value_sum;
 }
 /*}}}*/
