/*!\file GetVectorFromControlInputsx
 * \brief retrieve vector from inputs in elements
 */

#include "./GetVectorFromControlInputsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void GetVectorFromControlInputsx(Vector<IssmDouble>** pvector, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters,const char* data){/*{{{*/

	int  num_controls;
	int* N = NULL;
	int* M = NULL;
	int* control_type = NULL;

	/*Retrieve some parameters*/
	parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);
	parameters->FindParam(&N,NULL,ControlInputSizeNEnum);
	parameters->FindParam(&M,NULL,ControlInputSizeMEnum);

	/*1. Get vector size*/
	int size = 0;
	for(int i=0;i<num_controls;i++) size+=M[i]*N[i];

	/*2. Allocate vector*/
	Vector<IssmDouble>* vector=new Vector<IssmDouble>(size);

	/*3. Populate vector*/
	int offset = 0;
	for(int i=0;i<num_controls;i++){

		/*Is the control a Param?*/
		if(IsParamEnum(control_type[i])){
			parameters->GetVectorFromControl(vector,control_type[i],i,N[i],data,offset);
		}
		else if(IsInputEnum(control_type[i])){
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->GetVectorFromControlInputs(vector,control_type[i],i,N[i],data,offset);
			}
		}
		else{
			_error_("not supported yet");
		}
		offset += M[i]*N[i];
	}
	vector->Assemble();

	/*Assign output pointers:*/
	xDelete<int>(control_type);
	xDelete<int>(M);
	xDelete<int>(N);
	*pvector=vector;
}/*}}}*/
void GetVectorFromControlInputsx( IssmDouble** pvector,int *pN, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters, const char* data){/*{{{*/

	/*intermediary: */
	int N;
	Vector<IssmDouble>* vec_vector=NULL;

	/*Get PETSc vector*/
	GetVectorFromControlInputsx( &vec_vector, elements,nodes, vertices, loads, materials, parameters,data);

	/*Serialize*/
	vec_vector->GetSize(&N);
	IssmDouble* vector=vec_vector->ToMPISerial();
	delete vec_vector;

	/*Assign output pointers:*/
	*pvector=vector;
	if(pN) *pN=N;
}/*}}}*/

/*For autodiff, we sometimes need to cast our vectors to passive*/
#ifdef _HAVE_AD_
void GetPassiveVectorFromControlInputsx(IssmPDouble** pvector,int* pN, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters, const char* data){/*{{{*/

	/*Get active vector first*/
	Vector<IssmDouble> *activevector = NULL;
	IssmPDouble        *vector       = NULL;
	int                size;

	/*Retrieve some parameters*/
	GetVectorFromControlInputsx(&activevector, elements,nodes, vertices, loads, materials, parameters,data);

	/*Serialize vector*/
	activevector->GetSize(&size);
	IssmDouble* dactivevector=activevector->ToMPISerial();

	/*Cast to passive*/
	vector=xNew<IssmPDouble>(size);
	for(int i=0;i<size;i++) vector[i] = reCast<IssmPDouble>(dactivevector[i]);

	/*Assign output pointers:*/
	delete activevector;
	xDelete<IssmDouble>(dactivevector);
	*pvector=vector;
	if(pN) *pN=size;

}/*}}}*/
#else
void GetPassiveVectorFromControlInputsx(IssmPDouble** pvector,int* pN, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters, const char* data){/*{{{*/

	GetVectorFromControlInputsx(pvector,pN,elements,nodes, vertices, loads, materials, parameters,data);
}/*}}}*/
#endif
