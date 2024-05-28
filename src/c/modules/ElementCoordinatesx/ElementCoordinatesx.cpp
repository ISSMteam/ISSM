/*!\file ElementCoordinatesx
 * \brief: compute a vector xe,ye and ze of element centroid coordinates by
 * marching through all our elements.
 */

#include "./ElementCoordinatesx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void ElementCoordinatesx( IssmDouble** pxe, IssmDouble** pye, IssmDouble** pze, IssmDouble** pareae, Elements* elements,bool spherical) { /*{{{*/

	/*figure out how many vertices we have: */
	int numberofelements=elements->NumberOfElements();

	Vector<IssmDouble>* vxe=new Vector<IssmDouble>(numberofelements);
	Vector<IssmDouble>* vye=new Vector<IssmDouble>(numberofelements);
	Vector<IssmDouble>* vze=new Vector<IssmDouble>(numberofelements);
	Vector<IssmDouble>* vareae=new Vector<IssmDouble>(numberofelements);

	/*march through our elements: */
	for(Object* & object : elements->objects){
		Element* element=(Element*)object;
		element->ElementCoordinates(vxe,vye,vze,vareae,spherical);
	}

	/*Assemble*/
	vxe->Assemble();
	vye->Assemble();
	vze->Assemble();
	vareae->Assemble();

	/*serialize: */
	IssmDouble* xe=vxe->ToMPISerial();
	IssmDouble* ye=vye->ToMPISerial();
	IssmDouble* ze=vze->ToMPISerial();
	IssmDouble* areae=vareae->ToMPISerial();

	/*Free resources: */
	delete vxe;
	delete vye;
	delete vze;
	delete vareae;

	/*output: */
	if(pxe) *pxe=xe;
	else xDelete<IssmDouble>(xe);
	if(pye) *pye=ye;
	else xDelete<IssmDouble>(ye);
	if(pze) *pze=ze;
	else xDelete<IssmDouble>(ze);
	if(pareae) *pareae=areae;
	else xDelete<IssmDouble>(areae);

} /*}}}*/
void ElementCoordinatesx( IssmDouble** plonge, IssmDouble** plate, IssmDouble** pareae, Elements* elements) { /*{{{*/

	/*figure out how many vertices we have: */
	int numberofelements=elements->NumberOfElements();

	Vector<IssmDouble>* vlonge=new Vector<IssmDouble>(numberofelements);
	Vector<IssmDouble>* vlate=new Vector<IssmDouble>(numberofelements);
	Vector<IssmDouble>* vareae=new Vector<IssmDouble>(numberofelements);

	/*march through our elements: */
	for(Object* & object : elements->objects){
		Element* element=(Element*)object;
		element->ElementCoordinates(vlonge,vlate,vareae);
	}

	/*Assemble*/
	vlonge->Assemble();
	vlate->Assemble();
	vareae->Assemble();

	/*serialize: */
	IssmDouble* longe=vlonge->ToMPISerial();
	IssmDouble* late=vlate->ToMPISerial();
	IssmDouble* areae=vareae->ToMPISerial();

	/*Free resources: */
	delete vlonge;
	delete vlate;
	delete vareae;

	/*output: */
	if(plonge) *plonge=longe;
	else xDelete<IssmDouble>(longe);
	if(plate) *plate=late;
	else xDelete<IssmDouble>(late);
	if(pareae) *pareae=areae;
	else xDelete<IssmDouble>(areae);

} /*}}}*/
