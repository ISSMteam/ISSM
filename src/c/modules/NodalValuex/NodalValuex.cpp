/*!\file NodalValuex
 * \brief: compute value at certain node
 */

#include "./NodalValuex.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void NodalValuex( IssmDouble* pnodalvalue, int natureofdataenum,Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters){

	IssmDouble value;
	int        index;
	int        found,sumfound,cpu_found,cpu;

	/*retrieve element we are interested in: */
	parameters->FindParam(&index,IndexEnum);

	/*This is the vertex id for which we want to collect the data. Go through elements, and for each 
	 *element, figure out  if they hold the vertex, and the data. If so, return it: */
	cpu_found=-1;
	found=0;

	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		found=element->NodalValue(&value,index,natureofdataenum);
		if(found){
			cpu_found=IssmComm::GetRank();
			break;
		}
	}

	/*Broadcast whether we found the element: */
	ISSM_MPI_Allreduce(&found,&sumfound,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
	if(!sumfound)_error_("could not find element with vertex with id " << index << " to compute nodal value " << EnumToStringx(natureofdataenum));

	/*Broadcast and plug into response: */
	ISSM_MPI_Allreduce ( &cpu_found,&cpu,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());
	ISSM_MPI_Bcast(&value,1,ISSM_MPI_DOUBLE,cpu,IssmComm::GetComm()); 

	*pnodalvalue=value;
}
