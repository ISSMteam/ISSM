/*
 * UpdateElementsAndMaterialsControl:
 */

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"
#include "../../MeshPartitionx/MeshPartitionx.h"
#include "../ModelProcessorx.h"

void	UpdateElementsAndMaterialsDakota(Elements* elements,Inputs* inputs,Materials* materials, IoModel* iomodel){

	/*recover parameters: */
	bool dakota_analysis;
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");

	if(dakota_analysis) iomodel->FetchDataToInput(inputs,elements,"md.geometry.hydrostatic_ratio",GeometryHydrostaticRatioEnum,0.);
}
