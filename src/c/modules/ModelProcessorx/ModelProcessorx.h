/* \file ModelProcessorx.h
 * \brief  Header file for model processor
 */

#ifndef _MODEL_PROCESSORX_H_
#define _MODEL_PROCESSORX_H_

#include "../../classes/classes.h"
#include "../../analyses/analyses.h"

void ModelProcessorx(Elements** pelements, Nodes*** pnodes, Vertices** pvertices, Materials** pmaterials, Constraints*** pconstraints, Loads*** ploads, Parameters** pparameters,Inputs** pinputs,IoModel* iomodel,FILE* toolkitfile, char* rootpath,const int solution_type,const int nummodels,const int* analysis_type_listh);

/*Creation of fem datasets: general drivers*/
void CreateElements(Elements* elements,IoModel* iomodel,int nummodels);
void CreateMaterials(Elements* elements,Inputs* inputs,Materials* materials,IoModel* iomodel,int nummodels);
void CreateVertices(Elements* elements,Vertices* vertices,IoModel* iomodel,int solution_type,bool isamr=false);
void CreateParameters(Parameters*parameters,IoModel* iomodel,char* rootpath,FILE* toolkitfile,const int solution_type);
void CreateParametersAutodiff(Parameters* parameters,IoModel* iomodel);
void CreateParametersControl(Parameters* parameters,IoModel* iomodel,int solution_type);
void CreateParametersDakota(Parameters* parameters,IoModel* iomodel,char* rootpath);
void CreateOutputDefinitions(Elements* elements, Parameters* parameters,Inputs* inputs,IoModel* iomodel);
void UpdateElementsAndMaterialsControl(Elements* elements,Parameters* parameters,Inputs* inputs,Materials* materials, IoModel* iomodel);
void UpdateElementsAndMaterialsControlAD(Elements* elements,Parameters* parameters,Inputs* inputs,Materials* materials, IoModel* iomodel);
void UpdateElementsAndMaterialsDakota(Elements* elements,Inputs* inputs,Materials* materials, IoModel* iomodel);
void UpdateElementsTransient(Elements* elements,Parameters* parameters,Inputs* inputs,IoModel* iomodel);
void UpdateParametersTransient(Parameters* parameters,IoModel* iomodel);
void CreateNodes(Nodes*nodes, IoModel* iomodel,int analysis,int finite_element,bool isamr=false,int approximation=NoneApproximationEnum,int* approximations=NULL);

/*partitioning: */
void ElementsAndVerticesPartitioning(IoModel* iomodel);
void FacesPartitioning(IoModel* iomodel);
void EdgesPartitioning(IoModel* iomodel);

/*Mesh properties*/
void CreateEdges(IoModel* iomodel);
void CreateFaces(IoModel* iomodel);
void CreateFaces3d(IoModel* iomodel);
void EdgeOnBoundaryFlags(bool** pflags,IoModel* iomodel);

/*Connectivity*/
void CreateSingleNodeToElementConnectivity(IoModel* iomodel);
void CreateNumberNodeToElementConnectivity(IoModel* iomodel);
#endif
