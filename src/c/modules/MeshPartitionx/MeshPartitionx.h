/*!\file:  MeshPartitionx.h
 * \brief  header file for partitioning module.
 */ 

#ifndef _MESHPARTITIONX_H
#define _MESHPARTITIONX_H

#include "../../shared/shared.h"

/* local prototypes: */
template <class doubletype> 
int MeshPartitionx(int** pepart,int** pnpart,int numberofelements,int numberofnodes,int* elements,
		int numberofelements2d,int numberofnodes2d,doubletype* elements2d,int* vweights,int numlayers,int elements_width, int meshelementtype,int num_procs){

	int noerr=1;
	int i,j;

	/*Metis partitioning: */
	int* epart=NULL;
	int* npart=NULL;
	int* index=NULL;

	int* epart2d=NULL;
	int* npart2d=NULL;
	int* index2d=NULL;
	int  count=0;

	switch(meshelementtype){
		case TriaEnum:
		case TetraEnum:
			epart=xNew<int>(numberofelements);
			npart=xNew<int>(numberofnodes);
			index=xNew<int>(elements_width*numberofelements);
			for (i=0;i<numberofelements;i++){
				for (j=0;j<elements_width;j++){
					*(index+elements_width*i+j)=(*(elements+elements_width*i+j))-1; //-1 for C indexing in Metis
				}
			}

			/*Partition using Metis:*/
			if (num_procs>1){
#ifdef _HAVE_METIS_
				METIS_PartMeshNodalPatch(numberofelements,numberofnodes,index,vweights,num_procs,epart, npart);
#else
				_error_("metis has not beed installed. Cannot run with more than 1 cpu");
#endif
			}
			else if (num_procs==1){
				/*METIS does not know how to deal with one cpu only!*/
				for (i=0;i<numberofelements;i++) epart[i]=0;
				for (i=0;i<numberofnodes;i++)    npart[i]=0;
			}
			else _error_("At least one processor is required");
			break;
		case PentaEnum:
			/*We have a 3d mesh, made of a regularly extruded 2d mesh. We first partition the 2d mesh, then we extrude the partition: */

			/*First build concatenated 2d mesh  from 2d_coll and 2d_noncoll: */
			epart2d=xNew<int>(numberofelements2d);
			npart2d=xNew<int>(numberofnodes2d); 
			index2d=xNew<int>(3*numberofelements2d);

			for (i=0;i<numberofelements2d;i++){
				for (j=0;j<3;j++){
					*(index2d+3*i+j)=reCast<int>(*(elements2d+3*i+j))-1; //-1 for C indexing in Metis
				}
			}

			/*Partition using Metis:*/
			if (num_procs>1){
#ifdef _HAVE_METIS_
				METIS_PartMeshNodalPatch(numberofelements2d,numberofnodes2d,index2d,vweights,num_procs,epart2d,npart2d);
#else
				_error_("metis has not beed installed. Cannot run with more than 1 cpu");
#endif
			}
			else if (num_procs==1){
				/*METIS does not know how to deal with one cpu only!*/
				for (i=0;i<numberofelements2d;i++) epart2d[i]=0;
				for (i=0;i<numberofnodes2d;i++)    npart2d[i]=0;
			}
			else _error_("At least one processor is required");

			/*Extrude epart2d to epart, using numlayers: */
			epart=xNew<int>(numberofelements);

			count=0;
			for(i=0;i<(numlayers-1);i++){
				for(j=0;j<numberofelements2d;j++){
					epart[count]=epart2d[j];
					count++;
				}
			}

			/*Extrude npart2d to npart, using numlayers: */
			npart=xNew<int>(numberofnodes);

			count=0;
			for(i=0;i<(numlayers);i++){
				for(j=0;j<numberofnodes2d;j++){
					npart[count]=npart2d[j];
					count++;
				}
			}
			break;
		default:
			_error_("mesh type "<<EnumToStringx(meshelementtype)<<" not supported yet");
	}

	/*Assign output pointer:*/
	*pepart=epart;
	*pnpart=npart;

	/*Free resources: */
	xDelete<int>(index);
	xDelete<int>(epart2d);
	xDelete<int>(npart2d);
	xDelete<int>(index2d);
	return noerr;
}
#endif /* _MESHPARTITIONX_H */
