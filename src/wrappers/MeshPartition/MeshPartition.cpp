/*!\file:  MeshPartition.cpp
 * \brief: partition mesh according to number of areas, using Metis library.
*/

#include "./MeshPartition.h"

void MeshPartitionUsage(void){/*{{{*/
	_printf_("   usage:\n");
	_printf_("   [element_partitioning,node_partitioning]=MeshPartition(md.mesh,numpartitions)");
	_printf_("   where:\n");
	_printf_("      element_partitioning is a vector of partitioning area numbers, for every element.\n");
	_printf_("      node_partitioning is a vector of partitioning area numbers, for every node.\n");
	_printf_("\n");
}/*}}}*/
WRAPPER(MeshPartition_python){

	/* required input: */
	int  numberofelements;
	int  numberofvertices;
	int *elements = NULL;
	int  elements_width;
	int numberofelements2d=0;
	int numberofvertices2d=0;
	int* elements2d=NULL;
	int numberoflayers;
	int numareas=1;

	/* output: */
	int    *int_element_partitioning = NULL;
	int    *int_node_partitioning    = NULL;
	double *element_partitioning     = NULL;
	double *node_partitioning        = NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	CHECKARGUMENTS(NLHS,NRHS,&MeshPartitionUsage);

	/*Fetch data: */
	FetchData(&numberofvertices,NUMBEROFVERTICES);
	FetchData(&elements,&numberofelements,&elements_width,ELEMENTS);
	FetchData(&numberofvertices2d,NUMBEROFVERTICES2D);
	FetchData(&elements2d,&numberofelements2d,NULL,ELEMENTS2D);
	FetchData(&numberoflayers,NUMBEROFLAYERS);
	FetchData(&numareas,NUMAREAS);

	/*Get mesh element type and convert to Enum*/
	char* meshtype_str = NULL;
	FetchData(&meshtype_str,MESHELEMENTTYPE);
	int meshelementtype = StringToEnumx(meshtype_str); 
	xDelete<char>(meshtype_str);

	/*Run partitioning algorithm based on a "clever" use of the Metis partitioner: */
	MeshPartitionx(&int_element_partitioning,&int_node_partitioning,numberofelements,numberofvertices,elements,
		numberofelements2d,numberofvertices2d,elements2d,NULL,numberoflayers,elements_width,meshelementtype,numareas);

	/*Post process node_partitioning and element_partitioning to be in double format. Metis needed them in int* format: */
	element_partitioning=xNew<double>(numberofelements);
	for(int i=0;i<numberofelements;i++){
		element_partitioning[i]=(double)int_element_partitioning[i]+1; //Metis indexing from 0, matlab from 1.
	}

	node_partitioning=xNew<double>(numberofvertices);
	for(int i=0;i<numberofvertices;i++){
		node_partitioning[i]=(double)int_node_partitioning[i]+1; //Metis indexing from 0, matlab from 1.
	}

	/*Write data:*/
	WriteData(ELEMENTPARTITIONING,element_partitioning,numberofelements);
	WriteData(NODEPARTITIONING,node_partitioning,numberofvertices);

	/*Free resources:*/
	xDelete<int>(elements);
	xDelete<int>(elements2d);
	xDelete<int>(int_element_partitioning);
	xDelete<int>(int_node_partitioning);
	xDelete<double>(element_partitioning);
	xDelete<double>(node_partitioning);

	/*end module: */
	MODULEEND();
}
