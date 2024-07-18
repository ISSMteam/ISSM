/*!\file:  MeshPartitionx.h
 * \brief  header file for partitioning module.
 */ 

#ifndef _MESHPARTITIONX_H
#define _MESHPARTITIONX_H

int MeshPartitionx(int** pepart,int** pnpart,int numberofelements,int numberofnodes,int* elements,
		int numberofelements2d,int numberofnodes2d,int* elements2d,int* vweights,int numlayers,int elements_width, int meshelementtype,int num_procs);

#endif /* _MESHPARTITIONX_H */
