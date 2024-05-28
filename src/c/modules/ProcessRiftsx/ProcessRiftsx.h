/*!\file:  ProcessRiftsx.h
 * \brief header file for ProcessRifts module
 */ 

#ifndef _PROCESSRIFTX_H
#define _PROCESSRIFTX_H

class RiftStruct;

void ProcessRiftsx(int** pindex,int* pnel,double** px,double** py,int* pnods,int** psegments,int** psegmentmarkers,int *pnum_seg,RiftStruct **priftstruct);

#endif  /* _PROCESSRIFTX_H*/
