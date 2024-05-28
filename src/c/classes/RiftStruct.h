/*!\file:  RiftStruct.h
 * \brief place holder for Rift Structure in mex module
 */ 

#ifndef _RIFTSTRUCT_H_
#define _RIFTSTRUCT_H_

class RiftStruct{

	public:
		int      numrifts;
		int    **riftspairs;
		double **riftspenaltypairs;
		int     *riftsnumpairs;
		int     *riftsnumpenaltypairs;
		int     *riftsnumsegments;
		int    **riftssegments;
		int     *riftstips;
		double **state;

		RiftStruct();
		RiftStruct(int numrifts_in,int *riftsnumsegments_in,int **riftssegments_in,int *riftsnumpairs_in,int **riftspairs_in,int *riftsnumpenaltypairs_in,double **riftspenaltypairs_in,int* riftstips_in);
		~RiftStruct();
};

#endif
