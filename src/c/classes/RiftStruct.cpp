/*!\file RiftStruct.c
 * \brief: implementation of the RiftStruct object
 */

#include "./classes.h"
#include "../shared/Enum/Enum.h"
#include "../shared/shared.h"

/*RiftStruct constructors and destructor*/
RiftStruct::RiftStruct(void){/*{{{*/

	this->numrifts             = 0;
	this->riftsnumsegments     = NULL;
	this->riftssegments        = NULL;
	this->riftsnumpairs        = NULL;
	this->riftspairs           = NULL;
	this->riftsnumpenaltypairs = NULL;
	this->riftspenaltypairs    = NULL;
	this->riftstips            = NULL;
	this->state                = NULL;

}/*}}}*/
RiftStruct::RiftStruct(int numrifts_in,int *riftsnumsegments_in,int**riftssegments_in,int *riftsnumpairs_in,int**riftspairs_in,int *riftsnumpenaltypairs_in,double **riftspenaltypairs_in,int * riftstips_in){/*{{{*/

	int i;

	/*numrifts*/
	this->numrifts = numrifts_in;
	if(!numrifts_in) return;

	/*riftsnumsegments*/
	_assert_(riftsnumsegments_in);
	this->riftsnumsegments=xNew<int>(numrifts_in);
	xMemCpy<int>(this->riftsnumsegments,riftsnumsegments_in,numrifts_in);

	/*riftssegments*/
	_assert_(riftssegments_in);
	this->riftssegments=xNew<int*>(numrifts_in);
	for(i=0;i<numrifts_in;i++){
		this->riftssegments[i]=xNew<int>(riftsnumsegments_in[i]*3);
		xMemCpy<int>(this->riftssegments[i],riftssegments_in[i],riftsnumsegments_in[i]*3);
	}

	/*riftsnumpairs*/
	_assert_(riftsnumpairs_in);
	this->riftsnumpairs=xNew<int>(numrifts_in);
	xMemCpy<int>(this->riftsnumpairs,riftsnumpairs_in,numrifts_in);

	/*riftspairs*/
	_assert_(riftspairs_in);
	this->riftspairs=xNew<int*>(numrifts_in);
	for(i=0;i<numrifts_in;i++){
		this->riftspairs[i]=xNew<int>(riftsnumpairs_in[i]*2);
		xMemCpy<int>(this->riftspairs[i],riftspairs_in[i],riftsnumpairs_in[i]*2);
	}

	/*riftsnumpenaltypairs*/
	_assert_(riftsnumpenaltypairs_in);
	this->riftsnumpenaltypairs=xNew<int>(numrifts_in);
	xMemCpy<int>(this->riftsnumpenaltypairs,riftsnumpenaltypairs_in,numrifts_in);

	/*riftspenaltypairs*/
	_assert_(riftspenaltypairs_in);
	this->riftspenaltypairs=xNew<double*>(numrifts_in);
	for(i=0;i<numrifts_in;i++){
		this->riftspenaltypairs[i]=xNew<double>(riftsnumpenaltypairs_in[i]*7);
		xMemCpy<double>(this->riftspenaltypairs[i],riftspenaltypairs_in[i],riftsnumpenaltypairs_in[i]*7);
	}

	/*riftstips*/
	_assert_(riftstips_in);
	this->riftstips=xNew<int>(2*numrifts_in);
	xMemCpy<int>(this->riftstips,riftstips_in,2*numrifts_in);

	/*state*/
	this->state=xNew<double*>(numrifts_in);
	for(i=0;i<numrifts_in;i++){
		this->state[i]=xNew<double>(riftsnumpenaltypairs_in[i]);
		for(int j=0;j<riftsnumpenaltypairs_in[i];j++) (this->state[i])[j]=FreeEnum;
	}

}/*}}}*/
RiftStruct::~RiftStruct(void){/*{{{*/

	xDelete<int>(this->riftsnumsegments);
	xDelete<int*>(this->riftssegments);
	xDelete<int>(this->riftsnumpairs);
	xDelete<int*>(this->riftspairs);
	xDelete<int>(this->riftsnumpenaltypairs);
	xDelete<double*>(this->riftspenaltypairs);
	xDelete<int>(this->riftstips);
	xDelete<double*>(this->state);

}/*}}}*/
