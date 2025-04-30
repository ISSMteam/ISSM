/*
 * \file GrdLoads.cpp
 * \brief: Implementation of GrdLoads class
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./GrdLoads.h"
#include "./SealevelGeometry.h"
using namespace std;
/*}}}*/

/*Object constructors and destructor*/
GrdLoads::GrdLoads(int nel,SealevelGeometry* slgeom){ /*{{{*/

	_assert_(slgeom);

	/*allocate:*/
	this->nactiveloads=0;
	this->vloads=new Vector<IssmDouble>(nel);
	for (int i=0;i<SLGEOM_NUMLOADS;i++) {
		this->vsubloads[i] = new Vector<IssmDouble>(slgeom->nbar[i]);
	}

	this->vsealevelloads=new Vector<IssmDouble>(nel);
	this->vsealevelloads->Set(0);vsealevelloads->Assemble();

	this->vsubsealevelloads=new Vector<IssmDouble>(slgeom->nbar[SLGEOM_OCEAN]);

	this->combined_loads=NULL;
	this->combined_loads_index=NULL;
	for (int i=0;i<SLGEOM_NUMLOADS;i++) {
		this->nactivesubloads[i]=0;
		this->combined_subloads[i]=NULL;
		this->combined_subloads_index[i]=NULL;
	}

	/*make sure all pointers that are not allocated are NULL:*/
	this->loads=NULL;
	for(int i=0;i<SLGEOM_NUMLOADS;i++) this->subloads[i]=NULL;
	this->sealevelloads=NULL;
	this->subsealevelloads=NULL;

}; /*}}}*/
GrdLoads::~GrdLoads(){ /*{{{*/

	delete vloads;
	xDelete<IssmDouble>(loads);
	delete vsealevelloads;
	xDelete<IssmDouble>(sealevelloads);
	delete vsubsealevelloads;
	xDelete<IssmDouble>(subsealevelloads);
	if (combined_loads) xDelete<IssmDouble>(combined_loads);
	if (combined_loads_index) xDelete<int>(combined_loads_index);
	for(int i=0;i<SLGEOM_NUMLOADS;i++){
		delete vsubloads[i];
		xDelete<IssmDouble>(subloads[i]);
		if (combined_subloads[i]) xDelete<IssmDouble>(combined_subloads[i]);
		if (combined_subloads_index[i]) xDelete<int>(combined_subloads_index[i]);
	}

}; /*}}}*/

void GrdLoads::BroadcastLoads(void){ /*{{{*/

	/*Initialize barycentre vectors, now that we know their size: */
	this->vloads->Assemble();
	for(int i=0;i<SLGEOM_NUMLOADS;i++){
		vsubloads[i]->Assemble();
	}

	/*Avoid leaks:*/
	if(loads) xDelete<IssmDouble>(loads);
	for(int i=0;i<SLGEOM_NUMLOADS;i++){
		if(subloads[i])xDelete<IssmDouble>(subloads[i]);
	}

	/*Serialize:*/
	loads=vloads->ToMPISerial();
	for(int i=0;i<SLGEOM_NUMLOADS;i++){
		subloads[i]=vsubloads[i]->ToMPISerial();
	}

} /*}}}*/
void GrdLoads::AssembleSealevelLoads(void){ /*{{{*/

	vsealevelloads->Assemble();
	vsubsealevelloads->Assemble();

} /*}}}*/
void GrdLoads::BroadcastSealevelLoads(void){ /*{{{*/

	/*Avoid leakds:*/
	if(sealevelloads)    xDelete<IssmDouble>(sealevelloads);
	if(subsealevelloads) xDelete<IssmDouble>(subsealevelloads);

	/*Serialize:*/
	sealevelloads    = vsealevelloads->ToMPISerial();
	subsealevelloads = vsubsealevelloads->ToMPISerial();

} /*}}}*/
void GrdLoads::SHDegree2Coefficients(IssmDouble* deg2coeff, FemModel* femmodel, SealevelGeometry* slgeom){ /*{{{*/

	IssmDouble re, S;
	int ylmindex, intj;
	IssmDouble deg2coeff_local[5];
	//IssmDouble area;

	femmodel->parameters->FindParam(&re,SolidearthPlanetRadiusEnum);

	for(int c=0;c<5;c++) deg2coeff_local[c]=0;

	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		//first, loads on centroids
		S=0;
		S+=loads[element->Sid()];

		if(sealevelloads) S+=sealevelloads[element->Sid()];
		if(S!=0){

			for (int c=0;c<5;c++){ //degree l=2 has 2*l+1=5 coefficients: 2,0; 2,1cos; 2,1sin; 2,2cos; 2,2sin
				ylmindex=(4+c)*slgeom->localnel+element->lid; // starting at index=l^2
				deg2coeff_local[c] += S/re/re*slgeom->Ylm[ylmindex];
			}
		}
		//add loads on subelement barycenters
		for (int i=0;i<SLGEOM_NUMLOADS;i++){
			if (slgeom->issubelement[i][element->lid]){
				intj=slgeom->subelementmapping[i][element->lid];
				S=0;
				S+=subloads[i][intj];
				if(i==SLGEOM_OCEAN && sealevelloads) S+=subsealevelloads[intj];
				if(S!=0){
					//area=slgeom->area_subel[i][intj];
					for (int c=0;c<5;c++){ //degree l=2 has 2*l+1=5 coefficients
						ylmindex=(4+c)*slgeom->localnel+element->lid; // starting at index=l^2
						deg2coeff_local[c] += S/re/re*slgeom->Ylm_subel[i][ylmindex];
					}
				}
			}
		}
	}

	//sum each degree 2 coefficient across all cpus to get global total
	ISSM_MPI_Reduce (&deg2coeff_local[0],&deg2coeff[0],1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&deg2coeff[0],1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	ISSM_MPI_Reduce (&deg2coeff_local[1],&deg2coeff[1],1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&deg2coeff[1],1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	ISSM_MPI_Reduce (&deg2coeff_local[2],&deg2coeff[2],1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&deg2coeff[2],1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

}; /*}}}*/ 
void GrdLoads::Combineloads(int nel,SealevelGeometry* slgeom){ /*{{{*/

	int e,l, nbar, ae;
	//Determine loads /*{{{*/
	nactiveloads=0;

	if (combined_loads) xDelete<IssmDouble>(combined_loads);
	if (combined_loads_index) xDelete<int>(combined_loads_index);

	//find non zero centroid loads, combine with sealevelloads
	if(sealevelloads){
		for (e=0;e<nel;e++){
			if (loads[e]+sealevelloads[e]!=0) nactiveloads++;
		}
	}
	else { 
		for (e=0;e<nel;e++){
			if (loads[e]!=0) nactiveloads++;
		}
	}

	combined_loads=xNewZeroInit<IssmDouble>(nactiveloads);
	combined_loads_index=xNewZeroInit<int>(nactiveloads);

	ae=0;
	if(sealevelloads){
		for (e=0;e<nel;e++){
			if (loads[e]+sealevelloads[e]!=0){
				combined_loads[ae]=loads[e]+sealevelloads[e];
				combined_loads_index[ae]=e;
				ae++;
			}
		}
	}
	else { 
		for (e=0;e<nel;e++){
			if (loads[e]!=0){
				combined_loads[ae]=loads[e];
				combined_loads_index[ae]=e;
				ae++;			
			}
		}
	}

	//subloads
	for(l=0;l<SLGEOM_NUMLOADS;l++){
		nactivesubloads[l]=0;
		nbar=slgeom->nbar[l];
		if (combined_subloads[l]) xDelete<IssmDouble>(combined_subloads[l]);
		if (combined_subloads_index[l]) xDelete<int>(combined_subloads_index[l]);

		//find non zero subelement loads, combine with subsealevelloads
		if(subsealevelloads && l==SLGEOM_OCEAN){
			for (e=0;e<nbar;e++){
				if (subloads[l][e]+subsealevelloads[e]!=0) nactivesubloads[l]++;
			}
		}
		else { 
			for (e=0;e<nbar;e++){
				if (subloads[l][e]!=0) nactivesubloads[l]++;;
			}
		}

		combined_subloads[l]=xNewZeroInit<IssmDouble>(nactivesubloads[l]);
		combined_subloads_index[l]=xNewZeroInit<int>(nactivesubloads[l]);

		ae=0;
		if(subsealevelloads && l==SLGEOM_OCEAN){
			for (e=0;e<nbar;e++){
				if (subloads[l][e]+sealevelloads[e]!=0){
					combined_subloads[l][ae]=subloads[l][e]+subsealevelloads[e];
					combined_subloads_index[l][ae]=e;
					ae++;
				}
			}
		}
		else { 
			for (e=0;e<nbar;e++){
				if (subloads[l][e]!=0){
					combined_subloads[l][ae]=subloads[l][e];
					combined_subloads_index[l][ae]=e;
					ae++;			
				}
			}
		}

		/*for(l=0;l<SLGEOM_NUMLOADS;l++){
			nbar=slgeom->nbar[l];
			for (e=0;e<nbar;e++){
				combined_subloads[l][e]=(subloads[l][e]);
			}
		}
		if(sealevelloads){
			nbar=slgeom->nbar[SLGEOM_OCEAN];
			for (e=0;e<nbar;e++){
				combined_subloads[SLGEOM_OCEAN][e]+=(subsealevelloads[e]);
			}
		}*/
	}
}; /*}}}*/
