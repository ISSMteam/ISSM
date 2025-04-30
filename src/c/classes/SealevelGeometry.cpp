/*
 * \file SealevelGeometry.cpp
 * \brief: Implementation of SealevelGeometry class
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./SealevelGeometry.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
SealevelGeometry::SealevelGeometry(int localnelin,int localnodsin){ /*{{{*/

	localnel=localnelin;

	for(int i=0;i<SLGEOM_NUMLOADS;i++){
		for (int j=0;j<SLMAXVERTICES;j++){
			LoadWeigths[i][j] = xNewZeroInit<IssmDouble>(localnel);
		}
		vlatbarycentre[i]  = NULL; //we don't know yet
		vlongbarycentre[i] = NULL;
		vareae_subel[i]    = NULL;
		latbarycentre[i]   = NULL; //we don't know yet
		longbarycentre[i]  = NULL;
		area_subel[i]      = NULL;

		LoadArea[i]      = xNewZeroInit<IssmDouble>(localnel);
		issubelement[i]  = xNewZeroInit<bool>(localnel);
      Ylm_subel[i]     = xNewZeroInit<IssmDouble>(localnel*9);
		subelementmapping[i] = NULL;
		nsubel[i] = 0;
		nbar[i]   = 0;
	}
	late      = xNew<IssmDouble>(localnel);
	longe     = xNew<IssmDouble>(localnel);
	isoceanin = xNew<bool>(localnel);
	lids      = xNew<int>(localnodsin);
	Ylm       = xNewZeroInit<IssmDouble>(localnel*9); // (degmax+1)^2 terms, degmax = 2

}; /*}}}*/
SealevelGeometry::~SealevelGeometry(){ /*{{{*/
	for(int i=0;i<SLGEOM_NUMLOADS;i++){
		for (int j=0;j<SLMAXVERTICES;j++){
			xDelete<IssmDouble>(LoadWeigths[i][j]);
		}
		xDelete<IssmDouble>(LoadArea[i]);
		xDelete<bool>(issubelement[i]);
		xDelete<int>(subelementmapping[i]);
		delete  vlatbarycentre[i];
		delete  vlongbarycentre[i];
		delete  vareae_subel[i];
		xDelete<IssmDouble>(latbarycentre[i]);
		xDelete<IssmDouble>(longbarycentre[i]);
		xDelete<IssmDouble>(area_subel[i]);
		xDelete<IssmDouble>(Ylm_subel[i]);
	}
	xDelete<IssmDouble>(Ylm);
	xDelete<IssmDouble>(late);
	xDelete<IssmDouble>(longe);
	xDelete<bool>(isoceanin);
	xDelete<int>(lids);
}; /*}}}*/

void SealevelGeometry::InitializeMappingsAndBarycentres(void){ /*{{{*/

	int dummy;
	bool fromlocalsize=true;
	int lower_row;

	for (int i=0;i<SLGEOM_NUMLOADS;i++){
		subelementmapping[i]=xNew<int>(localnel);
		#ifdef _HAVE_MPI_
		GetOwnershipBoundariesFromRange(&lower_row,&dummy,nsubel[i],IssmComm::GetComm());
		#else
		_error_("not supported without MIP ");
		#endif

		int count=0;
		for (int j=0;j<localnel;j++){
			if(issubelement[i][j]){
				subelementmapping[i][j]=lower_row+count;
				count++;
			}
		}
	}

	/*Initialize barycentre vectors, now that we know their size: */
	for (int i=0;i<SLGEOM_NUMLOADS;i++){
		vlatbarycentre[i]=new Vector<IssmDouble>(nsubel[i],fromlocalsize);
		vlongbarycentre[i]=new Vector<IssmDouble>(nsubel[i],fromlocalsize);
		vareae_subel[i]=new Vector<IssmDouble>(nsubel[i],fromlocalsize);
		vlatbarycentre[i]->GetSize(&nbar[i]);
	}

} /*}}}*/
void SealevelGeometry::Assemble(void){ /*{{{*/

	/*Initialize barycentre vectors, now that we know their size: */
	for (int i=0;i<SLGEOM_NUMLOADS;i++){
		vlatbarycentre[i]->Assemble();
		vlongbarycentre[i]->Assemble();
		vareae_subel[i]->Assemble();

		latbarycentre[i]=vlatbarycentre[i]->ToMPISerial();
		longbarycentre[i]=vlongbarycentre[i]->ToMPISerial();
		area_subel[i]=vareae_subel[i]->ToMPISerial();
	}

	/*Also, we'll need the barycentre associated areas:*/

} /*}}}*/
int SealevelGeometry::AlphaIndexEnum(int l){ /*{{{*/

	int output = -1;
	switch(l){
		case SLGEOM_OCEAN: output=SealevelchangeAlphaIndexOceanEnum; break;
		case SLGEOM_ICE:   output=SealevelchangeAlphaIndexIceEnum;   break;
		case SLGEOM_WATER: output=SealevelchangeAlphaIndexHydroEnum; break;
		default: _error_("not supported");
	}
	return output;

} /*}}}*/
int SealevelGeometry::AzimuthIndexEnum(int l){ /*{{{*/

	int output = -1;
	switch(l){
		case SLGEOM_OCEAN: output=SealevelchangeAzimuthIndexOceanEnum; break;
		case SLGEOM_ICE:   output=SealevelchangeAzimuthIndexIceEnum;   break;
		case SLGEOM_WATER: output=SealevelchangeAzimuthIndexHydroEnum; break;
		default: _error_("not supported");
	}
	return output;

} /*}}}*/
void SealevelGeometry::BuildSphericalHarmonics(){ /*{{{*/
	//builds spherical harmonics functions for degrees 0, 1, 2 on centroids/barycenters
	//0: used for global average
	//1: used for geocenter motion
	//2: used for rotational feedback
	int intj, count;

	IssmDouble YlmNorm[9];

	//YlmNormalization: N^2=(2*l+1)/4/pi * factorial(l-m)/factorial(l+m) if m==0
	//             : 2*N^2 if m>0
	// such that integral(Ylm * Ylm *YlmNorm dS) = 1 on the unit sphere.
	YlmNorm[0] = (0.25/M_PI);     // Y00
	YlmNorm[1] = (0.75/M_PI);     // Y10
	YlmNorm[2] = (0.75/M_PI);     // Y11c
	YlmNorm[3] = YlmNorm[2];      // Y11s
	YlmNorm[4] = (1.25/M_PI);     // Y20
	YlmNorm[5] = (1.25/3./M_PI);  // Y21c
	YlmNorm[6] = YlmNorm[5];      // Y21s
	YlmNorm[7] = (1.25/12./M_PI); // Y22c
	YlmNorm[8] = YlmNorm[7];      // Y22s

	for (int e=0;e<localnel;e++){
		IssmDouble lat = late[e]*M_PI/180.;
		IssmDouble lon = longe[e]*M_PI/180.;
		Ylm[0*localnel+e] = 1.0 *YlmNorm[0]; //Y00
		Ylm[1*localnel+e] = sin(lat)*YlmNorm[1]; //Y10
		Ylm[2*localnel+e] = cos(lat)*cos(lon)*YlmNorm[2]; //Y11cos
		Ylm[3*localnel+e] = cos(lat)*sin(lon)*YlmNorm[3]; //Y11sin

		//Ylm[4*localnel+e] = 0.25 - 0.75*cos(2.0*lat) ; //Y20
		Ylm[4*localnel+e] = (1.5*pow(sin(lat),2.)-0.5)*YlmNorm[4]; //Y20
		Ylm[5*localnel+e] = 1.5*sin(2.*lat)*cos(lon)*YlmNorm[5]; //Y21cos
		Ylm[6*localnel+e] = 1.5*sin(2.*lat)*sin(lon)*YlmNorm[6]; //Y21sin
		Ylm[7*localnel+e] = 1.5*(1.+cos(2.*lat))*cos(2.*lon)*YlmNorm[7]; //Y22cos
		Ylm[8*localnel+e] = 1.5*(1.+cos(2.*lat))*sin(2.*lon)*YlmNorm[8]; //Y22sin
	}

	for (int i=0;i<SLGEOM_NUMLOADS;i++){
		for (int e=0;e<localnel;e++){
			if (issubelement[i][e]){
				intj=subelementmapping[i][e];
				IssmDouble lat=latbarycentre[i][intj]*M_PI/180.;
				IssmDouble lon=longbarycentre[i][intj]*M_PI/180.;
				Ylm_subel[i][0*localnel+e] = 1.0*YlmNorm[0]; //Y00

				Ylm_subel[i][1*localnel+e] = sin(lat)*YlmNorm[1]; //Y10
				Ylm_subel[i][2*localnel+e] = cos(lat)*cos(lon)*YlmNorm[2]; //Y11cos
				Ylm_subel[i][3*localnel+e] = cos(lat)*sin(lon)*YlmNorm[3]; //Y11sin

				Ylm_subel[i][4*localnel+e] = (1.5*pow(sin(lat),2.)-0.5)*YlmNorm[4]; //Y20
				Ylm_subel[i][5*localnel+e] = 1.5*sin(2.*lat)*cos(lon)*YlmNorm[5]; //Y21cos
				Ylm_subel[i][6*localnel+e] = 1.5*sin(2.*lat)*sin(lon)*YlmNorm[6]; //Y21sin
				Ylm_subel[i][7*localnel+e] = 1.5*(1.+cos(2.*lat))*cos(2.*lon)*YlmNorm[7]; //Y22cos
				Ylm_subel[i][8*localnel+e] = 1.5*(1.+cos(2.*lat))*sin(2.*lon)*YlmNorm[8]; //Y22sin
			}
		}
	}
} /*}}}*/
