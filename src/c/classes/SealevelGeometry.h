/*!\file SealevelGeometry.h
 * \brief: header file for SealevelMask geometry
 */

#ifndef _SEALEVELGEOMETRY_H_
#define _SEALEVELGEOMETRY_H_

/*Headers:*/
#define SLGEOM_NUMLOADS 3
#define SLGEOM_OCEAN 0 
#define SLGEOM_ICE 1 
#define SLGEOM_WATER 2
#define SLMAXVERTICES 3

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"

class SealevelGeometry{ 

	public: 

		int         localnel;
		IssmDouble* LoadWeigths[SLGEOM_NUMLOADS][SLMAXVERTICES];
		IssmDouble* LoadArea[SLGEOM_NUMLOADS];
		Vector<IssmDouble>* vlatbarycentre[SLGEOM_NUMLOADS];
		Vector<IssmDouble>* vlongbarycentre[SLGEOM_NUMLOADS];
		Vector<IssmDouble>* vareae_subel[SLGEOM_NUMLOADS];
		IssmDouble* latbarycentre[SLGEOM_NUMLOADS];
		IssmDouble* longbarycentre[SLGEOM_NUMLOADS];
		IssmDouble* area_subel[SLGEOM_NUMLOADS];
		IssmDouble* late;
		IssmDouble* longe;
		IssmDouble* Ylm;
		IssmDouble* Ylm_subel[SLGEOM_NUMLOADS];
		IssmDouble* YlmNorm[9];
		bool* isoceanin;
		bool*       issubelement[SLGEOM_NUMLOADS]; 
		int*        subelementmapping[SLGEOM_NUMLOADS];
		int         nsubel[SLGEOM_NUMLOADS];
		int         nbar[SLGEOM_NUMLOADS];
		int*        lids; 

		SealevelGeometry(int localnel,int localnods);
		~SealevelGeometry();
		void InitializeMappingsAndBarycentres(void);
		void Assemble(void);
		int AlphaIndexEnum(int l);
		int AzimuthIndexEnum(int l);
		void BuildSphericalHarmonics(void);
};
#endif  /* _SEALEVELGEOMETRY_H_ */
