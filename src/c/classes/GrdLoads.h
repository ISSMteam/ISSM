/*!\file GrdLoads.h
 * \brief: header file for GrdLoads
 */

#ifndef _SEALEVELGRDLOADS_H_
#define _SEALEVELGRDLOADS_H_

/*Headers:*/
#include "./SealevelGeometry.h"
#include "../toolkits/toolkits.h"

class GrdLoads{ 

	public: 

		int nactiveloads;
		int nactivesubloads[SLGEOM_NUMLOADS];
		Vector<IssmDouble>* vloads=NULL;
		IssmDouble*         loads=NULL;
		Vector<IssmDouble>* vsubloads[SLGEOM_NUMLOADS];
		IssmDouble*         subloads[SLGEOM_NUMLOADS];
		Vector<IssmDouble>* vsealevelloads=NULL;
		IssmDouble*         sealevelloads=NULL;
		Vector<IssmDouble>* vsubsealevelloads=NULL;
		IssmDouble*         subsealevelloads=NULL;
		int*         	    combined_loads_index=NULL;
		int*         	    combined_subloads_index[SLGEOM_NUMLOADS];
		IssmDouble*         combined_loads=NULL;
		IssmDouble*         combined_subloads[SLGEOM_NUMLOADS];

		GrdLoads(int nel, SealevelGeometry* slgeom);
		~GrdLoads();

		void AssembleSealevelLoads(void);
		void BroadcastLoads(void);
		void BroadcastSealevelLoads(void);
		void SHDegree2Coefficients(IssmDouble* deg2coeff, FemModel* femmodel, SealevelGeometry* slgeom);
		void Combineloads(int nel, SealevelGeometry* slgeom);
};
#endif  /* _SEALEVELGRDLOADS_H_ */
