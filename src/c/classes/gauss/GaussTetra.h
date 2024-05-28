/*!\file GaussTetra.h
 * \brief: header file for node object
 */

#ifndef _GAUSSTETRA_H_
#define _GAUSSTETRA_H_

/*Headers:*/
#include "../../shared/Numerics/types.h"
#include "./Gauss.h"

class GaussTetra: public Gauss{

	private:
		int         numgauss;   /*Total number of gauss points*/
		int         ig;         /*Current gauss point index*/
		IssmDouble* weights;
		IssmDouble* coords1;
		IssmDouble* coords2;
		IssmDouble* coords3;
		IssmDouble* coords4;

	public:
		IssmDouble coord1;
		IssmDouble coord2;
		IssmDouble coord3;
		IssmDouble coord4;

	public:

		/*GaussTetra constructors, destructors*/
		GaussTetra();
		GaussTetra(int order);
		GaussTetra(int index1,int index2,int index3,int order);
		~GaussTetra();

		/*Methods*/
		bool next(void);
		void Echo(void);
		int  Enum(void);
		void GaussPoint(int ig);
		void GaussNode(int finitelement,int iv);
		void GaussVertex(int iv);
		void Reset(void);
		void SynchronizeGaussBase(Gauss* gauss);
};
#endif
