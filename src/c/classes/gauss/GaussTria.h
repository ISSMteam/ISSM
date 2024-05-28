/*!\file GaussTria.h
 * \brief: header file for node object
 */

#ifndef _GAUSSTRIA_H_
#define _GAUSSTRIA_H_

/*Headers:*/
#include "../../shared/Numerics/types.h"
#include "./Gauss.h"

class GaussTria: public Gauss{

	private:
		int         numgauss;   /*Total number of gauss points*/
		int         ig;         /*Current gauss point index*/
		IssmDouble* weights;
		IssmDouble* coords1;
		IssmDouble* coords2;
		IssmDouble* coords3;

	public:
		IssmDouble coord1;
		IssmDouble coord2;
		IssmDouble coord3;

	public:

		/*GaussTria constructors, destructors*/
		GaussTria();
		GaussTria(int order);
		GaussTria(int index1,int index2,int order);
		GaussTria(int index,IssmDouble r1, IssmDouble r2,bool maintlyfloating,int order);
		GaussTria(int index,IssmDouble r1, IssmDouble r2,int order);
		GaussTria(IssmDouble r1, IssmDouble r2,int order);
		GaussTria(IssmDouble area_coordinates[2][3],int order);
		~GaussTria();

		/*Methods*/
		bool next(void);
		void Echo(void);
		int  Enum(void);
		void GaussEdgeCenter(int index1,int index2);
		void GaussFromCoords(IssmDouble x1,IssmDouble y1,IssmDouble* xyz_list);
		void GaussPoint(int ig);
		void GaussNode(int finitelement,int iv);
		void GaussVertex(int iv);
		void Reset(void);
		void SynchronizeGaussBase(Gauss* gauss);
};
#endif  /* _GAUSSTRIA_H_ */
