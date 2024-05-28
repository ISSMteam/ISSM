/*!\file GaussPenta.h
 * \brief: header file for node object
 */

#ifndef _GAUSSPENTA_H_
#define _GAUSSPENTA_H_

/*Headers:*/
#include "../../shared/Numerics/types.h"
#include "./Gauss.h"
class GaussTria;

class GaussPenta: public Gauss{

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

		/*GaussPenta constructors, destructors*/
		GaussPenta();
		GaussPenta(int order_horiz,int order_vert);
		GaussPenta(int index1, int index2,int order);
		GaussPenta(int index1, int index2, int index3, int order);
		GaussPenta(int index1, int index2, int index3, int index4,int order_horiz,int order_vert);
		GaussPenta(int index,IssmDouble r1, IssmDouble r2,bool maintlyfloating,int order);
		GaussPenta(IssmDouble area_coordinates[4][3],int order_horiz,int order_vert);
		GaussPenta(IssmDouble area_coordinates[2][3],int order_horiz);
		~GaussPenta();

		/*Methods*/
		bool next(void);
		void Echo(void);
		int  Enum(void);
		void GaussFaceTria(int index1, int index2, int index3, int order);
		void GaussNode(int finitelement,int iv);
		void GaussPoint(int ig);
		void GaussVertex(int iv);
		void Reset(void);
		void SynchronizeGaussBase(Gauss* gauss);
};
#endif
