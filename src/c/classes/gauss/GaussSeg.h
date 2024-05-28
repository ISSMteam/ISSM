/*!\file GaussSeg.h
 * \brief: header file for node object
 */

#ifndef _GAUSSSEG_H_
#define _GAUSSSEG_H_

/*Headers:*/
#include "../../shared/Numerics/types.h"
#include "./Gauss.h"

class GaussSeg: public Gauss{

	private:
		int         numgauss;   /*Total number of gauss points*/
		int         ig;         /*Current gauss point index*/
		IssmDouble* weights;    /*List of weights*/
		IssmDouble* coords1;

	public:
		IssmDouble coord1;

	public:

		/*GaussSeg constructors, destructors*/
		GaussSeg();
		GaussSeg(int order);
		GaussSeg(IssmDouble position);
		~GaussSeg();

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
