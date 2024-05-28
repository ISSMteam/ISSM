/*!\file Gauss.h
 * \brief: header file for vvirtual Gauss object
 */

#ifndef _GAUSS_H_
#define _GAUSS_H_

class Gauss{

	public: 
		IssmDouble   weight;

		virtual      ~Gauss(){};
		virtual void Echo(void)=0;
		virtual void Reset(void)=0;
		virtual bool next(void)=0;
		virtual int  Enum(void)=0;
		virtual void GaussNode(int finitelement,int iv)=0;
		virtual void GaussPoint(int ig)=0;
		virtual void GaussVertex(int iv)=0;
		virtual void SynchronizeGaussBase(Gauss* gauss)=0;

};
#endif
