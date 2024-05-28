/*!\file:  Material.h
 * \brief abstract class for Material object
 */ 

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

/*Headers:*/
/*{{{*/
class Inputs;
template <class doubletype> class Vector;
#include "../../datastructures/datastructures.h"
#include "../../toolkits/toolkits.h"
class Element;
class Elements;
class Gauss;
class Input;
class Input;
/*}}}*/

class Material: public Object{

	public: 
		virtual ~Material(){};

		/*Numerics*/
		virtual void       Configure(Elements* elements)=0;
		virtual Material*  copy2(Element* element)=0;
		virtual IssmDouble GetA(Gauss* gauss)=0;
		virtual IssmDouble GetAbar(Gauss* gauss)=0;
		virtual IssmDouble GetB(Gauss* gauss)=0;
		virtual IssmDouble GetBbar(Gauss* gauss)=0;
		virtual IssmDouble GetD(Gauss* gauss)=0;
		virtual IssmDouble GetDbar(Gauss* gauss)=0;
		virtual IssmDouble GetN()=0;
		virtual void       GetViscosity(IssmDouble* pviscosity,IssmDouble epseff,Gauss* gauss)=0;
		virtual void       GetViscosityBar(IssmDouble* pviscosity,IssmDouble epseff,Gauss* gauss)=0;
		virtual void       GetViscosityComplement(IssmDouble* pviscosity_complement, IssmDouble* pepsilon,Gauss* gauss)=0;
		virtual void       GetViscosityDComplement(IssmDouble* pviscosity_complement, IssmDouble* pepsilon,Gauss* gauss)=0;
		virtual void       GetViscosityDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon,Gauss* gauss)=0;
		virtual void       GetViscosity_B(IssmDouble* pviscosity,IssmDouble epseff,Gauss* gauss)=0;
		virtual void       GetViscosity_D(IssmDouble* pviscosity,IssmDouble epseff,Gauss* gauss)=0;
		virtual void       GetViscosity2dDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon,Gauss* gauss)=0;
		virtual bool       IsDamage()=0;
		virtual bool       IsEnhanced()=0;
		virtual void       ResetHooks()=0;

		virtual void       ViscosityFSDerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss)=0;
		virtual void       ViscosityHODerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss)=0;
		virtual void       ViscositySSADerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss)=0;
		virtual void       ViscosityFS(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input)=0;
		virtual void       ViscosityHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input)=0;
		virtual void       ViscosityMOLHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input)=0;
		virtual void       ViscosityMOLHOAdjoint(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input)=0;
		virtual void       ViscosityL1L2(IssmDouble* pviscosity,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* surf)=0;
		virtual void       ViscositySSA(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input)=0;
		virtual void       ViscosityBFS(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input,IssmDouble epseff)=0;
		virtual void       ViscosityBHO(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,IssmDouble epseff)=0;
		virtual void       ViscosityBSSA(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,IssmDouble epseff)=0;

};
#endif
