/*!\file Matestar.h
 * \brief: header file for matice object
 */

#ifndef MATESTAR_H_
#define MATESTAR_H_

/*Headers:*/
/*{{{*/
#include "./Material.h"
#include "../Hook.h"
class IoModel;
class Elements;
class Element;
class Loads;
class Nodes;
class Vertices;
class Materials;
class Parameters;
class Gauss;
class Input;
/*}}}*/

class Matestar: public Material{

	private: 
		int      mid;
		Hook    *helement;
		Element *element;

	public:
		/*Matestar constructors, destructors: {{{*/
		Matestar();
		Matestar(int mid,int i, IoModel* iomodel);
		~Matestar();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void  Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Material virtual functions resolution: {{{*/
		void   Configure(Elements* elements);
		Material*  copy2(Element* element);
		void       GetViscosity(IssmDouble* pviscosity, IssmDouble eps_eff,Gauss* gauss);
		void       GetViscosityBar(IssmDouble* pviscosity, IssmDouble eps_eff,Gauss* gauss);
		void       GetViscosityComplement(IssmDouble* pviscosity_complement, IssmDouble* pepsilon,Gauss* gauss);
		void       GetViscosityDComplement(IssmDouble*, IssmDouble*,Gauss* gauss);
		void       GetViscosityDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon,Gauss* gauss);
		void       GetViscosity_B(IssmDouble* pviscosity, IssmDouble eps_eff,Gauss* gauss);
		void       GetViscosity_D(IssmDouble* pviscosity, IssmDouble eps_eff,Gauss* gauss);
		void       GetViscosity2dDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon,Gauss* gauss);
		IssmDouble GetA(Gauss* gauss);
		IssmDouble GetAbar(Gauss* gauss);
		IssmDouble GetB(Gauss* gauss);
		IssmDouble GetBbar(Gauss* gauss);
		IssmDouble GetD(Gauss* gauss);
		IssmDouble GetDbar(Gauss* gauss);
		IssmDouble GetEc(Gauss* gauss);
		IssmDouble GetEcbar(Gauss* gauss);
		IssmDouble GetEs(Gauss* gauss);
		IssmDouble GetEsbar(Gauss* gauss);
		IssmDouble GetN();
		bool       IsDamage();
		bool       IsEnhanced(){_error_("not supported");};
		void       ResetHooks();
		void       SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);

		void       ViscosityFSDerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss);
		void       ViscosityHODerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss);
		void       ViscositySSADerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss);

		void       ViscosityFS(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input);
		void       ViscosityHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void       ViscosityMOLHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input);
		void       ViscosityMOLHOAdjoint(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input);
		void       ViscosityL1L2(IssmDouble* pviscosity,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* surf);
		void       ViscositySSA(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void       ViscosityBFS(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input,IssmDouble eps_eff);
		void       ViscosityBHO(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,IssmDouble eps_eff);
		void       ViscosityBSSA(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,IssmDouble eps_eff);
		/*}}}*/
		IssmDouble GetViscosityGeneral(IssmDouble vx,IssmDouble vy,IssmDouble vz,IssmDouble* dvx,IssmDouble* dvy,IssmDouble* dvz,IssmDouble eps_eff,bool isdepthaveraged,Gauss* gauss);
		IssmDouble GetViscosity_BGeneral(IssmDouble vx,IssmDouble vy,IssmDouble vz,IssmDouble* dvx,IssmDouble* dvy,IssmDouble* dvz,IssmDouble eps_eff,bool isdepthaveraged,Gauss* gauss);
};

#endif  /* _MATESTAR_H_ */
