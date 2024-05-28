/*!\file Matlitho.h
 * \brief: header file for matlitho object
 */

#ifndef _MATLITHO_H_
#define _MATLITHO_H_

/*Headers:*/
/*{{{*/
#include "./Material.h"
class IoModel;
/*}}}*/

class Matlitho: public Material{

	public: 
		int	      mid;
		int          numlayers;
		IssmDouble*  radius;
		IssmDouble*  viscosity;
		IssmDouble*  lame_lambda;
		IssmDouble*  lame_mu;
		IssmDouble*  burgers_viscosity;
		IssmDouble*  burgers_mu;
		IssmDouble*  ebm_alpha;
		IssmDouble*  ebm_delta;
		IssmDouble*  ebm_taul;
		IssmDouble*  ebm_tauh;
		IssmDouble*  density;
		int*         rheologymodel;
		bool*        issolid;

		Matlitho();
		Matlitho(int matlitho_id, IoModel* iomodel, bool* issolid_in, int* rheo_in);
		~Matlitho();
		void SetMid(int matlitho_mid);

		/*Object virtual functions definitions:{{{ */
		Object *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*Material virtual functions resolution: {{{*/
		Material*  copy2(Element* element){_error_("not implemented");};
		void       Configure(Elements* elements);
		void       GetViscosity(IssmDouble* pviscosity,IssmDouble eps_eff,Gauss* gauss){_error_("not supported");};
		void       GetViscosityBar(IssmDouble* pviscosity,IssmDouble eps_eff,Gauss* gauss){_error_("not supported");};
		void       GetViscosityComplement(IssmDouble* pviscosity_complement, IssmDouble* pepsilon,Gauss* gauss){_error_("not supported");};
		void       GetViscosityDComplement(IssmDouble* pviscosity_complement, IssmDouble* pepsilon,Gauss* gauss){_error_("not supported");};
		void       GetViscosityDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon,Gauss* gauss){_error_("not supported");};
		void       GetViscosity_B(IssmDouble* pviscosity,IssmDouble eps_eff,Gauss* gauss){_error_("not supported");};
		void       GetViscosity_D(IssmDouble* pviscosity,IssmDouble eps_eff,Gauss* gauss){_error_("not supported");};
		void       GetViscosity2dDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon,Gauss* gauss){_error_("not supported");};
		IssmDouble GetA(Gauss* gauss){_error_("not supported");};
		IssmDouble GetAbar(Gauss* gauss){_error_("not supported");};
		IssmDouble GetB(Gauss* gauss){_error_("not supported");};
		IssmDouble GetBbar(Gauss* gauss){_error_("not supported");};
		IssmDouble GetD(Gauss* gauss){_error_("not supported");};
		IssmDouble GetDbar(Gauss* gauss){_error_("not supported");};
		IssmDouble GetN(){_error_("not supported");};
		bool       IsDamage(){_error_("not supported");};
		bool       IsEnhanced(){_error_("not supported");};
		void       ResetHooks();

		void       ViscosityFSDerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss){_error_("not supported");};
		void       ViscosityHODerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss){_error_("not supported");};
		void       ViscositySSADerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon,Gauss* gauss){_error_("not supported");};

		void       ViscosityFS(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){_error_("not supported");};
		void       ViscosityHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){_error_("not supported");};
		void       ViscosityMOLHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input){_error_("not supported");};
		void       ViscosityMOLHOAdjoint(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input){_error_("not supported");};
		void       ViscosityL1L2(IssmDouble* pviscosity,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* surf){_error_("not supported");};
		void       ViscositySSA(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){_error_("not supported");};
		void       ViscosityBFS(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input,IssmDouble epseff){_error_("not supported");};
		void       ViscosityBHO(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,IssmDouble epseff){_error_("not supported");};
		void       ViscosityBSSA(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,IssmDouble epseff){_error_("not supported");};

		/*}}}*/

};

#endif  /* _MATLITHO_H_ */
