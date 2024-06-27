/*!\file Friction.h
 * \brief: header file for friction object
 */

#ifndef _FRICTION_H_
#define _FRICTION_H_

/*Headers:*/
class Inputs;
class Elements;
class Parameters;
class IoModel;
class GaussPenta;
class GaussTria;

class Friction{

	public:
		Element    *element;
		int         law;
		int         domaintype;
		int         linearize;
		IssmPDouble apply_dim;
		Input      *vx_input;
		Input      *vy_input;
		Input      *vz_input;
		IssmDouble *alpha2_list;
		IssmDouble *alpha2_complement_list;

		/*methods: */
		Friction();
		Friction(Element* element_in);
		Friction(Element* element_in, int dim);
		Friction(Element* element_in, IssmPDouble dim);
		~Friction();

		void  Echo(void);
		void  GetAlphaComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlphaHydroComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlphaTempComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlphaBuddComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlphaSchoofComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlphaCoulomb2Complement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlphaRegCoulombComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlphaWeertmanComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlpha2(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Coulomb(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Coulomb2(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Hydro(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Josh(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Shakti(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Temp(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Budd(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2WaterLayer(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Weertman(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2WeertmanTemp(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2PISM(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Schoof(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2RegCoulomb(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2RegCoulomb2(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Tsai(IssmDouble* palpha2,Gauss* gauss);

		IssmDouble EffectivePressure(Gauss* gauss);
		IssmDouble IcePressure(Gauss* gauss);
		IssmDouble SubglacialWaterPressure(Gauss* gauss);
		IssmDouble VelMag(Gauss* gauss);
		void GetBasalSlidingSpeeds(IssmDouble* pvx, Gauss* gauss);
		void GetBasalSlidingSpeeds(IssmDouble* pvx, IssmDouble* pvy, Gauss* gauss);
		void GetBasalSlidingSpeeds(IssmDouble* pvx, IssmDouble* pvy, IssmDouble* pvz, Gauss* gauss);
};

/*Friction related IO*/
void FrictionUpdateParameters(Parameters* parameters,IoModel* iomodel);
void FrictionUpdateInputs(Elements* elements,Inputs* inputs,IoModel* iomodel);

#endif  /* _FRICTION_H_ */
