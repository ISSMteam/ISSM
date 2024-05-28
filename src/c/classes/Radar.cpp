/*!\file Radar.cpp
 * \brief: Radar Object
 */

/*Headers:*/
#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./classes.h"
#include "./Elements/Element.h"
#include "./Elements/Elements.h"
#include "./ExternalResults/ExternalResult.h"
#include "./ExternalResults/Results.h"
#include "../datastructures/datastructures.h"
#include "./FemModel.h"
#include "../classes/Params/Parameters.h"
#include "../classes/gauss/Gauss.h"
#include "./Radar.h"

/*Element macros*/
#define NUMVERTICES   6
#define NUMVERTICES2D 3

/*Radar constructors, destructors :*/
Radar::Radar(){/*{{{*/
	this->definitionenum = -1;
	this->name = NULL;
}
/*}}}*/
Radar::Radar(char* in_name, int in_definitionenum){/*{{{*/
	this->definitionenum=in_definitionenum;
	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);
}
/*}}}*/
Radar::~Radar(){/*{{{*/
	if(this->name)xDelete(this->name);
}
/*}}}*/
/*Object virtual function resolution: */
Object* Radar::copy() {/*{{{*/
	Radar* out =new Radar(this->name,this->definitionenum);
	return (Object*)out;
}
/*}}}*/
void Radar::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Radar::Echo(void){/*{{{*/
	_printf_(" Radar: " << name << " " << this->definitionenum << "\n");
}
/*}}}*/
int Radar::Id(void){/*{{{*/
	return -1;
}
/*}}}*/
void Radar::Marshall(MarshallHandle* marshallhandle){/*{{{*/
	_error_("not implemented yet!");
}
/*}}}*/
int Radar::ObjectEnum(void){/*{{{*/
	return RadarEnum;
}
/*}}}*/
/*Definition virtual function resolutoin: */
int Radar::DefinitionEnum(){/*{{{*/
	return this->definitionenum;
}
/*}}}*/
char* Radar::Name(){/*{{{*/
	char* name2=xNew<char>(strlen(this->name)+1);
	xMemCpy(name2,this->name,strlen(this->name)+1);
	return name2;
}
/*}}}*/
IssmDouble Radar::Response(FemModel* femmodel){/*{{{*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		this->ComputeRadarAttenuation(element);
		this->ComputeRadarPower(element);
	}
	return 0.;
}
	/*}}}*/
void Radar::ComputeRadarAttenuation(Element* element){/*{{{*/
	//int         numvertices = element->GetNumberOfVertices();
	IssmDouble  eps0=8.85418782e-12;
	IssmDouble  eps_ice=3.17;
	IssmDouble  e=2.7183;
	IssmDouble  c=299792458;
	IssmDouble  k=1.3806488e-23;
	IssmDouble  eVtoJ=1.6e-19;
	IssmDouble  Tr_M07=252.1500;
	IssmDouble  Tr_W97=258.1500;
	IssmDouble  sig_ice_M07=9.2;
	IssmDouble  sig_ice_W97=9;
	IssmDouble  cond_H_M07=3.2;
	IssmDouble  cond_H_W97=4;
	IssmDouble  cond_Cl_M07=0.43;
	IssmDouble  cond_Cl_W97=0.55;
	IssmDouble  cond_NH_M07=0.8;
	IssmDouble  cond_NH_W97=1;
	IssmDouble  E_M07=8.1600e-20;
	IssmDouble  E_W97=9.2800e-20;
	IssmDouble  E_H_M07=3.2000e-20;
	IssmDouble  E_H_W97=3.3600e-20;
	IssmDouble  E_Cl_M07=3.0400e-20;
	IssmDouble  E_Cl_W97=3.6800e-20;
	IssmDouble  E_NH=3.6800e-20;
	IssmDouble  mol_H_hol=1.6;
	IssmDouble  mol_H_lgp=0.2;
	IssmDouble  mol_Cl_hol=0.4;
	IssmDouble  mol_Cl_lgp=1.8;
	IssmDouble  mol_NH_hol=0.5;
	IssmDouble  mol_NH_lgp=0.4;
	IssmDouble  mol_H, mol_Cl, mol_NH;
	IssmDouble  attenuation_rate_macgregor[NUMVERTICES];
	IssmDouble  attenuation_rate_wolff[NUMVERTICES];
	IssmDouble  temperature, ice_period, attenuation_rate_M07_pureice, attenuation_rate_M07_H, attenuation_rate_M07_Cl, attenuation_rate_M07_NH;
	IssmDouble  attenuation_rate_W97_pureice, attenuation_rate_W97_H, attenuation_rate_W97_Cl, attenuation_rate_W97_NH;
	IssmDouble  m1, m2, m3, m4, m5, m6, m7, m8;
	IssmDouble  w2, w3, w4, w5, w6, w7, w8;
	GaussPenta* gauss=NULL;

	/*Retrieve all inputs we will be needing: */
	Input* temp_input=element->GetInput(TemperatureEnum); _assert_(temp_input);
	Input* ice_period_input=element->GetInput(RadarIcePeriodEnum); _assert_(ice_period_input);

	/* Start looping on the number of vertices: */
	gauss=new GaussPenta();

	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/*Get ice temperature: */
		temp_input->GetInputValue(&temperature,gauss);
		ice_period_input->GetInputValue(&ice_period,gauss);

		if(ice_period>0){;
			mol_H=mol_H_hol;
			mol_Cl=mol_Cl_hol;
			mol_NH=mol_NH_hol;
			}
		else{
			mol_H=mol_H_lgp;
			mol_Cl=mol_Cl_lgp;
			mol_NH=mol_NH_lgp;
		}

		/*Compute M07 radar conductivity constant: */
		m1=(10*log10(e))/(1000*eps0*sqrt(eps_ice)*c);
		m2=E_M07/k;
		m3=cond_H_M07*mol_H;
		m4=E_H_M07/k;
		m5=cond_Cl_M07*mol_Cl;
		m6=E_Cl_M07/k;
		m7=cond_NH_M07*mol_NH;
		m8=E_NH/k;

		/*Compute MacGregor (M07) attenuation rate: */
		attenuation_rate_M07_pureice=m1*sig_ice_M07*exp(m2*((1/Tr_M07)-(1/temperature)));
		attenuation_rate_M07_H=m3*exp(m4*((1/Tr_M07)-(1/temperature)));
		attenuation_rate_M07_Cl=m5*exp(m6*((1/Tr_M07)-(1/temperature)));
		attenuation_rate_M07_NH=m7*exp(m8*((1/Tr_M07)-(1/temperature)));
		attenuation_rate_macgregor[iv]=attenuation_rate_M07_pureice+attenuation_rate_M07_H+attenuation_rate_M07_Cl+attenuation_rate_M07_NH;

	   /*Compute Wolff (W97) radar conductivity constant: */
		w2=E_W97/k;
		w3=cond_H_W97*mol_H;
		w4=E_H_W97/k;
		w5=cond_Cl_W97*mol_Cl;
		w6=E_Cl_W97/k;
		w7=cond_NH_W97*mol_NH;
		w8=E_NH/k;

		/*Compute Wolff attenuation rate: */
		attenuation_rate_W97_pureice=m1*sig_ice_W97*exp(w2*((1/Tr_W97)-(1/temperature)));
		attenuation_rate_W97_H=w3*exp(w4*((1/Tr_W97)-(1/temperature)));
		attenuation_rate_W97_Cl=w5*exp(w6*((1/Tr_W97)-(1/temperature)));
		attenuation_rate_W97_NH=w7*exp(w8*((1/Tr_W97)-(1/temperature)));
		attenuation_rate_wolff[iv]=attenuation_rate_W97_pureice+attenuation_rate_W97_H+attenuation_rate_W97_Cl+attenuation_rate_W97_NH;
	}

		/*Add Attenuation rate results into inputs*/
	   element->AddInput(RadarAttenuationMacGregorEnum,&attenuation_rate_macgregor[0],P1Enum);
		element->AddInput(RadarAttenuationWolffEnum,&attenuation_rate_wolff[0],P1Enum);

		/*Clean up*/
		delete gauss;

}/*}}}*/
void Radar::ComputeRadarPower(Element* element){/*{{{*/

	IssmDouble  *xyz_list=NULL;
	IssmDouble  power_M07[NUMVERTICES];
	IssmDouble  power_W97[NUMVERTICES];
	IssmDouble  depth[NUMVERTICES];
	IssmDouble  aircraft_elev=0.5;
	IssmDouble  eps_ice=3.15;
	IssmDouble  t_tp=273.15;         /* triple point temperature [K] */
	IssmDouble  p_tp=611.73;         /* water pressure [Pa] */
	IssmDouble  gamma=7.4200e-07; /* Clausius-Clapeyron constant [K/kPa] */
	IssmDouble  attenuation_rate_macgregor, attenuation_rate_wolff, attenuation_total_M07, attenuation_total_W97;
	IssmDouble  thickness, surface, z, temperature, geometric_loss, reflectivity;
	IssmDouble  rho_ice, gravity, pressure, pressure_melting_pt, frozen_temp, basal_temp, basal_pmp;
	GaussPenta* gauss=NULL;

	/* Get node coordinates*/
	element->GetVerticesCoordinates(&xyz_list);
	Input *atten_input_M07 = element->GetInput(RadarAttenuationMacGregorEnum); _assert_(atten_input_M07);
	Input *atten_input_W97 = element->GetInput(RadarAttenuationWolffEnum);     _assert_(atten_input_W97);
	Input *surf_input      = element->GetInput(SurfaceEnum);                   _assert_(surf_input);
	Input *thick_input     = element->GetInput(ThicknessEnum);                 _assert_(thick_input);
	Input *temp_input      = element->GetInput(TemperatureEnum);               _assert_(temp_input);

	/* Start looping on the number of vertices: */
	gauss=new GaussPenta();
	for (int iv=0;iv<NUMVERTICES;iv++){
			gauss->GaussVertex(iv);

			/*Get all the inputs: */
			atten_input_M07->GetInputValue(&attenuation_rate_macgregor,gauss);
			atten_input_W97->GetInputValue(&attenuation_rate_wolff,gauss);
			thick_input->GetInputValue(&thickness,gauss);
			temp_input->GetInputValue(&temperature,gauss);
			surf_input->GetInputValue(&surface,gauss);

			/*Compute depth below the ice surface: */
			z=xyz_list[3*iv+2];
			depth[iv]=(surface-z)/1e3;

			/*Compute total attenuation: */
			attenuation_total_M07=attenuation_rate_macgregor*depth[iv];
			attenuation_total_W97=attenuation_rate_wolff*depth[iv];

			/*Compute geometric loss: */
			geometric_loss=10*log10((depth[iv]+aircraft_elev)/sqrt(eps_ice));

			/*Compute radar power: */
			power_M07[iv]=-geometric_loss-attenuation_total_M07;
			power_W97[iv]=-geometric_loss-attenuation_total_W97;

			/*Identify basal elements: */
			if(element->IsOnBase() && iv<NUMVERTICES2D){

				/*Compute pressure melting point: */
				rho_ice=element->FindParam(MaterialsLatentheatEnum);
				gravity=element->FindParam(ConstantsGEnum);
				pressure=rho_ice*gravity*thickness;
				pressure_melting_pt=t_tp-gamma*(pressure-p_tp);

				if((temperature-pressure_melting_pt)<=-1){
					reflectivity=-40;
					}
				else if((temperature-pressure_melting_pt)>-1 && (temperature-pressure_melting_pt)<0){
					reflectivity=0;
					}
				else{
					reflectivity=70;
					}
				power_M07[iv]=power_M07[iv]+reflectivity;
				power_W97[iv]=power_W97[iv]+reflectivity;
			}
		}

	    /*Add power results into inputs*/
			element->AddInput(RadarPowerMacGregorEnum,&power_M07[0],P1Enum);
			element->AddInput(RadarPowerWolffEnum,&power_W97[0],P1Enum);

		/*Clean up and return*/
		delete gauss;
}/*}}}*/
