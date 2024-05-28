/*!\file:  GiaDeflectionCoreArgs.h
 * \brief place holder for arguments to the GiaDeflectionCoreArgs routine
 */ 

#ifndef _GIADEFLECTIONCOREARGS_H_
#define _GIADEFLECTIONCOREARGS_H_

struct GiaDeflectionCoreArgs{

	/*inputs: */
	IssmDouble currenttime; 
	IssmDouble* hes; //loading history (in ice thickness)
	int numtimes; //loading history length
	IssmDouble ri; //radial distance from center of disk to vertex  i
	IssmDouble re; //radius of disk
	IssmDouble* times; //loading history times

	/*gia material parameters: */
	IssmDouble lithosphere_density;
	IssmDouble lithosphere_shear_modulus;
	IssmDouble lithosphere_thickness;
	IssmDouble mantle_density;
	IssmDouble mantle_shear_modulus;
	IssmDouble mantle_viscosity;

	/*gia solution parameters: */
	int iedge;

	/*ice properties: */
	IssmDouble rho_ice;

	/*constants: */
	IssmDouble yts;

	/*debug info: */
	int        idisk; //id of the element we are running the gia code in.

};

#endif
