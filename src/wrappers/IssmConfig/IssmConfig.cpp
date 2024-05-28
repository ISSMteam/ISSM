/*\file IssmConfig.c
 *\brief: get configuration names
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./IssmConfig.h"

void IssmConfigUsage(void){/*{{{*/
	_printf0_("\n");
	_printf0_("   usage: " << __FUNCT__ << "value = IssmConfig('string');\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(IssmConfig_python){

	/*input/output*/
	char       *name     = NULL;
	bool        isstring = false;
	IssmDouble  value    = 0.;
	char       *svalue   = NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments: */
	CHECKARGUMENTS(NLHS,NRHS,&IssmConfigUsage);

	/*Fetch inputs: */
	FetchData(&name,NAME);

	/*Core*/
	if(strcmp(name,"_HAVE_MPI_")==0){
		#ifdef _HAVE_MPI_
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_HAVE_PETSC_MPI_")==0){
		#ifdef _HAVE_PETSC_MPI_
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_HAVE_DAKOTA_")==0){
		#ifdef _HAVE_DAKOTA_
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_HAVE_MUMPS_")==0){
		#ifdef _HAVE_MUMPS_
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_HAVE_GSL_")==0){
		#ifdef _HAVE_GSL_
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_HAVE_TAO_")==0){
		#ifdef _HAVE_TAO_
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_HAVE_M1QN3_")==0){
		#ifdef _HAVE_M1QN3_
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_HAVE_ADOLC_")==0){
		#ifdef _HAVE_ADOLC_
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_HAVE_CODIPACK_")==0){
		#ifdef _HAVE_CODIPACK_
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_HAVE_PETSC_")==0){
		#ifdef _HAVE_PETSC_
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_PETSC_MAJOR_")==0){
		#ifdef _PETSC_MAJOR_
		value = IssmDouble(_PETSC_MAJOR_);
		#else
		_error_("_PETSC_MAJOR_ not found in config.h");
		#endif
	}
	else if(strcmp(name,"_PETSC_MINOR_")==0){
		#ifdef _PETSC_MINOR_
		value = IssmDouble(_PETSC_MINOR_);
		#else
		_error_("_PETSC_MINOR_ not found in config.h");
		#endif
	}
	else if(strcmp(name,"_PETSC_HAVE_CUDA_")==0){
		value = 0.;
		#ifdef PETSC_HAVE_CUDA
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_PETSC_HAVE_PASTIX_")==0){
		value = 0.;
		#ifdef PETSC_HAVE_PASTIX
		value = 1.;
		#endif
	}
	else if(strcmp(name,"_DAKOTA_VERSION_")==0){
		#ifdef _DAKOTA_VERSION_
		isstring = true;
		svalue =xNew<char>(strlen(_DAKOTA_VERSION_)+1);
		xMemCpy<char>(svalue,_DAKOTA_VERSION_,(strlen(_DAKOTA_VERSION_)+1));
		#else
		_error_("_DAKOTA_VERSION_ not found in config.h");
		#endif
	}
	else if(strcmp(name,"ISSM_PREFIX")==0){
		isstring = true;
		svalue =xNew<char>(strlen(ISSM_PREFIX)+1);
		xMemCpy<char>(svalue,ISSM_PREFIX,(strlen(ISSM_PREFIX)+1));
	}
	else if(strcmp(name,"PACKAGE_NAME")==0){
		isstring = true;
		svalue =xNew<char>(strlen(PACKAGE_NAME)+1);
		xMemCpy<char>(svalue,PACKAGE_NAME,(strlen(PACKAGE_NAME)+1));
	}
	else if(strcmp(name,"PACKAGE_VERSION")==0){
		isstring = true;
		svalue =xNew<char>(strlen(PACKAGE_VERSION)+1);
		xMemCpy<char>(svalue,PACKAGE_VERSION,(strlen(PACKAGE_VERSION)+1));
	}
	else if(strcmp(name,"PACKAGE_URL")==0){
		isstring = true;
		svalue =xNew<char>(strlen(PACKAGE_URL)+1);
		xMemCpy<char>(svalue,PACKAGE_URL,(strlen(PACKAGE_URL)+1));
	}
	else if(strcmp(name,"PACKAGE_BUGREPORT")==0){
		isstring = true;
		svalue =xNew<char>(strlen(PACKAGE_BUGREPORT)+1);
		xMemCpy<char>(svalue,PACKAGE_BUGREPORT,(strlen(PACKAGE_BUGREPORT)+1));
	}
	else if(strcmp(name,"PACKAGE_BUILD_DATE")==0){
		isstring = true;
		svalue =xNew<char>(strlen(PACKAGE_BUILD_DATE)+1);
		xMemCpy<char>(svalue,PACKAGE_BUILD_DATE,(strlen(PACKAGE_BUILD_DATE)+1));
	}
	else if(strcmp(name,"HOST_OS")==0){
		isstring = true;
		svalue =xNew<char>(strlen(HOST_OS)+1);
		xMemCpy<char>(svalue,HOST_OS,(strlen(HOST_OS)+1));
	}
	else if(strcmp(name,"USER_NAME")==0){
		isstring = true;
		svalue =xNew<char>(strlen(USER_NAME)+1);
		xMemCpy<char>(svalue,USER_NAME,(strlen(USER_NAME)+1));
	}
	else if(strcmp(name,"HOST_VENDOR")==0){
		isstring = true;
		svalue =xNew<char>(strlen(HOST_VENDOR)+1);
		xMemCpy<char>(svalue,HOST_VENDOR,(strlen(HOST_VENDOR)+1));
	}
	else if(strcmp(name,"HOST_OS")==0){
		isstring = true;
		svalue =xNew<char>(strlen(HOST_OS)+1);
		xMemCpy<char>(svalue,HOST_OS,(strlen(HOST_OS)+1));
	}
	else if(strcmp(name,"HOST_ARCH")==0){
		isstring = true;
		svalue =xNew<char>(strlen(HOST_ARCH)+1);
		xMemCpy<char>(svalue,HOST_ARCH,(strlen(HOST_ARCH)+1));
	}
	else{
		_error_("variable " << name << " not supported yet");
	}

	/* output: */
	if(isstring)
	 WriteData(SVALUE,svalue);
	else
	 WriteData(VALUE,value);

	/*Clean up*/
	xDelete<char>(name);
	xDelete<char>(svalue);

	/*end module: */
	MODULEEND();
}
