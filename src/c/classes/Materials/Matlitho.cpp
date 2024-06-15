/*!\file Matlitho.c
 * \brief: implementation of the Matlitho object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*Matlitho constructors and destructor*/
Matlitho::Matlitho(){/*{{{*/
	this->numlayers         = 0;
	this->radius            = NULL;
	this->viscosity         = NULL;
	this->lame_lambda       = NULL;
	this->lame_mu           = NULL;
	this->burgers_viscosity = NULL;
	this->burgers_mu        = NULL;
	this->ebm_alpha         = NULL;
	this->ebm_delta         = NULL;
	this->ebm_taul          = NULL;
	this->ebm_tauh          = NULL;
	this->density           = NULL;
	this->rheologymodel     = NULL;
	this->issolid           = NULL;
	return;
}
/*}}}*/
Matlitho::Matlitho(int matlitho_mid, IoModel* iomodel, bool* issolid_in, int* rheo_in){/*{{{*/

	this->mid=matlitho_mid;
	iomodel->FindConstant(&this->numlayers,"md.materials.numlayers");

	this->radius=xNew<IssmDouble>(this->numlayers+1);
	xMemCpy<IssmDouble>(this->radius, iomodel->Data("md.materials.radius"),this->numlayers+1);

	this->viscosity=xNew<IssmDouble>(this->numlayers);
	xMemCpy<IssmDouble>(this->viscosity, iomodel->Data("md.materials.viscosity"),this->numlayers);

	this->lame_lambda=xNew<IssmDouble>(this->numlayers);
	xMemCpy<IssmDouble>(this->lame_lambda, iomodel->Data("md.materials.lame_lambda"),this->numlayers);

	this->lame_mu=xNew<IssmDouble>(this->numlayers);
	xMemCpy<IssmDouble>(this->lame_mu, iomodel->Data("md.materials.lame_mu"),this->numlayers);

	this->burgers_viscosity=xNew<IssmDouble>(this->numlayers);
	xMemCpy<IssmDouble>(this->burgers_viscosity, iomodel->Data("md.materials.burgers_viscosity"),this->numlayers);

	this->burgers_mu=xNew<IssmDouble>(this->numlayers);
	xMemCpy<IssmDouble>(this->burgers_mu, iomodel->Data("md.materials.burgers_mu"),this->numlayers);

	this->ebm_alpha=xNew<IssmDouble>(this->numlayers);
	xMemCpy<IssmDouble>(this->ebm_alpha, iomodel->Data("md.materials.ebm_alpha"),this->numlayers);

	this->ebm_delta=xNew<IssmDouble>(this->numlayers);
	xMemCpy<IssmDouble>(this->ebm_delta, iomodel->Data("md.materials.ebm_delta"),this->numlayers);

	this->ebm_taul=xNew<IssmDouble>(this->numlayers);
	xMemCpy<IssmDouble>(this->ebm_taul, iomodel->Data("md.materials.ebm_taul"),this->numlayers);

	this->ebm_tauh=xNew<IssmDouble>(this->numlayers);
	xMemCpy<IssmDouble>(this->ebm_tauh, iomodel->Data("md.materials.ebm_tauh"),this->numlayers);

	this->density=xNew<IssmDouble>(this->numlayers);
	xMemCpy<IssmDouble>(this->density, iomodel->Data("md.materials.density"),this->numlayers);

	this->rheologymodel=xNew<int>(this->numlayers);
	xMemCpy<int>(this->rheologymodel, rheo_in, this->numlayers);

	this->issolid=xNew<bool>(this->numlayers);
	xMemCpy<bool>(this->issolid, issolid_in, this->numlayers);
}
/*}}}*/
Matlitho::~Matlitho(){/*{{{*/

	xDelete<IssmDouble>(radius);
	xDelete<IssmDouble>(viscosity);
	xDelete<IssmDouble>(lame_lambda);
	xDelete<IssmDouble>(lame_mu);
	xDelete<IssmDouble>(burgers_viscosity);
	xDelete<IssmDouble>(burgers_mu);
	xDelete<IssmDouble>(ebm_alpha);
	xDelete<IssmDouble>(ebm_delta);
	xDelete<IssmDouble>(ebm_taul);
	xDelete<IssmDouble>(ebm_tauh);
	xDelete<IssmDouble>(density);
	xDelete<int>(rheologymodel);
	xDelete<bool>(issolid);

	return;
}
/*}}}*/
void Matlitho::SetMid(int matlitho_mid){/*{{{*/
	this->mid=matlitho_mid;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Matlitho::copy() {/*{{{*/

	/*Output*/
	Matlitho* matlitho;

	/*Initialize output*/
	matlitho=new Matlitho(*this);

	/*copy fields: */
	matlitho->mid=this->mid;
	matlitho->numlayers=this->numlayers;
	if(matlitho->numlayers){
		matlitho->radius=xNew<IssmDouble>(this->numlayers+1); xMemCpy<IssmDouble>(matlitho->radius, this->radius,this->numlayers+1);
		matlitho->viscosity=xNew<IssmDouble>(this->numlayers); xMemCpy<IssmDouble>(matlitho->viscosity, this->viscosity,this->numlayers);
		matlitho->lame_lambda=xNew<IssmDouble>(this->numlayers); xMemCpy<IssmDouble>(matlitho->lame_lambda, this->lame_lambda,this->numlayers);
		matlitho->lame_mu=xNew<IssmDouble>(this->numlayers); xMemCpy<IssmDouble>(matlitho->lame_mu, this->lame_mu,this->numlayers);
		matlitho->burgers_viscosity=xNew<IssmDouble>(this->numlayers); xMemCpy<IssmDouble>(matlitho->burgers_viscosity, this->burgers_viscosity,this->numlayers);
		matlitho->burgers_mu=xNew<IssmDouble>(this->numlayers); xMemCpy<IssmDouble>(matlitho->burgers_mu, this->burgers_mu,this->numlayers);
		matlitho->ebm_alpha=xNew<IssmDouble>(this->numlayers); xMemCpy<IssmDouble>(matlitho->ebm_alpha, this->ebm_alpha,this->numlayers);
		matlitho->ebm_delta=xNew<IssmDouble>(this->numlayers); xMemCpy<IssmDouble>(matlitho->ebm_delta, this->ebm_delta,this->numlayers);
		matlitho->ebm_taul=xNew<IssmDouble>(this->numlayers); xMemCpy<IssmDouble>(matlitho->ebm_taul, this->ebm_taul,this->numlayers);
		matlitho->ebm_tauh=xNew<IssmDouble>(this->numlayers); xMemCpy<IssmDouble>(matlitho->ebm_tauh, this->ebm_tauh,this->numlayers);
		matlitho->density=xNew<IssmDouble>(this->numlayers); xMemCpy<IssmDouble>(matlitho->density, this->density,this->numlayers);
		matlitho->rheologymodel=xNew<int>(this->numlayers); xMemCpy<int>(matlitho->rheologymodel, this->rheologymodel,this->numlayers);
		matlitho->issolid=xNew<bool>(this->numlayers); xMemCpy<bool>(matlitho->issolid, this->issolid,this->numlayers);
	}

	return matlitho;
}
/*}}}*/
void Matlitho::DeepEcho(void){/*{{{*/

	this->Echo();
}		
/*}}}*/
void Matlitho::Echo(void){/*{{{*/

	_printf_("Matlitho:\n");
	_printf_("   mid: " << mid << "\n");
	_printf_("   numlayers: " << numlayers << "\n");
	_printf_("layer radius viscosity lame_lambda lame_mu burgers_viscosity burgers_mu ebm_alpha ebm_delta ebm_taul ebm_tauh density rheologymodel issolid\n");
	for (int i=0;i<numlayers;i++){
		_printf_(i << " " << radius[i] << " " << viscosity[i] << " " << lame_lambda[i] << " " << lame_mu[i] << " " << burgers_viscosity[i] << " " << burgers_mu[i] << " " << ebm_alpha[i] << " " << ebm_delta[i] << " " << ebm_taul[i] << " " << ebm_tauh[i] << " " << density[i] << " " << rheologymodel[i] << " " << issolid[i]);
	}
	return;
}
/*}}}*/
int  Matlitho::Id(void){ return mid; }/*{{{*/
/*}}}*/
void Matlitho::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = MatlithoEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->numlayers);
	if(numlayers) { 
		marshallhandle->call(this->radius,numlayers+1);
		marshallhandle->call(this->viscosity,numlayers);
		marshallhandle->call(this->lame_lambda,numlayers);
		marshallhandle->call(this->lame_mu,numlayers);
		marshallhandle->call(this->burgers_viscosity,numlayers);
		marshallhandle->call(this->burgers_mu,numlayers);
		marshallhandle->call(this->ebm_alpha,numlayers);
		marshallhandle->call(this->ebm_delta,numlayers);
		marshallhandle->call(this->ebm_taul,numlayers);
		marshallhandle->call(this->ebm_tauh,numlayers);
		marshallhandle->call(this->density,numlayers);
		marshallhandle->call(this->rheologymodel,numlayers);
		marshallhandle->call(this->issolid,numlayers);
	}
	else{
		radius            = NULL;
		viscosity         = NULL;
		lame_lambda       = NULL;
		lame_mu           = NULL;
		burgers_viscosity = NULL;
		burgers_mu        = NULL;
		ebm_alpha         = NULL;
		ebm_delta         = NULL;
		ebm_taul          = NULL;
		ebm_tauh          = NULL;
		density           = NULL;
		rheologymodel     = NULL;
		issolid           = NULL;
	}

}
/*}}}*/
int  Matlitho::ObjectEnum(void){/*{{{*/

	return MatlithoEnum;

}
/*}}}*/

/*Matlitho management: */
void       Matlitho::Configure(Elements* elementsin){/*{{{*/
	/*don't do anything, we don't have a hook to an element! As there is only 
	 * one Matlitho object!*/
}
/*}}}*/
void       Matlitho::ResetHooks(){/*{{{*/
	/*don't do anything, we don't have a hook to an element! As there is only 
	 * one Matlitho object!*/
	return;
}
/*}}}*/
