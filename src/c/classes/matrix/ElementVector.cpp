/*!\file ElementVector.cpp
 * \brief: implementation of the ElementVector object, used to plug values from element into global load
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>
#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*ElementVector constructors and destructor*/
ElementVector::ElementVector(){/*{{{*/

	this->nrows=0;
	this->values=NULL;
	this->fglobaldoflist=NULL;

}
/*}}}*/
ElementVector::ElementVector(ElementVector* pe1, ElementVector* pe2){/*{{{*/

	/*If one of the two matrix is NULL, we copy the other one*/
	if(!pe1 && !pe2){
		_error_("Two input element matrices are NULL");
	}
	else if(!pe1){
		this->Init(pe2);
		return;
	}
	else if(!pe2){
		this->Init(pe1);
		return;
	}

	/*Initialize itransformation matrix pe[P[i]] = pe2[i]*/
	int* P=xNew<int>(pe2->nrows);

	/*1: Get the new numbering of pe2 and get size of the new matrix*/
	int gsize=pe1->nrows;
	for(int i=0;i<pe2->nrows;i++){
		bool found=false;
		for(int j=0;j<pe1->nrows;j++){
			if(pe2->gglobaldoflist[i]==pe1->gglobaldoflist[j]){
				found=true; P[i]=j; break;
			}
		}
		if(!found){
			P[i]=gsize; gsize++;
		}
	}

	/*2: Initialize static fields*/
	this->nrows=gsize;

	/*Gset and values*/
	this->gglobaldoflist=xNew<int>(this->nrows);
	this->fglobaldoflist=xNew<int>(this->nrows);
	this->values=xNewZeroInit<IssmDouble>(this->nrows);
	for(int i=0;i<pe1->nrows;i++){
		this->values[i] += pe1->values[i];
		this->gglobaldoflist[i]=pe1->gglobaldoflist[i];
		this->fglobaldoflist[i]=pe1->fglobaldoflist[i];
	}
	for(int i=0;i<pe2->nrows;i++){
		this->values[P[i]] += pe2->values[i];
		this->gglobaldoflist[P[i]]=pe2->gglobaldoflist[i];
		this->fglobaldoflist[P[i]]=pe2->fglobaldoflist[i];
	}

	/*clean-up*/
	xDelete<int>(P);
}
/*}}}*/
ElementVector::ElementVector(ElementVector* pe1, ElementVector* pe2,ElementVector* pe3){/*{{{*/

	/*Concatenate all matrices*/
	ElementVector* pe12 =new ElementVector(pe1,pe2);
	ElementVector* pe123=new ElementVector(pe12,pe3);

	/*Initialize current object with this matrix*/
	this->Init(pe123);

	/*clean-up*/
	delete pe12;
	delete pe123;
}
/*}}}*/
ElementVector::ElementVector(Node** nodes,int numnodes,Parameters* parameters,int approximation){/*{{{*/

	/*get Vector size and properties*/
	this->nrows=GetNumberOfDofs(nodes,numnodes,GsetEnum,approximation);

	/*fill values with 0: */
	this->values=xNewZeroInit<IssmDouble>(this->nrows);

	/*dof list*/
	this->gglobaldoflist=GetGlobalDofList(nodes,numnodes,GsetEnum,approximation);
	this->fglobaldoflist=GetGlobalDofList(nodes,numnodes,FsetEnum,approximation);
}
/*}}}*/
ElementVector::~ElementVector(){/*{{{*/
	xDelete<IssmDouble>(this->values);
	xDelete<int>(this->gglobaldoflist);
	xDelete<int>(this->fglobaldoflist);
}
/*}}}*/

/*ElementVector specific routines: */
void ElementVector::AddToGlobal(Vector<IssmDouble>* pf){/*{{{*/

	/*In debugging mode, check consistency (no NaN, and values not too big)*/
	this->CheckConsistency();

	/*do we have any component in the F set?*/
	bool isfset= false;
	for(int i=0;i<this->nrows;i++){
		if(this->fglobaldoflist[i]>=0){
			isfset = true;
			break;
		}
	}

	if(isfset){
		pf->SetValues(this->nrows,this->fglobaldoflist,this->values,ADD_VAL);
	}

}
/*}}}*/
void ElementVector::CheckConsistency(void){/*{{{*/
	/*Check element matrix values, only in debugging mode*/
#ifdef _ISSM_DEBUG_ 
	for (int i=0;i<this->nrows;i++){
		if (xIsNan<IssmDouble>(this->values[i])) _error_("NaN found in Element Vector");
		if (xIsInf<IssmDouble>(this->values[i])) _error_("Inf found in Element Vector");
		if (fabs( this->values[i])>1.e+50) _error_("Element Vector values exceeds 1.e+50");
	}
#endif
}
/*}}}*/
void ElementVector::Echo(void){/*{{{*/

	int i;

	_printf_("Element Vector echo:\n");
	_printf_("   nrows: " << nrows << "\n");
	_printf_("   values:\n");
	for(i=0;i<nrows;i++) _printf_(setw(4) << right << i << ": " << setw(10) << values[i] << "\n");

	_printf_("   gglobaldoflist (" << gglobaldoflist << "): ");
	if(gglobaldoflist) for(i=0;i<nrows;i++) _printf_(" " << gglobaldoflist[i] );
	_printf_(" \n");

	_printf_("   fglobaldoflist (" << fglobaldoflist << "): ");
	if(fglobaldoflist) for(i=0;i<nrows;i++) _printf_(" " << fglobaldoflist[i] );
	_printf_(" \n");
}
/*}}}*/
void ElementVector::Init(ElementVector* pe){/*{{{*/

	_assert_(pe);

	this->nrows =pe->nrows;

	this->values=xNew<IssmDouble>(this->nrows);
	xMemCpy<IssmDouble>(this->values,pe->values,this->nrows);

	this->gglobaldoflist=xNew<int>(this->nrows);
	xMemCpy<int>(this->gglobaldoflist,pe->gglobaldoflist,this->nrows);

	this->fglobaldoflist=xNew<int>(this->nrows);
	xMemCpy<int>(this->fglobaldoflist,pe->fglobaldoflist,this->nrows);
}
/*}}}*/
void ElementVector::InsertIntoGlobal(Vector<IssmDouble>* pf){/*{{{*/

	/*do we have any component in the F set?*/
	bool isfset= false;
	for(int i=0;i<this->nrows;i++){
		if(this->fglobaldoflist[i]>=0){
			isfset = true;
			break;
		}
	}

	if(isfset){
		/*add local values into global  vector, using the fglobaldoflist: */
		pf->SetValues(this->nrows,this->fglobaldoflist,this->values,INS_VAL);
	}

}
/*}}}*/
void ElementVector::SetValue(IssmDouble scalar){/*{{{*/

	for(int i=0;i<this->nrows;i++)this->values[i]=scalar;

}
/*}}}*/
void ElementVector::StaticCondensation(ElementMatrix* Ke,int bsize,int* bindices){/*{{{*/
	/* 
	 * | Kii  Kib | | Ui |    |Fi|
	 * | Kbi  Kbb | | Ub |  = |Fb|
	 *
	 * Kii Ui + Kib Ub = Fi
	 * Kbi Ui + Kbb Ub = Fb
	 *
	 * We want to remove Ub from the equation:
	 *
	 * Kii Ui + Kib inv(Kbb) (Fb - Kbi Ui) = Fi
	 *
	 * which gives:
	 *
	 * (Kii - Kib inv(Kbb) Kbi) Ui = Fi - Kib inv(Kbb) Fb
	 */

	/*Checks in debugging mode*/
	_assert_(bsize>0 && bsize<this->nrows && this->values && Ke); 
	_assert_(this->nrows==Ke->nrows);

	/*Intermediaries*/
	int         counter,i,j,isize;
	IssmDouble *Fb          = NULL;
	IssmDouble *Fi          = NULL;
	IssmDouble *Kib         = NULL;
	IssmDouble *Kbb         = NULL;
	IssmDouble *Kbbinv      = NULL;
	IssmDouble *Ftemp       = NULL;
	int        *iindices    = NULL;
	bool        flag;

	/*Get new sizes and indices*/
	isize    = this->nrows - bsize;
	iindices = xNew<int>(isize);
	counter  = 0;
	for(i=0;i<this->nrows;i++){
		flag = true;
		for(j=0;j<bsize;j++){
			if(i==bindices[j]){
				flag = false;
				break;
			}
		}
		if(flag){
			_assert_(counter<isize);
			iindices[counter++] = i;
		}
	}
	_assert_(counter == isize);

	/*Get submatrices*/
	Kib = xNew<IssmDouble>(isize*bsize);
	Kbb = xNew<IssmDouble>(bsize*bsize);
	Fb  = xNew<IssmDouble>(bsize);
	Fi  = xNew<IssmDouble>(isize);
	for(i=0;i<isize;i++) for(j=0;j<bsize;j++) Kib[i*bsize+j] = Ke->values[iindices[i]*Ke->nrows + bindices[j]];
	for(i=0;i<bsize;i++) for(j=0;j<bsize;j++) Kbb[i*bsize+j] = Ke->values[bindices[i]*Ke->nrows + bindices[j]];
	for(i=0;i<bsize;i++) Fb[i] = this->values[bindices[i]];
	for(i=0;i<isize;i++) Fi[i] = this->values[iindices[i]];

	/*Invert Kbb*/
	Kbbinv = xNew<IssmDouble>(bsize*bsize);
	switch(bsize){
		case 1:
			Kbbinv[0] = 1./Kbb[0];
			break;
		case 2:
			Matrix2x2Invert(Kbbinv,Kbb);
			break;
		case 3:
			Matrix3x3Invert(Kbbinv,Kbb);
			break;
		default:
			MatrixInverse(Kbbinv,bsize,bsize,NULL,0,NULL);
			break;
	}

	/*Calculate  Kib inv(Kbb) Fb*/
	Ftemp = xNew<IssmDouble>(isize);
	TripleMultiply(Kib,isize,bsize,0, Kbbinv,bsize,bsize,0, Fb,bsize,1,0, Ftemp,0);

	/*New Pe*/
	for(i=0;i<isize;i++) Ftemp[i] = Fi[i] - Ftemp[i];

	/*Update matrix values*/
	for(i=0;i<this->nrows;i++) this->values[i]=0.;
	for(i=0;i<isize;i++){
		this->values[iindices[i]] = Ftemp[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(Kib);
	xDelete<IssmDouble>(Kbb);
	xDelete<IssmDouble>(Kbbinv);
	xDelete<IssmDouble>(Fb);
	xDelete<IssmDouble>(Fi);
	xDelete<IssmDouble>(Ftemp);
	xDelete<int>(iindices);
	return;
}
/*}}}*/
