/*!\file ElementMatrix.cpp
 * \brief: implementation of the ElementMatrix object, used to plug values from element into global stiffness matrix
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

/*ElementMatrix constructors and destructor*/
ElementMatrix::ElementMatrix(){/*{{{*/

	this->nrows=0;
	this->values=NULL;
	this->fglobaldoflist=NULL;
	this->sglobaldoflist=NULL;
}
/*}}}*/
ElementMatrix::ElementMatrix(ElementMatrix* Ke){/*{{{*/

	if(!Ke) _error_("Input Element Matrix is a NULL pointer");
	this->Init(Ke);
	return;
}
/*}}}*/
ElementMatrix::ElementMatrix(ElementMatrix* Ke1, ElementMatrix* Ke2){/*{{{*/

	/*intermediaries*/
	int i,j,counter;

	/*If one of the two matrix is NULL, we copy the other one*/
	if(!Ke1 && !Ke2){
		_error_("Two input element matrices are NULL");
	}
	else if(!Ke1){
		this->Init(Ke2);
		return;
	}
	else if(!Ke2){
		this->Init(Ke1);
		return;
	}

	/*General Case: Ke1 and Ke2 are not empty*/

	/*Initialize itransformation matrix Ke[P[i]] = Ke2[i]*/
	int* P=xNew<int>(Ke2->nrows);

	/*1: Get the new numbering of Ke2 and get size of the new matrix*/
	int gsize=Ke1->nrows;
	for(i=0;i<Ke2->nrows;i++){
		bool found=false;
		for(j=0;j<Ke1->nrows;j++){
			if(Ke2->gglobaldoflist[i]==Ke1->gglobaldoflist[j]){
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
	this->sglobaldoflist=xNew<int>(this->nrows);
	this->values=xNewZeroInit<IssmDouble>(this->nrows*this->nrows);
	for(i=0;i<Ke1->nrows;i++){
		for(j=0;j<Ke1->nrows;j++){
			this->values[i*this->nrows+j] += Ke1->values[i*Ke1->nrows+j];
		}
		this->gglobaldoflist[i]=Ke1->gglobaldoflist[i];
		this->fglobaldoflist[i]=Ke1->fglobaldoflist[i];
		this->sglobaldoflist[i]=Ke1->sglobaldoflist[i];
	}
	for(i=0;i<Ke2->nrows;i++){
		for(j=0;j<Ke2->nrows;j++){
			this->values[P[i]*this->nrows+P[j]] += Ke2->values[i*Ke2->nrows+j];
		}
		this->gglobaldoflist[P[i]]=Ke2->gglobaldoflist[i];
		this->fglobaldoflist[P[i]]=Ke2->fglobaldoflist[i];
		this->sglobaldoflist[P[i]]=Ke2->sglobaldoflist[i];
	}

	/*clean-up*/
	xDelete<int>(P);
}
/*}}}*/
ElementMatrix::ElementMatrix(ElementMatrix* Ke1, ElementMatrix* Ke2,ElementMatrix* Ke3){/*{{{*/

	/*Concatenate all matrices*/
	ElementMatrix* Ke12 =new ElementMatrix(Ke1,Ke2);
	ElementMatrix* Ke123=new ElementMatrix(Ke12,Ke3);

	/*Initialize current object with this matrix*/
	this->Init(Ke123);

	/*clean-up*/
	delete Ke12;
	delete Ke123;
}
/*}}}*/
ElementMatrix::ElementMatrix(Node** nodes,int numnodes,Parameters* parameters,int approximation){/*{{{*/

	/*get Matrix size and properties*/
	this->nrows=GetNumberOfDofs(nodes,numnodes,GsetEnum,approximation);

	/*fill values with 0: */
	this->values=xNewZeroInit<IssmDouble>(this->nrows*this->nrows);

	this->gglobaldoflist=GetGlobalDofList(nodes,numnodes,GsetEnum,approximation);
	this->fglobaldoflist=GetGlobalDofList(nodes,numnodes,FsetEnum,approximation);
	this->sglobaldoflist=GetGlobalDofList(nodes,numnodes,SsetEnum,approximation);
}
/*}}}*/
ElementMatrix::~ElementMatrix(){/*{{{*/

	xDelete<IssmDouble>(this->values);
	xDelete<int>(this->gglobaldoflist);
	xDelete<int>(this->fglobaldoflist);
	xDelete<int>(this->sglobaldoflist);
}
/*}}}*/

/*ElementMatrix specific routines: */
void ElementMatrix::AddDiagonalToGlobal(Vector<IssmDouble>* pf){/*{{{*/

	IssmDouble* localvalues=NULL;

	/*Check that pf is not NULL*/
	_assert_(pf); 

	/*In debugging mode, check consistency (no NaN, and values not too big)*/
	this->CheckConsistency();

	/*do we have any component in the F set?*/
	int fsize = 0;
	for(int i=0;i<this->nrows;i++){
		if(this->fglobaldoflist[i]>=0) fsize++;
	}

	if(fsize){
		/*first, retrieve values that are in the f-set from the g-set values matrix: */
		localvalues=xNew<IssmDouble>(this->nrows);
		for(int i=0;i<this->nrows;i++){
			localvalues[i] = this->values[this->nrows*i + i];
		}

		/*add local values into global  matrix, using the fglobaldoflist: */
		pf->SetValues(fsize,this->fglobaldoflist,localvalues,ADD_VAL);

		/*Free resources:*/
		xDelete<IssmDouble>(localvalues);
	}

}
/*}}}*/
void ElementMatrix::AddToGlobal(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/

	/*Check that Kff has been alocated in debugging mode*/
	_assert_(Kff);

	/*If Kfs is not provided, call the other function*/
	if(!Kfs){
		this->AddToGlobal(Kff);
		return;
	}

	/*In debugging mode, check consistency (no NaN, and values not too big)*/
	this->CheckConsistency();

	/*do we have any component in the F or S set?*/
   bool is_fset= false;
	bool is_sset= false;
   for(int i=0;i<this->nrows;i++){
      if(this->fglobaldoflist[i]>=0){
         is_fset = true;
      }
		else{
			_assert_(this->sglobaldoflist[i]>=0);
			is_sset = true;
		}
   }

	/*only use row dofs to add values into global matrices: */
	if(is_fset){
		Kff->SetValues(this->nrows,this->fglobaldoflist,this->nrows,this->fglobaldoflist,this->values,ADD_VAL);
	}
	if(is_fset && is_sset){
		Kfs->SetValues(this->nrows,this->fglobaldoflist,this->nrows,this->sglobaldoflist,this->values,ADD_VAL);
	}
}
/*}}}*/
void ElementMatrix::AddToGlobal(Matrix<IssmDouble>* Jff){/*{{{*/

	/*Check that Jff is not NULL*/
	_assert_(Jff); 

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
		Jff->SetValues(this->nrows,this->fglobaldoflist,this->nrows,this->fglobaldoflist,this->values,ADD_VAL);
	}
}
/*}}}*/
void ElementMatrix::CheckConsistency(void){/*{{{*/
	/*Check element matrix values, only in debugging mode*/
	#ifdef _ISSM_DEBUG_ 
	for (int i=0;i<this->nrows;i++){
		for(int j=0;j<this->nrows;j++){
			if (xIsNan<IssmDouble>(this->values[i*this->nrows+j])) _error_("NaN found in Element Matrix");
			if (xIsInf<IssmDouble>(this->values[i*this->nrows+j])) _error_("Inf found in Element Matrix");
			if (fabs(this->values[i*this->nrows+j])>1.e+50) _error_("Element Matrix values exceeds 1.e+50");
		}
	}
	#endif
}
/*}}}*/
void ElementMatrix::Echo(void){/*{{{*/

	int i,j;
	_printf_("Element Matrix echo:\n");
	_printf_("   nrows: " << this->nrows << "\n");

	_printf_("   values:\n");
	for(i=0;i<nrows;i++){
		_printf_(setw(4) << right << i << ": ");
		for(j=0;j<nrows;j++) _printf_( " " << setw(11) << setprecision (5) << right << values[i*nrows+j]);
		_printf_("\n");
	}

	_printf_("   gglobaldoflist (" << gglobaldoflist << "): ");
	if(gglobaldoflist) for(i=0;i<nrows;i++) _printf_(" " << gglobaldoflist[i]); _printf_("\n");

	_printf_("   fglobaldoflist  (" << fglobaldoflist << "): ");
	for(i=0;i<nrows;i++)_printf_(" " << fglobaldoflist[i]); _printf_(" \n");

	_printf_("   sglobaldoflist  (" << sglobaldoflist << "): ");
	for(i=0;i<nrows;i++)_printf_(" " << sglobaldoflist[i]); _printf_(" \n");
}
/*}}}*/
bool ElementMatrix::HasDof(int dof,int set){/*{{{*/

	if(set==FsetEnum){
		for(int i=0;i<this->nrows;i++) if(this->fglobaldoflist[i] == dof) return true;
	}
	else if(set==GsetEnum){
		for(int i=0;i<this->nrows;i++) if(this->gglobaldoflist[i] == dof) return true;
	}
	else if(set==SsetEnum){
		for(int i=0;i<this->nrows;i++) if(this->sglobaldoflist[i] == dof) return true;
	}
	else{
		_error_("not supported yet");
	}

	return false;
}
/*}}}*/
void ElementMatrix::Init(ElementMatrix* Ke){/*{{{*/

	_assert_(Ke);
	_assert_(this);

	this->nrows =Ke->nrows;

	this->values=xNew<IssmDouble>(this->nrows*this->nrows);
	xMemCpy<IssmDouble>(this->values,Ke->values,this->nrows*this->nrows);

	this->gglobaldoflist=xNew<int>(this->nrows);
	xMemCpy<int>(this->gglobaldoflist,Ke->gglobaldoflist,this->nrows);
	this->fglobaldoflist=xNew<int>(this->nrows);
	xMemCpy<int>(this->fglobaldoflist,Ke->fglobaldoflist,this->nrows);
	this->sglobaldoflist=xNew<int>(this->nrows);
	xMemCpy<int>(this->sglobaldoflist,Ke->sglobaldoflist,this->nrows);
}
/*}}}*/
void ElementMatrix::Lump(void){/*{{{*/

	for(int i=0;i<this->nrows;i++){
		for(int j=0;j<this->nrows;j++){
			if(i!=j){
				this->values[i*this->nrows+i] += this->values[i*this->nrows+j];
				this->values[i*this->nrows+j]  = 0.;
			}
		}
	}

	return;
}
/*}}}*/
void ElementMatrix::Transpose(void){/*{{{*/

	/*Intermediaries*/
	ElementMatrix* Ke_copy=new ElementMatrix(this);

	/*Transpose values*/
	for (int i=0;i<this->nrows;i++) for(int j=0;j<this->nrows;j++) this->values[i*this->nrows+j]=Ke_copy->values[j*Ke_copy->nrows+i];

	/*Clean up and return*/
	delete Ke_copy;
	return;
}
/*}}}*/
void ElementMatrix::StaticCondensation(int bsize,int* bindices){/*{{{*/
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
	_assert_(bsize>0 && bsize<this->nrows && this->values); 

	/*Intermediaries*/
	int         counter,i,j,isize;
	IssmDouble *Kii         = NULL;
	IssmDouble *Kib         = NULL;
	IssmDouble *Kbi         = NULL;
	IssmDouble *Kbb         = NULL;
	IssmDouble *Kbbinv      = NULL;
	IssmDouble *Ktemp       = NULL;
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
	Kii = xNew<IssmDouble>(isize*isize);
	Kib = xNew<IssmDouble>(isize*bsize);
	Kbi = xNew<IssmDouble>(bsize*isize);
	Kbb = xNew<IssmDouble>(bsize*bsize);
	for(i=0;i<isize;i++) for(j=0;j<isize;j++) Kii[i*isize+j] = this->values[iindices[i]*this->nrows + iindices[j]];
	for(i=0;i<isize;i++) for(j=0;j<bsize;j++) Kib[i*bsize+j] = this->values[iindices[i]*this->nrows + bindices[j]];
	for(i=0;i<bsize;i++) for(j=0;j<isize;j++) Kbi[i*isize+j] = this->values[bindices[i]*this->nrows + iindices[j]];
	for(i=0;i<bsize;i++) for(j=0;j<bsize;j++) Kbb[i*bsize+j] = this->values[bindices[i]*this->nrows + bindices[j]];

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

	/*Calculate  Kib inv(Kbb) Kbi*/
	Ktemp = xNew<IssmDouble>(isize*isize);
	TripleMultiply(Kib,isize,bsize,0, Kbbinv,bsize,bsize,0, Kbi,bsize,isize,0, Ktemp,0);

	/*New Ke*/
	for(i=0;i<isize*isize;i++) Ktemp[i] = Kii[i] - Ktemp[i];

	/*Update matrix values*/
	for(i=0;i<this->nrows*this->nrows;i++) this->values[i]=0.;
	for(i=0;i<isize;i++){
		for(j=0;j<isize;j++){
			this->values[iindices[i]*this->nrows + iindices[j]] = Ktemp[i*isize+j];
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(Kii);
	xDelete<IssmDouble>(Kib);
	xDelete<IssmDouble>(Kbi);
	xDelete<IssmDouble>(Kbb);
	xDelete<IssmDouble>(Kbbinv);
	xDelete<IssmDouble>(Ktemp);
	xDelete<int>(iindices);
	return;
}
/*}}}*/
