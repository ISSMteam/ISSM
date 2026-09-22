/*!\file SystemMatricesx
 * \brief retrieve vector from inputs in elements
 */

#include "./SystemMatricesx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../AllocateSystemMatricesx/AllocateSystemMatricesx.h"

#ifdef _HAVE_GPU_HO_
/*Add GPU header files here*/
#endif

void SystemMatricesx(Matrix<IssmDouble>** pKff, Matrix<IssmDouble>** pKfs, Vector<IssmDouble>** ppf, Vector<IssmDouble>** pdf, IssmDouble* pkmax,FemModel* femmodel, bool isAllocated, bool isGPU){

	/*intermediary: */
	int      M,N;
	int      analysisenum;
	Element *element = NULL;
	Load    *load    = NULL;

	/*output: */
	Matrix<IssmDouble> *Kff  = NULL;
	Matrix<IssmDouble> *Kfs  = NULL;
	Vector<IssmDouble> *pf   = NULL;
	Vector<IssmDouble> *df   = NULL;
	IssmDouble          kmax = 0;

	/*Display message*/
	if(VerboseModule()) _printf0_("   Allocating matrices");

	/*retrieve parameters: */
	femmodel->parameters->FindParam(&analysisenum,AnalysisTypeEnum);
	Analysis* analysis = EnumToAnalysis(analysisenum);

	/*Check if there are penalties*/
	bool ispenalty = femmodel->loads->IsPenalty();

	/*First, we might need to do a dry run to get kmax if penalties are employed*/
	if(ispenalty){

		/*Allocate Kff_temp*/
		Matrix<IssmDouble> *Kff_temp = NULL;
		AllocateSystemMatricesx(&Kff_temp,NULL,NULL,NULL,femmodel);

		/*Get complete stiffness matrix without penalties*/
		for(Object* & object : femmodel->elements->objects){
			element = xDynamicCast<Element*>(object);
			if(!element->AnyFSet() && analysisenum!=StressbalanceAnalysisEnum) continue;
			ElementMatrix* Ke = analysis->CreateKMatrix(element);
			ElementVector* pe = analysis->CreatePVector(element);
			element->ReduceMatrices(Ke,pe);
			if(Ke) Ke->AddToGlobal(Kff_temp,NULL);
			delete Ke;
			delete pe;
		}

		for(Object* & object : femmodel->loads->objects){
			load = xDynamicCast<Load*>(object);
			load->CreateKMatrix(Kff_temp,NULL);
		}
		Kff_temp->Assemble();

		/*Now, figure out maximum value of stiffness matrix, so that we can penalize it correctly: */
		kmax=Kff_temp->Norm(NORM_INF);
		delete Kff_temp;
	}

	/*Allocate stiffness matrices and load vector*/
	if(isAllocated) {
		Kff  = *pKff;
		Kfs  = *pKfs;
		pf   = *ppf;
		df   = *pdf;
	}
	else {
		AllocateSystemMatricesx(&Kff,&Kfs,&df,&pf,femmodel);
	}

	/*Display size*/
	if(VerboseModule()){
		Kff->GetSize(&M,&N);
		_printf0_(" (Kff stiffness matrix size: "<<M<<" x "<<N<<")\n");
	}

	if(VerboseModule()) _printf0_("   Assembling matrices\n");

	if(isGPU){
		#ifdef _HAVE_GPU_HO_
		/*Create Metadata structure and add as parameter*/
		GPUHOParam* gpuho_param = new GPUHOParam();
		femmodel->parameters->SetParam(gpuho_param, GPUHOParamEnum);

		/*Determine real-element metadata array/vector sizes first*/
		for(Object* & object : femmodel->elements->objects){
			element = xDynamicCast<Element*>(object);
			if(!element->IsIceInElement() || !element->AnyFSet()) continue;
			element->GetHOMetadataArraySizes();
		}

		#else
		_error_("cannot run GPU code, GPU_HO not enabled");
		#endif
	}

	/*Fill stiffness matrix and load vector from elements*/
	for(Object* & object : femmodel->elements->objects){
      element = xDynamicCast<Element*>(object);
		if(!element->AnyFSet() && analysisenum!=StressbalanceAnalysisEnum) continue;
		ElementMatrix* Ke = analysis->CreateKMatrix(element);
		ElementVector* pe = analysis->CreatePVector(element);
		element->ReduceMatrices(Ke,pe);
		if(Ke) Ke->AddToGlobal(Kff,Kfs);
		if(pe){
			pe->AddToGlobal(pf);
		}
		delete Ke;
		delete pe;
	}

	/*Fill stiffness matrix and load vector from loads*/
	for(Object* & object : femmodel->loads->objects){
      load = xDynamicCast<Load*>(object);
		load->CreateKMatrix(Kff,Kfs);
		load->CreatePVector(pf);
	}

	/*Now deal with penalties (only in loads)*/
	if(ispenalty){
		for(Object* & object : femmodel->loads->objects){
			load = xDynamicCast<Load*>(object);
			load->PenaltyCreateKMatrix(Kff,Kfs,kmax);
			load->PenaltyCreatePVector(pf,kmax);
		}
	}

	/*Create dof vector for stiffness matrix preconditioning*/
	if(pdf){
	for(Object* & object : femmodel->elements->objects){
      element = xDynamicCast<Element*>(object);
		ElementVector* de=analysis->CreateDVector(element);
		if(de) de->InsertIntoGlobal(df);
		delete de;
		}
	}

	/*Assemble matrices and vector*/
	Kff->Assemble();
	Kfs->Assemble();
	pf->Assemble();
	df->Assemble();
	//Kff->AllocationInfo();
	//Kfs->AllocationInfo();

	/*cleanu up and assign output pointers: */
	delete analysis;
	if(pKff) *pKff=Kff;
	else      delete Kff;
	if(pKfs) *pKfs=Kfs;
	else      delete Kfs;
	if(ppf)  *ppf=pf;
	else      delete pf;
	if(pdf)  *pdf=df;
	else      delete df;
	if(pkmax) *pkmax=kmax;
}
