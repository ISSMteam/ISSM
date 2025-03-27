/*!\file: sealevelchange_core.cpp
 * \brief: core of the sea-level change solution 
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../classes/Inputs/TriaInput.h"
#include "../classes/Inputs/TransientInput.h"
#include "../classes/Inputs/DatasetInput.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

/*support routines local definitions:{{{*/
void TransferForcing(FemModel* femmodel,int forcingenum);
void TransferSealevel(FemModel* femmodel,int forcingenum);
bool slcconvergence(IssmDouble* RSLg,IssmDouble* RSLg_old,IssmDouble eps_rel,IssmDouble eps_abs, IssmDouble totaloceanarea,FemModel* femmodel);
IssmDouble  SealevelloadsOceanAverage(GrdLoads* loads, Vector<IssmDouble>* oceanareas, Vector<IssmDouble>* subelementoceanareas, IssmDouble totaloceanarea);
void PolarMotion(IssmDouble* m, FemModel* femmodel,GrdLoads* loads, SealevelGeometry* slgeom, bool computefuture);
void SealevelchangeUpdateViscousTimeSeries(FemModel* femmodel);
void ConserveOceanMass(FemModel* femmodel,GrdLoads* loads, IssmDouble offset, SealevelGeometry* slgeom);
void ivins_deformation_core(FemModel* femmodel);
IssmDouble* CombineLoads(IssmDouble* load,IssmDouble* subload,FemModel* femmodel, SealevelGeometry* slgeom,int loadtype,int nel);
void slc_geometry_cleanup(SealevelGeometry* slgeom, FemModel* femmodel);
/*}}}*/

/*main cores:*/
void              sealevelchange_core(FemModel* femmodel){ /*{{{*/

	SealevelGeometry* slgeom=NULL;

	/*Start profiler*/
	femmodel->profiler->Start(SLRCORE);

	/*Parameters, variables:*/
	bool save_results;

	/*Retrieve parameters:*/
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);

	/*Verbose: */
	if(VerboseSolution()) _printf0_("   computing sea level change\n");

	/*set SLR configuration: */
	femmodel->SetCurrentConfiguration(SealevelchangeAnalysisEnum);

	/*Run coupler input transfer:*/
	couplerinput_core(femmodel);

	/*run geometry core: */
	slgeom=sealevelchange_geometry(femmodel);

	/*any external forcings?:*/
	solidearthexternal_core(femmodel);

	/*Run geodetic:*/
	grd_core(femmodel,slgeom);

	/*Run steric core for sure:*/
	dynstr_core(femmodel);

	/*Run coupler output transfer: */
	coupleroutput_core(femmodel);

	/*Save results: */
	if(save_results){
		int     numoutputs;
		char **requested_outputs = NULL;
		if(VerboseSolution()) _printf0_("   saving results\n");
		femmodel->parameters->FindParam(&requested_outputs,&numoutputs,SealevelchangeRequestedOutputsEnum);
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
		if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}
	}

	/*End profiler*/
	femmodel->profiler->Stop(SLRCORE);

	/*Free resources:*/
	slc_geometry_cleanup(slgeom, femmodel);
}
/*}}}*/
void              solidearthexternal_core(FemModel* femmodel){ /*{{{*/

	/*variables:*/
	Vector<IssmDouble> *bedrock  = NULL; 
	Vector<IssmDouble> *bedrock_rate = NULL;
	Vector<IssmDouble> *bedrockeast  = NULL; 
	Vector<IssmDouble> *bedrockeast_rate = NULL;
	Vector<IssmDouble> *bedrocknorth  = NULL; 
	Vector<IssmDouble> *bedrocknorth_rate = NULL;
	Vector<IssmDouble> *geoid= NULL; 
	Vector<IssmDouble> *geoid_rate= NULL; 
	int horiz=0;
	int modelid=-1;
	int  isexternal=0;

	/*parameters: */
	IssmDouble          dt;

	/*Retrieve parameters:*/
	femmodel->parameters->FindParam(&horiz,SolidearthSettingsHorizEnum);
	femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	femmodel->parameters->FindParam(&isexternal,SolidearthIsExternalEnum); 

	/*Early return:*/
	if (!isexternal)return;

	/*Verbose: */
	if(VerboseSolution()) _printf0_("	  computing external solid earth contributions\n");

	/*Retrieve geoid viscous and elastic rates, bedrock uplift viscous and elastic rates + steric rate, as vectors:*/
	GetVectorFromInputsx(&bedrock,femmodel,BedEnum,VertexSIdEnum);
	GetVectorFromInputsx(&geoid,femmodel,SealevelEnum,VertexSIdEnum); //In ISSM, Sealevel is absolute.
	if(horiz){
		GetVectorFromInputsx(&bedrockeast,femmodel,BedEastEnum,VertexSIdEnum);
		GetVectorFromInputsx(&bedrocknorth,femmodel,BedNorthEnum,VertexSIdEnum);
	}

	GetVectorFromInputsx(&geoid_rate,femmodel,SolidearthExternalGeoidRateEnum,VertexSIdEnum);
	GetVectorFromInputsx(&bedrock_rate,femmodel,SolidearthExternalDisplacementUpRateEnum,VertexSIdEnum);
	if(horiz){
		GetVectorFromInputsx(&bedrockeast_rate,femmodel,SolidearthExternalDisplacementEastRateEnum,VertexSIdEnum);
		GetVectorFromInputsx(&bedrocknorth_rate,femmodel,SolidearthExternalDisplacementNorthRateEnum,VertexSIdEnum);
	}

	/*compute: sea level change = initial sea level + (N_gia_rate+N_esa_rate)  * dt + steric_rate + dynamic_rate dt*/
	geoid->AXPY(geoid_rate,dt);
	bedrock->AXPY(bedrock_rate,dt);
	if(horiz){
		bedrockeast->AXPY(bedrockeast_rate,dt);
		bedrocknorth->AXPY(bedrocknorth_rate,dt);
	}

	/*update element inputs:*/
	InputUpdateFromVectorx(femmodel,bedrock,BedEnum,VertexSIdEnum);	
	InputUpdateFromVectorx(femmodel,geoid,SealevelEnum,VertexSIdEnum);	
	if(horiz){
		InputUpdateFromVectorx(femmodel,bedrockeast,BedEastEnum,VertexSIdEnum);	
		InputUpdateFromVectorx(femmodel,bedrocknorth,BedNorthEnum,VertexSIdEnum);	
	}

	/*Free resources:*/	
	delete bedrock; delete bedrock_rate;
	delete geoid; delete geoid_rate;
	if(horiz){
		delete bedrockeast; delete bedrockeast_rate;
		delete bedrocknorth; delete bedrocknorth_rate;
	}
}
/*}}}*/
void              couplerinput_core(FemModel* femmodel){  /*{{{*/

	/*Be very careful here, everything is well thought through, do not remove 
	 * without taking big risks:*/

	/*parameters:*/
	int  iscoupling;
	int  modelid,earthid;
	int  count,frequency;
	int  horiz;

	/*retrieve more parameters:*/
	femmodel->parameters->FindParam(&iscoupling,IsSlcCouplingEnum);
	femmodel->parameters->FindParam(&frequency,SolidearthSettingsRunFrequencyEnum);
	femmodel->parameters->FindParam(&count,SealevelchangeRunCountEnum);

	if(iscoupling){
		femmodel->parameters->FindParam(&modelid,ModelIdEnum);
		femmodel->parameters->FindParam(&earthid,EarthIdEnum);
	}
	else{
		/* we are here, we are not running in a coupler, so we will indeed compute SLR,
		 * so make sure we are identified as being the Earth.:*/
		modelid=1; earthid=1; 
	}

	/*if we are carrying loads but are not yet computing grd core, accumulate them and skip 
	 * the rest: */
	if (count!=frequency){
		if (count>frequency){
			count=1;
			femmodel->parameters->SetParam(count,SealevelchangeRunCountEnum); 
		}
		return;
	}

	/*Basins are supposed to accumulate loads and hand them over to the Earth
	  for slr computations every "frequency" time steps. If we are here, we
	  have reached the "frequency"'th time step, and we are going to pick up
	  the old loads, the current loads, and send them to the Earth for slr
	  computations.  So the Earth is never supposed to compute loads. Except,
	  when we are running the Earth as a giant basin (ex: running a
	  Peltier-style GIA model, or running a "Snow Ball" Earth) which only
	  happens when we are not coupled, at which point we are sending these
	  loads to ourselves (hence the convoluted condition).
	  */

	/*transer loads from basins to Earth for grd core computations. :*/
	if(iscoupling){

		/*transfer ice thickness change load from basins to earth: */
		TransferForcing(femmodel,DeltaIceThicknessEnum);
		TransferForcing(femmodel,DeltaBottomPressureEnum);
		TransferForcing(femmodel,DeltaTwsEnum);
		TransferForcing(femmodel,MaskOceanLevelsetEnum);
		TransferForcing(femmodel,MaskIceLevelsetEnum);

		/*transfer external forcings back to Earth:*/
		TransferSealevel(femmodel,BedEnum);
		TransferSealevel(femmodel,SealevelEnum);
		if(horiz){
			TransferSealevel(femmodel,BedEastEnum);
			TransferSealevel(femmodel,BedNorthEnum);
		}
	}

}; /*}}}*/
void              grd_core(FemModel* femmodel, SealevelGeometry* slgeom) { /*{{{*/

	/*variables:{{{*/
	int nel;
	BarystaticContributions* barycontrib=NULL;
	GenericParam<BarystaticContributions*>* barycontribparam=NULL;
	IssmDouble polarmotionvector[3]={0};

	GrdLoads*              loads=NULL;
	IssmDouble*    oldsealevelloads=NULL;
	Vector<IssmDouble>*    oceanareas=NULL;
	IssmDouble             totaloceanarea;
	Vector<IssmDouble>*    subelementoceanareas=NULL;
	IssmDouble             oceanaverage;
	bool                   scaleoceanarea=false;
	IssmDouble             rho_water;

	IssmDouble           eps_rel;
	IssmDouble           eps_abs;
	int                  max_nonlinear_iterations;
	int                  iterations=0;
	int                  step;
	IssmDouble           time; 

	int  modelid,earthid;
	int  horiz;
	int  count,frequency,iscoupling;
	int  grd=0;
	int  grdmodel; 
	int  sealevelloading=0;
	bool sal=false;
	bool viscous=false;
	bool rotation=false;
	bool planethasocean=false;
	bool computefuture=false;
	IssmDouble*           sealevelpercpu=NULL;

	/*}}}*/

	/*retrieve parameters:{{{*/
	femmodel->parameters->FindParam(&grd,SolidearthSettingsGRDEnum); 
	femmodel->parameters->FindParam(&grdmodel,GrdModelEnum);
	femmodel->parameters->FindParam(&frequency,SolidearthSettingsRunFrequencyEnum);
	femmodel->parameters->FindParam(&count,SealevelchangeRunCountEnum);
	femmodel->parameters->FindParam(&step,StepEnum);
	femmodel->parameters->FindParam(&time,TimeEnum);
	femmodel->parameters->FindParam(&sealevelloading,SolidearthSettingsSealevelLoadingEnum);
	femmodel->parameters->FindParam(&max_nonlinear_iterations,SolidearthSettingsMaxiterEnum);
	femmodel->parameters->FindParam(&viscous,SolidearthSettingsViscousEnum);
	femmodel->parameters->FindParam(&rotation,SolidearthSettingsRotationEnum);
	femmodel->parameters->FindParam(&planethasocean,SolidearthSettingsGrdOceanEnum);
	/*}}}*/

	/*only run if grd was requested, if we are the earth, and we have reached
	 * the necessary number of time steps dictated by :*/
	if(!grd)            return;
	if(count!=frequency)return;
	femmodel->parameters->FindParam(&iscoupling,IsSlcCouplingEnum);
	if(iscoupling){
		femmodel->parameters->FindParam(&modelid,ModelIdEnum);
		femmodel->parameters->FindParam(&earthid,EarthIdEnum);
		if(modelid!=earthid)return;
	}
	/*branch directly to Ivins deformation core if requested:*/
	if(grdmodel==IvinsEnum){
		ivins_deformation_core(femmodel);
		return;
	}

	/*Verbose: */
	if(VerboseSolution()) _printf0_("	  computing GRD patterns\n");

	/*retrieve parameters: {{{*/ 
	femmodel->parameters->FindParam(&scaleoceanarea,SolidearthSettingsOceanAreaScalingEnum);
	barycontribparam = xDynamicCast<GenericParam<BarystaticContributions*>*>(femmodel->parameters->FindParamObject(BarystaticContributionsEnum));
	barycontrib=barycontribparam->GetParameterValue();
	femmodel->parameters->FindParam(&rho_water,MaterialsRhoSeawaterEnum);
	femmodel->parameters->FindParam(&eps_rel,SolidearthSettingsReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,SolidearthSettingsAbstolEnum);
	femmodel->parameters->FindParam(&horiz,SolidearthSettingsHorizEnum);
	femmodel->parameters->FindParam(&sal,SolidearthSettingsSelfAttractionEnum);
	/*}}}*/

	/*initialize loads and sea level loads:*/
	femmodel->parameters->FindParam(&nel,MeshNumberofelementsEnum);

	loads=new GrdLoads(nel,slgeom);
	subelementoceanareas=new Vector<IssmDouble>(slgeom->nbar[SLGEOM_OCEAN]);
	oceanareas=new Vector<IssmDouble>(nel);
	sealevelpercpu=xNewZeroInit<IssmDouble>(femmodel->vertices->Size());

	if(VerboseSolution()) _printf0_("	  starting  GRD convolutions\n");

	/*update viscous time series to keep up with time stepping:*/
	SealevelchangeUpdateViscousTimeSeries(femmodel);

	/*buildup loads: */
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->SealevelchangeBarystaticLoads(loads, barycontrib,slgeom); 
	}

	//broadcast loads 
	loads->BroadcastLoads();

	/*skip computation of sea level equation if requested, which means sea level loads should be zeroed */
	if(!sealevelloading){
		loads->sealevelloads=xNewZeroInit<IssmDouble>(nel);
		loads->subsealevelloads=xNewZeroInit<IssmDouble>(slgeom->nbar[SLGEOM_OCEAN]);
		PolarMotion(&polarmotionvector[0],femmodel,loads, slgeom, computefuture=true);
		goto deformation;
	}

	if(VerboseSolution()) _printf0_("	  converging GRD convolutions\n");
	for(;;){

		//compute polar motion:
		PolarMotion(&polarmotionvector[0],femmodel,loads,slgeom,computefuture=false);

		oldsealevelloads=xNewZeroInit<IssmDouble>(nel);
		if (loads->sealevelloads){
			xMemCpy<IssmDouble>(oldsealevelloads,loads->sealevelloads,nel);
		}

		/*convolve load and sealevel loads on oceans:*/
		loads->Combineloads(nel,slgeom); //This combines loads and sealevelloads into a single vector 
		for(Object* & object : femmodel->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			element->SealevelchangeConvolution(sealevelpercpu, loads, polarmotionvector,slgeom);
		}

		/*retrieve sea level average  and ocean area:*/
		for(Object* & object : femmodel->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			element->SealevelchangeOceanAverage(loads, oceanareas, subelementoceanareas, sealevelpercpu, slgeom);
		}

		loads->AssembleSealevelLoads();

		/*compute ocean areas:*/
		if(!loads->sealevelloads){ //first time in the loop
			oceanareas->Assemble(); 
			subelementoceanareas->Assemble();
			oceanareas->Sum(&totaloceanarea); _assert_(totaloceanarea>0.);
			if(scaleoceanarea) totaloceanarea=3.619e+14; // use true ocean area, m^2
		}

		//Conserve ocean mass: 
		oceanaverage=SealevelloadsOceanAverage(loads, oceanareas,subelementoceanareas, totaloceanarea);
		ConserveOceanMass(femmodel,loads,barycontrib->Total()/totaloceanarea -oceanaverage,slgeom);

		//broadcast sea level loads 
		loads->BroadcastSealevelLoads();

		if (!sal) {xDelete<IssmDouble>(oldsealevelloads); break;}

		//convergence?
		if(slcconvergence(loads->sealevelloads,oldsealevelloads,eps_rel,eps_abs,totaloceanarea,femmodel)){
			xDelete<IssmDouble>(oldsealevelloads); break;
		}

		//early return?
		if(iterations>=max_nonlinear_iterations){
			xDelete<IssmDouble>(oldsealevelloads); break;
		}
		iterations++;
		xDelete<IssmDouble>(oldsealevelloads);
	}

	//recompute polar motion one final time, this time updating viscous stacks for future time steps
	if (viscous)	PolarMotion(&polarmotionvector[0],femmodel,loads, slgeom, computefuture=true);

	deformation:

	if(VerboseSolution()) _printf0_("	  deformation GRD convolutions\n");

	/*convolve loads and sea level loads to get the deformation:*/
	loads->Combineloads(nel,slgeom); //This combines loads and sealevelloads into a single vector 
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->SealevelchangeDeformationConvolution(sealevelpercpu, loads, polarmotionvector,slgeom);
	}

	if(VerboseSolution()) _printf0_("	  updating GRD fields\n");

	if(planethasocean){
		if (!sealevelloading){ //we haven't done so before, so we need to compute the ocean average and area
			if(loads->sealevelloads)xDelete<IssmDouble>(loads->sealevelloads); loads->sealevelloads=NULL; //needed to trigger the calculation of areas
			/*retrieve sea level average  and ocean area:*/
			for(Object* & object : femmodel->elements->objects){
				Element* element = xDynamicCast<Element*>(object);
				element->SealevelchangeOceanAverage(loads, oceanareas, subelementoceanareas, sealevelpercpu, slgeom);
			}
			loads->sealevelloads=xNewZeroInit<IssmDouble>(nel);
			loads->AssembleSealevelLoads();
			/*compute ocean areas:*/
			oceanareas->Assemble(); 
			subelementoceanareas->Assemble();
			oceanareas->Sum(&totaloceanarea); _assert_(totaloceanarea>0.);
			if(scaleoceanarea) totaloceanarea=3.619e+14; // use true ocean area, m^2

			//Conserve ocean mass
			//Note that here we create sea-level loads but they will not generate GRD as we have already run all the convolutions
			oceanaverage=SealevelloadsOceanAverage(loads, oceanareas,subelementoceanareas, totaloceanarea);
			ConserveOceanMass(femmodel,loads,barycontrib->Total()/totaloceanarea - oceanaverage,slgeom);
		}

		femmodel->inputs->Shift(SealevelGRDEnum,barycontrib->Total()/rho_water/totaloceanarea- oceanaverage/rho_water);

		//cumulate barystatic contributions and save to results: 
		barycontrib->Cumulate(femmodel->parameters);
		barycontrib->Save(femmodel->results,femmodel->parameters,totaloceanarea);
		barycontrib->Reset();
	}

	if (rotation) {
		femmodel->results->AddResult(new GenericExternalResult<IssmDouble>(femmodel->results->Size()+1,SealevelchangePolarMotionXEnum,polarmotionvector[0],step,time));
		femmodel->results->AddResult(new GenericExternalResult<IssmDouble>(femmodel->results->Size()+1,SealevelchangePolarMotionYEnum,polarmotionvector[1],step,time));
		femmodel->results->AddResult(new GenericExternalResult<IssmDouble>(femmodel->results->Size()+1,SealevelchangePolarMotionZEnum,polarmotionvector[2],step,time));
	}

	/*Update bedrock motion and geoid:*/
	femmodel->inputs->AXPY(1,SealevelGRDEnum,SealevelEnum);
	femmodel->inputs->AXPY(1,BedGRDEnum,BedEnum);
	if(horiz){
		femmodel->inputs->AXPY(1,BedEastGRDEnum,BedEastEnum);
		femmodel->inputs->AXPY(1,BedNorthGRDEnum, BedNorthEnum);
	}

	/*Free resources:*/
	delete loads;
	delete subelementoceanareas;
	delete oceanareas;
	xDelete<IssmDouble>(sealevelpercpu); 
}
/*}}}*/
void              dynstr_core(FemModel* femmodel){ /*{{{*/

	/*variables:*/
	Vector<IssmDouble> *sealevel  = NULL; 
	Vector<IssmDouble> *deltadsl  = NULL; 
	Vector<IssmDouble> *deltastr = NULL;

	/*parameters: */
	int  step;
	bool isocean=false;
	IssmDouble time;

	IssmDouble cumgmtslc=0;
	IssmDouble cumbslc=0;
	IssmDouble cumgmslc=0;
	IssmDouble gmtslc=0;

	/*early return if we have no ocean transport:*/
	femmodel->parameters->FindParam(&isocean,TransientIsoceantransportEnum);
	if(!isocean)return;

	/*Verbose: */
	if(VerboseSolution()) _printf0_("	  computing steric and dynamic sea level change\n");

	/*Retrieve sealevel and add steric + dynamic rates:*/
	GetVectorFromInputsx(&sealevel,femmodel,SealevelEnum,VertexSIdEnum);
	GetVectorFromInputsx(&deltadsl,femmodel,DeltaDslEnum,VertexSIdEnum);
	GetVectorFromInputsx(&deltastr,femmodel,DeltaStrEnum,VertexSIdEnum);

	/*compute: sea level change = initial sea level + steric + dynamic*/
	sealevel->AXPY(deltadsl,1);
	sealevel->AXPY(deltastr,1);

	/*cumulate thermal steric rate:*/
	femmodel->parameters->FindParam(&cumgmtslc,CumGmtslcEnum); 
	femmodel->parameters->FindParam(&cumbslc,CumBslcEnum); 

	gmtslc=deltastr->Norm(NORM_TWO);
	cumgmtslc+=gmtslc;
	cumgmslc=cumbslc+cumgmtslc;

	femmodel->parameters->SetParam(cumgmtslc,CumGmtslcEnum);
	femmodel->parameters->SetParam(cumgmslc,CumGmslcEnum);

	/*Outputs some metrics:*/
	femmodel->parameters->FindParam(&step,StepEnum);
	femmodel->parameters->FindParam(&time,TimeEnum);

	femmodel->results->AddResult(new GenericExternalResult<IssmDouble>(femmodel->results->Size()+1,GmtslcEnum,gmtslc,step,time));
	femmodel->results->AddResult(new GenericExternalResult<IssmDouble>(femmodel->results->Size()+1,CumGmtslcEnum,cumgmtslc,step,time));
	femmodel->results->AddResult(new GenericExternalResult<IssmDouble>(femmodel->results->Size()+1,CumGmslcEnum,cumgmslc,step,time));

	/*update element inputs:*/
	InputUpdateFromVectorx(femmodel,sealevel,SealevelEnum,VertexSIdEnum);	

	/*Free resources:*/	
	delete sealevel;
	delete deltadsl;
	delete deltastr;
}
/*}}}*/
void              coupleroutput_core(FemModel* femmodel){  /*{{{*/

	/*parameters:*/
	int iscoupling;
	int horiz=0;
	int count, frequency;

	/*retrieve more parameters:*/
	femmodel->parameters->FindParam(&iscoupling,IsSlcCouplingEnum);
	femmodel->parameters->FindParam(&horiz,SolidearthSettingsHorizEnum);
	femmodel->parameters->FindParam(&count,SealevelchangeRunCountEnum);
	femmodel->parameters->FindParam(&frequency,SolidearthSettingsRunFrequencyEnum);

	count++;
	femmodel->parameters->SetParam(count,SealevelchangeRunCountEnum); 

	if(iscoupling){
		/*transfer sea level back to ice caps:*/
		TransferSealevel(femmodel,SealevelEnum);
		TransferSealevel(femmodel,BedEnum);
		if(horiz){
			TransferSealevel(femmodel,BedNorthEnum);
			TransferSealevel(femmodel,BedEastEnum);
		}
	}

}; /*}}}*/
void              ivins_deformation_core(FemModel* femmodel){ /*{{{*/

	int  gsize;
	Vector<IssmDouble> *bedup  = NULL; 
	Vector<IssmDouble> *beduprate= NULL; 
	IssmDouble          *xx     = NULL;
	IssmDouble          *yy     = NULL;

	if(VerboseSolution()) _printf0_("	  computing vertical deformation using Ivins model. \n");

	/*find size of vectors:*/
	gsize      = femmodel->nodes->NumberOfDofs(GsetEnum);

	/*Find the litho material to be used by all the elements:*/
	Matlitho* matlitho=NULL;
	for (Object* & object: femmodel->materials->objects){
		Material* material=xDynamicCast<Material*>(object);
		if(material->ObjectEnum()==MatlithoEnum){
			matlitho=xDynamicCast<Matlitho*>(material);
			break;
		}
	}

	/*initialize vectors:*/
	bedup = new Vector<IssmDouble>(gsize);
	beduprate = new Vector<IssmDouble>(gsize);

	/*retrieve geometric information: */
	VertexCoordinatesx(&xx,&yy,NULL,femmodel->vertices); 

	/*Go through elements, and add contribution from each element to the deflection vector wg:*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->GiaDeflection(bedup,beduprate, matlitho, xx,yy);
	}

	/*Assemble parallel vector:*/
	beduprate->Assemble();
	bedup->Assemble();

	/*Save results:*/
	InputUpdateFromVectorx(femmodel,bedup,BedGRDEnum,VertexSIdEnum);
	femmodel->inputs->AXPY(1,BedGRDEnum,BedEnum);

	/*Free resources: */
	xDelete<IssmDouble>(xx);
	xDelete<IssmDouble>(yy);
	delete beduprate;
	delete bedup;
}
/*}}}*/
void              sealevelchange_initialgeometry(FemModel* femmodel) {  /*{{{*/

	/*Geometry core where we compute geometrical kernels and weights:*/

	/*parameters: */
	IssmDouble *xxe    = NULL;
	IssmDouble *yye    = NULL;
	IssmDouble *zze    = NULL;
	IssmDouble* areae  = NULL;
	int  nel;
	int* lids=NULL;
	int* n_activevertices=NULL;
	int  grdmodel=0;
	bool geometrydone;

	/*retrieve parameters:*/
	femmodel->parameters->FindParam(&grdmodel,GrdModelEnum);
	nel=femmodel->elements->NumberOfElements();

	/*did we already do this? if so, skip :*/
	femmodel->parameters->FindParam(&geometrydone,SealevelchangeGeometryDoneEnum);
	if (geometrydone){
		if(VerboseSolution()) _printf0_("	  initial sea level geometrical already computed, skipping.\n");
		return;
	}

	/*early return?:*/
	if(grdmodel!=ElasticEnum) return;

	/*Verbose: */
	if(VerboseSolution()) _printf0_("	  computing initial sea level geometrical kernels and weights.\n");

	/*recover x,y,z and areas from elements: */
	ElementCoordinatesx(&xxe,&yye,&zze,&areae,femmodel->elements);

	/*Compute element ids, used to speed up computations in convolution phase:{{{*/
	lids=xNew<int>(femmodel->vertices->Size());
	n_activevertices = xNew<int>(nel);
	//initialize lids to -1, vertex count to 3
	for (int v=0; v<femmodel->vertices->Size();v++) lids[v]=-1;
	for (int e=0; e<nel;e++) n_activevertices[e]=3;

	for(Object* & object : femmodel->elements->objects){
		Element*   element=xDynamicCast<Element*>(object);
		for(int i=0;i<3;i++){
			// if lids where we are looking points to an element id (.i.e. not -1) then we are about to claim that element's vertex
			// and need to lower the number of vertices it is in charge of
			if (lids[element->vertices[i]->lid] !=-1){
				n_activevertices[lids[element->vertices[i]->lid]]-=1;
			}
			lids[element->vertices[i]->lid]=element->lid;
		}
	}

	/*}}}*/

	/*Run sealevel geometry routine in elements:*/
	for(Object* & object : femmodel->elements->objects){
		Element*   element=xDynamicCast<Element*>(object);
		element->SealevelchangeGeometryInitial(xxe,yye,zze,areae,lids,n_activevertices);
	}

	femmodel->parameters->AddObject(new DoubleVecParam(XxeEnum,xxe,nel));
	femmodel->parameters->AddObject(new DoubleVecParam(YyeEnum,yye,nel));
	femmodel->parameters->AddObject(new DoubleVecParam(ZzeEnum,zze,nel));
	femmodel->parameters->AddObject(new DoubleVecParam(AreaeEnum,areae,nel));

	#ifdef _ISSM_DEBUG_
	femmodel->results->AddResult(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,XxeEnum,xxe,nel,1,1,1));
	femmodel->results->AddResult(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,YyeEnum,yye,nel,1,1,1));
	femmodel->results->AddResult(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,ZzeEnum,zze,nel,1,1,1));
	femmodel->results->AddResult(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,AreaeEnum,areae,nel,1,1,1));
	#endif
	
	geometrydone=true;
	femmodel->parameters->SetParam(geometrydone,SealevelchangeGeometryDoneEnum);

	xDelete<IssmDouble>(xxe);
	xDelete<IssmDouble>(yye);
	xDelete<IssmDouble>(zze);
	xDelete<IssmDouble>(areae);
	xDelete<int>(lids);
	xDelete<int>(n_activevertices);

	return;

}/*}}}*/
SealevelGeometry* sealevelchange_geometry(FemModel* femmodel) {  /*{{{*/

	/*Geometry core where we compute updates to the Green function kernels and weights, dependent 
	 * on the evolution of levelsets: */

	/*parameters: */
	IssmDouble *xxe    = NULL;
	IssmDouble *yye    = NULL;
	IssmDouble *zze    = NULL;
	IssmDouble* areae  = NULL;

	int nel;
	int grdmodel=0;
	int isgrd=0;
	int count, frequency;
	SealevelGeometry* slgeom=NULL;

	/*early return?:*/
	femmodel->parameters->FindParam(&grdmodel,GrdModelEnum);
	femmodel->parameters->FindParam(&isgrd,SolidearthSettingsGRDEnum);
	femmodel->parameters->FindParam(&count,SealevelchangeRunCountEnum);
	femmodel->parameters->FindParam(&frequency,SolidearthSettingsRunFrequencyEnum);
	if(grdmodel!=ElasticEnum || !isgrd) return NULL;
	if(count!=frequency) return NULL;

	/*retrieve parameters:*/
	femmodel->parameters->FindParam(&xxe,&nel,XxeEnum);
	femmodel->parameters->FindParam(&yye,&nel,YyeEnum);
	femmodel->parameters->FindParam(&zze,&nel,ZzeEnum);
	femmodel->parameters->FindParam(&areae,&nel,AreaeEnum);

	/*initialize SealevelloadMasks structure: */
	slgeom=new SealevelGeometry(femmodel->elements->Size(),femmodel->vertices->Size());

	/*Verbose: */
	if(VerboseSolution()) _printf0_("	  computing sea level geometrical kernel and weight updates.\n");

	/*Run sealevel geometry routine for elements with full loading:*/
	for(Object* & object : femmodel->elements->objects){
		Element*   element=xDynamicCast<Element*>(object);
		element->SealevelchangeGeometryCentroidLoads(slgeom,xxe,yye,zze,areae);
	}

	/*Initialize fractional loading mapping: */
	slgeom->InitializeMappingsAndBarycentres();

	/*Run sealevel geometry routine for elements with fractional loading:*/
	for(Object* & object : femmodel->elements->objects){
		Element*   element=xDynamicCast<Element*>(object);
		element->SealevelchangeGeometrySubElementLoads(slgeom,areae);
	}

	/*Assemble barycentres of fraction loading elements:*/
	slgeom->Assemble();

	/*Create fractional green function kernels: */
	for(Object* & object : femmodel->elements->objects){
		Element*   element=xDynamicCast<Element*>(object);
		element->SealevelchangeGeometrySubElementKernel(slgeom);
	}

	/*Compute spherical harmonic functions for spatial integrations of the loads*/
	slgeom->BuildSphericalHarmonics();

	femmodel->parameters->AddObject(new DoubleVecParam(XxeEnum,xxe,nel));
	femmodel->parameters->AddObject(new DoubleVecParam(YyeEnum,yye,nel));
	femmodel->parameters->AddObject(new DoubleVecParam(ZzeEnum,zze,nel));
	femmodel->parameters->AddObject(new DoubleVecParam(AreaeEnum,areae,nel));

	/*Free resources:*/
	xDelete<IssmDouble>(xxe);
	xDelete<IssmDouble>(yye);
	xDelete<IssmDouble>(zze);
	xDelete<IssmDouble>(areae);

	return slgeom;

}/*}}}*/
void              sealevelchange_finalize(FemModel* femmodel) {  /*{{{*/

	bool isuq=false;
	
	BarystaticContributions* barycontrib=NULL;
	GenericParam<BarystaticContributions*>* barycontribparam=NULL;
	
	femmodel->parameters->FindParam(&isuq,QmuIsdakotaEnum);

	if(isuq){
		//reset barycontrib object:
		barycontribparam = xDynamicCast<GenericParam<BarystaticContributions*>*>(femmodel->parameters->FindParamObject(BarystaticContributionsEnum));
		barycontrib=barycontribparam->GetParameterValue(); 
		barycontrib->Finalize();
	}
	else {
		/*Erase barycontrib object: */
		barycontribparam = xDynamicCast<GenericParam<BarystaticContributions*>*>(femmodel->parameters->FindParamObject(BarystaticContributionsEnum));
		barycontrib=barycontribparam->GetParameterValue();
		delete barycontrib;
	}

	return;

}/*}}}*/


void slc_geometry_cleanup(SealevelGeometry* slgeom, FemModel* femmodel){  /*{{{*/

	/*early return?:*/
	if(slgeom==NULL) return;

	int horiz;
	femmodel->parameters->FindParam(&horiz,SolidearthSettingsHorizEnum);
	for (int l=0;l<SLGEOM_NUMLOADS;l++){
		femmodel->inputs->DeleteInput(slgeom->AlphaIndexEnum(l));
		if(horiz) femmodel->inputs->DeleteInput(slgeom->AzimuthIndexEnum(l));
	}

	delete slgeom;
} /*}}}*/

/*subroutines:*/
bool slcconvergence(IssmDouble* RSLg,IssmDouble* RSLg_old,IssmDouble eps_rel,IssmDouble eps_abs, IssmDouble totaloceanarea, FemModel* femmodel){ /*{{{*/

	int nel;
	bool converged=true;
	IssmDouble ndS,nS, nS_old; 
	IssmDouble* dRSLg    = NULL;
	IssmDouble rho_water =0;

	femmodel->parameters->FindParam(&nel,MeshNumberofelementsEnum);
	femmodel->parameters->FindParam(&rho_water,MaterialsRhoSeawaterEnum);

	//compute norm(du) and norm(u) if requested
	dRSLg=xNewZeroInit<IssmDouble>(nel);

	ndS=0;
	nS=0;
	nS_old=0;

	for (int e=0;e<nel;e++){
		dRSLg[e]=(RSLg[e]-RSLg_old[e])/rho_water/totaloceanarea;
		ndS+=pow(dRSLg[e],2.0);
		nS+=pow(RSLg[e]/rho_water/totaloceanarea,2.0);
		nS_old+=pow(RSLg_old[e]/rho_water/totaloceanarea,2.0);
	}

	if (xIsNan<IssmDouble>(ndS)){
		_error_("convergence criterion is NaN (RSL_old=" << nS_old << " RSL=" << nS << ")");
	}

	if(!xIsNan<IssmDouble>(eps_rel)){
		if (xIsNan<IssmDouble>(nS_old)) _error_("convergence criterion is NaN! (check the initial RSL)");
	}

	//clean up
	xDelete<IssmDouble>(dRSLg);

	//print
	if(!xIsNan<IssmDouble>(eps_rel)){
		if((ndS/nS)<eps_rel){
			if(VerboseConvergence()) _printf0_(setw(50) << left << "              convergence criterion: norm(dS)/norm(S)" << ndS/nS*100 << " < " << eps_rel*100 << " %\n");
		}
		else{ 
			if(VerboseConvergence()) _printf0_(setw(50) << left << "              convergence criterion: norm(dS)/norm(S)" << ndS/nS*100 << " > " << eps_rel*100 << " %\n");
			converged=false;
		}
	}
	if(!xIsNan<IssmDouble>(eps_abs)){
		if(ndS<eps_abs){
			if(VerboseConvergence()) _printf0_(setw(50) << left << "              convergence criterion: norm(dS)" << ndS << " < " << eps_abs << " \n");
		}
		else{ 
			if(VerboseConvergence()) _printf0_(setw(50) << left << "              convergence criterion: norm(dS)" << ndS << " > " << eps_abs << " \n");
			converged=false;
		}
	}

	/*assign output*/
	return converged;

} /*}}}*/
IssmDouble  SealevelloadsOceanAverage(GrdLoads* loads, Vector<IssmDouble>* oceanareas, Vector<IssmDouble>* suboceanareas, IssmDouble totaloceanarea){ /*{{{*/

	IssmDouble sealevelloadsaverage;	
	IssmDouble subsealevelloadsaverage;	

	loads->vsealevelloads->Sum(&sealevelloadsaverage);
	loads->vsubsealevelloads->Sum(&subsealevelloadsaverage);

	return (sealevelloadsaverage+subsealevelloadsaverage)/totaloceanarea;
} /*}}}*/
void PolarMotion(IssmDouble* polarmotionvector, FemModel* femmodel,GrdLoads* loads, SealevelGeometry* slgeom, bool computefuture){ /*{{{*/
	//The purpose of this routine is to get the polar motion vector m=(m1, m2, m3) induced by the GrdLoads
	IssmDouble  S2coef[3];
	IssmDouble*	pmtf_col= NULL;
	IssmDouble*	pmtf_ortho   = NULL;
	IssmDouble*	pmtf_z   = NULL;
	IssmDouble* m1=NULL;
	IssmDouble* m2=NULL;
	IssmDouble* m3=NULL;
	IssmDouble* m1interp=NULL;
	IssmDouble* m2interp=NULL;
	IssmDouble* m3interp=NULL;

	IssmDouble  moi_e, moi_p, re;
	IssmDouble mhprefactor, mzprefactor;
	bool rotation=false;
	bool viscous=false;
	int nt=1;

	IssmDouble* viscoustimes=NULL;
	IssmDouble* viscouspolarmotion=NULL;
	int         viscousnumsteps;
	int         viscousindex=0; 
	int         dummy;
	IssmDouble  currenttime, final_time, lincoeff, timeacc;

	/*early return?:*/
	femmodel->parameters->FindParam(&rotation,SolidearthSettingsRotationEnum);
	if(!rotation)return;

	/*retrieve parameters: */
	femmodel->parameters->FindParam(&viscous,SolidearthSettingsViscousEnum);
	femmodel->parameters->FindParam(&moi_e,RotationalEquatorialMoiEnum);
	femmodel->parameters->FindParam(&moi_p,RotationalPolarMoiEnum);
	femmodel->parameters->FindParam(&pmtf_col,NULL,SealevelchangePolarMotionTransferFunctionColinearEnum);
	femmodel->parameters->FindParam(&pmtf_ortho,NULL,SealevelchangePolarMotionTransferFunctionOrthogonalEnum);
	femmodel->parameters->FindParam(&pmtf_z,NULL,SealevelchangePolarMotionTransferFunctionZEnum);
	femmodel->parameters->FindParam(&re,SolidearthPlanetRadiusEnum);

	if (viscous){
		femmodel->parameters->FindParam(&final_time,TimesteppingFinalTimeEnum);
		femmodel->parameters->FindParam(&timeacc,SolidearthSettingsTimeAccEnum);
		femmodel->parameters->FindParam(&viscousnumsteps,SealevelchangeViscousNumStepsEnum);
		femmodel->parameters->FindParam(&viscoustimes,NULL,SealevelchangeViscousTimesEnum);
		femmodel->parameters->FindParam(&viscousindex,SealevelchangeViscousIndexEnum);
		femmodel->parameters->FindParam(&currenttime,TimeEnum);
		femmodel->parameters->FindParam(&viscouspolarmotion,NULL,NULL,SealevelchangeViscousPolarMotionEnum);

		if (computefuture){
			nt = viscousnumsteps;
			m1interp=xNewZeroInit<IssmDouble>(nt);
			m2interp=xNewZeroInit<IssmDouble>(nt);
			m3interp=xNewZeroInit<IssmDouble>(nt);
		}
	}

	m1=xNewZeroInit<IssmDouble>(nt);
	m2=xNewZeroInit<IssmDouble>(nt);
	m3=xNewZeroInit<IssmDouble>(nt);

	//Isolate degree 2 load coefficients
	for (int i=0;i<3;i++) S2coef[i]=0;
	loads->SHDegree2Coefficients(&S2coef[0],femmodel,slgeom);

	//compute present (& future) polar motion from present loads
	mhprefactor=-4.*M_PI/5.*pow(re,4.)/(moi_p-moi_e);
	mzprefactor=8.*M_PI/15.*pow(re,4.)/moi_p;
	for (int tprime=0;tprime<nt;tprime++){
		m1[tprime] = mhprefactor * (pmtf_col[tprime] * S2coef[1] + pmtf_ortho[tprime] * S2coef[2]); //x-component
		m2[tprime] = mhprefactor * (pmtf_col[tprime] * S2coef[2] - pmtf_ortho[tprime] * S2coef[1]); //y-component
		m3[tprime] = mzprefactor *  pmtf_z[tprime]   * S2coef[0];				    //z-component
	}

	if(viscous){
		// we need to do up to 3 things (* = only if computefuture)
		// 1*: add new PM contribution to the viscous stack for future time steps
		// 2: collect viscous PM from past loads due at present-day and add it to PM[current_time]
		// 3*: subtract from viscous stack PM that has already been accounted for so we don't add it again at the next time step
		if(computefuture){ 
			if(viscoustimes[viscousindex]<final_time){
				lincoeff=(viscoustimes[viscousindex+1]-viscoustimes[viscousindex])/timeacc;
				for(int t=viscousindex;t<nt;t++){ //we resynchronize m from the relative time above to the absolute time where t=0 <=> beginning of the simulation
					if(t==viscousindex){
						m1interp[t]=  m1[0];
						m2interp[t]=  m2[0];
						m3interp[t]=  m3[0];
					}
					else{ //we reinterpolate PM on viscoustimes, so we can handle the case where we are running with adaptative/uneven time steps
						int tprime=t-viscousindex-1;
						m1interp[t]=  (1.0-lincoeff)*m1[tprime]+lincoeff*m1[tprime+1];
						m2interp[t]=  (1.0-lincoeff)*m2[tprime]+lincoeff*m2[tprime+1];
						m3interp[t]=  (1.0-lincoeff)*m3[tprime]+lincoeff*m3[tprime+1];
					}
				}
			}
		}
		/*update PM at present time using viscous stack at present time: */
		m1[0]+=viscouspolarmotion[0*nt+viscousindex];
		m2[0]+=viscouspolarmotion[1*nt+viscousindex]; 
		m3[0]+=viscouspolarmotion[2*nt+viscousindex]; 
		if(computefuture){ /*update viscous stack with future deformation from present load: */
			for(int t=nt-1;t>=viscousindex;t--){
				//offset viscous PM to remove all deformation that has already been added
				viscouspolarmotion[0*nt+t]+=m1interp[t]-m1interp[viscousindex]-viscouspolarmotion[0*nt+viscousindex];
				viscouspolarmotion[1*nt+t]+=m2interp[t]-m2interp[viscousindex]-viscouspolarmotion[1*nt+viscousindex];
				viscouspolarmotion[2*nt+t]+=m3interp[t]-m3interp[viscousindex]-viscouspolarmotion[2*nt+viscousindex];
			}
			// save updated viscous PM
			femmodel->parameters->SetParam(viscouspolarmotion,viscousnumsteps,3,SealevelchangeViscousPolarMotionEnum);
		}
	}

	/*Assign output pointers:*/
	polarmotionvector[0]=m1[0];
	polarmotionvector[1]=m2[0];
	polarmotionvector[2]=m3[0];

	/*Free allocations:*/
	xDelete<IssmDouble>(m1);
	xDelete<IssmDouble>(m2);
	xDelete<IssmDouble>(m3);
	if (viscous){
		if (computefuture){
			xDelete<IssmDouble>(m1interp);
			xDelete<IssmDouble>(m2interp);
			xDelete<IssmDouble>(m3interp);
		}
		xDelete<IssmDouble>(viscoustimes);
		xDelete<IssmDouble>(viscouspolarmotion);
	}
	xDelete<IssmDouble>(pmtf_col);
	xDelete<IssmDouble>(pmtf_ortho);
	xDelete<IssmDouble>(pmtf_z);

} /*}}}*/
void       SealevelchangeUpdateViscousTimeSeries(FemModel* femmodel){ /*{{{*/

	IssmDouble* viscouspolarmotion=NULL;
	IssmDouble* viscoustimes=NULL;
	int         viscousnumsteps;
	int         viscousindex=0; 
	int         newindex=0; 
	int         dummy;
	bool        viscous=false;
	bool        rotation=false;
	IssmDouble  currenttime;
	IssmDouble  lincoeff=0;

	femmodel->parameters->FindParam(&viscous,SolidearthSettingsViscousEnum);
	femmodel->parameters->FindParam(&rotation,SolidearthSettingsRotationEnum);

	if(viscous){
		femmodel->parameters->FindParam(&viscousnumsteps,SealevelchangeViscousNumStepsEnum);
		femmodel->parameters->FindParam(&viscoustimes,NULL,SealevelchangeViscousTimesEnum);
		femmodel->parameters->FindParam(&viscousindex,SealevelchangeViscousIndexEnum);
		femmodel->parameters->FindParam(&currenttime,TimeEnum);
		if (rotation) 	femmodel->parameters->FindParam(&viscouspolarmotion,NULL,NULL,SealevelchangeViscousPolarMotionEnum);

		bool foundtime=false;
		int offset=1; //handles the egde case where time found = max time in viscoustimes
		lincoeff=0;
		newindex=viscousnumsteps-2;

		for(int t=viscousindex;t<viscousnumsteps;t++){
			if (viscoustimes[t]>=currenttime){
				newindex=t-1;
				foundtime=true;
				lincoeff=(currenttime-viscoustimes[newindex])/(viscoustimes[t]-viscoustimes[newindex]);
				offset=0;
				break;
			}
		}

		if(rotation){
			int index=0;
			for (int i=0;i<3;i++){
				index=i*viscousnumsteps+newindex;
				viscouspolarmotion[index+offset]=(1-lincoeff)*viscouspolarmotion[index]+lincoeff*viscouspolarmotion[index+1];
			}
			femmodel->parameters->SetParam(viscouspolarmotion,viscousnumsteps,3,SealevelchangeViscousPolarMotionEnum);
		}

		/*update viscous inputs:*/
		for(Object* & object : femmodel->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			element->SealevelchangeUpdateViscousFields(lincoeff,newindex,offset);
		}

		viscoustimes[newindex]=currenttime;
		viscousindex=newindex+offset;

		femmodel->parameters->SetParam(viscousindex,SealevelchangeViscousIndexEnum);
		femmodel->parameters->SetParam(viscoustimes,viscousnumsteps,SealevelchangeViscousTimesEnum);

		/*free allocations:*/
		xDelete<IssmDouble>(viscoustimes);
		if (rotation) 	xDelete<IssmDouble>(viscouspolarmotion);
	}

}
void        ConserveOceanMass(FemModel* femmodel,GrdLoads* loads, IssmDouble offset, SealevelGeometry* slgeom){ /*{{{*/

	/*Shift sealevel loads by ocean average, only on ocean! :*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->SealevelchangeShift(loads, offset,slgeom);
	}
	loads->AssembleSealevelLoads();

} /*}}}*/
IssmDouble* CombineLoads(IssmDouble* load,IssmDouble* subload,FemModel* femmodel, SealevelGeometry* slgeom,int loadtype,int nel){ /*{{{*/

	//merges loads on centroids and subelements onto a single variable loadcopy
	int* indices=xNew<int>(nel);
	for(int i=0;i<nel;i++)indices[i]=i;

	Vector<IssmDouble>* vloadcopy=new Vector<IssmDouble>(nel);
	IssmDouble* loadcopy=xNew<IssmDouble>(nel);

	vloadcopy->SetValues(nel,indices,load,INS_VAL);
	vloadcopy->Assemble();

	if(subload){
		for (int i=0;i<femmodel->elements->Size();i++){
			if (slgeom->issubelement[loadtype][i]){
				int se= slgeom->subelementmapping[loadtype][i];
				IssmDouble subloadi=subload[se];
				Element* element=dynamic_cast<Element*>(femmodel->elements->GetObjectByOffset(i));
				vloadcopy->SetValue(element->Sid(),subloadi,ADD_VAL);
			}
		}
	}
	vloadcopy->Assemble();
	loadcopy=vloadcopy->ToMPISerial();

	return loadcopy;

} /*}}}*/

/*Coupling routines:*/
void TransferForcing(FemModel* femmodel,int forcingenum){ /*{{{*/

	/*forcing being transferred from models to earth: */
	IssmDouble** forcings=NULL;
	IssmDouble*  forcing=NULL; 
	Vector<IssmDouble>* forcingglobal=NULL; 
	IssmDouble* transfercount=NULL; 
	int*         nvs=NULL;
   int modelid,earthid,nummodels;

	/*transition vectors:*/
	IssmDouble** transitions=NULL;
	int          ntransitions; 
	int*         transitions_m=NULL;
	int*         transitions_n=NULL;
	int          nv;
	int          existforcing=0;         

	/*communicators:*/
	ISSM_MPI_Comm    tocomm;
	ISSM_MPI_Comm   *fromcomms = NULL;
	ISSM_MPI_Status  status;
	ISSM_MPI_Request send_request_1=ISSM_MPI_REQUEST_NULL;
   ISSM_MPI_Request send_request_2=ISSM_MPI_REQUEST_NULL;
	ISSM_MPI_Request send_request_3=ISSM_MPI_REQUEST_NULL;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&modelid,ModelIdEnum);
	femmodel->parameters->FindParam(&earthid,EarthIdEnum);
	femmodel->parameters->FindParam(&nummodels,NumModelsEnum);
	int my_rank=IssmComm::GetRank();


	/*retrieve the inter communicators that will be used to send data from each ice cap to the earth: */
	if(modelid==earthid){
		GenericParam<ISSM_MPI_Comm*>* parcoms = dynamic_cast<GenericParam<ISSM_MPI_Comm*>*>(femmodel->parameters->FindParamObject(IcecapToEarthCommEnum));
		if(!parcoms)_error_("TransferForcing error message: could not find IcecapToEarthComm communicator");
		fromcomms=parcoms->GetParameterValue();
	}
	else {
		GenericParam<ISSM_MPI_Comm>* parcom = dynamic_cast<GenericParam<ISSM_MPI_Comm>*>(femmodel->parameters->FindParamObject(IcecapToEarthCommEnum));
		if(!parcom)_error_("TransferForcing error message: could not find IcecapToEarthComm communicator");
		tocomm=parcom->GetParameterValue();
	}

	/*For each icecap, retrieve the forcing vector that will be sent to the earth model: */
	if(modelid!=earthid){
		nv=femmodel->vertices->NumberOfVertices();
		existforcing=reCast<int>(femmodel->inputs->Exist(forcingenum));
		if(existforcing){
			GetVectorFromInputsx(&forcing,femmodel,forcingenum,VertexSIdEnum);
			GetVectorFromInputsx(&transfercount,femmodel,CouplingTransferCountEnum,VertexSIdEnum);
			for (int i=0;i<nv;i++) forcing[i]/=transfercount[i]; //Divide forcing at this vertex by the number of icecaps that share it. This way we average the forcing when adding it into the earth model.
		}
	}

	/*Send the forcing to the earth model*/
	if(my_rank==0){
		if(modelid==earthid){
			forcings=xNew<IssmDouble*>(nummodels-1);
			nvs=xNew<int>(nummodels-1);
			for(int i=0;i<earthid;i++){
				ISSM_MPI_Recv(&existforcing, 1, ISSM_MPI_INT, 0,i, fromcomms[i], &status);
				if(existforcing){
					ISSM_MPI_Recv(nvs+i, 1, ISSM_MPI_INT, 0,i, fromcomms[i], &status);
					forcings[i]=xNew<IssmDouble>(nvs[i]);
					ISSM_MPI_Recv(forcings[i], nvs[i], ISSM_MPI_DOUBLE, 0,i, fromcomms[i], &status);
				}
				else{
					forcings[i]=NULL;
				}
			}

		}
		else{
			ISSM_MPI_Isend(&existforcing, 1, ISSM_MPI_INT, 0, modelid, tocomm,&send_request_1);
			if(existforcing){
				ISSM_MPI_Isend(&nv, 1, ISSM_MPI_INT, 0, modelid, tocomm, &send_request_2);
				ISSM_MPI_Isend(forcing, nv, ISSM_MPI_DOUBLE, 0, modelid, tocomm, &send_request_3);
			}
		}
	}

	/*On the earth model, consolidate all the forcings into one, and update the elements dataset accordingly: {{{*/
	if(modelid==earthid){

		/*Out of all the delta thicknesses, build one delta thickness vector made of all the ice cap contributions. 
		 *First, build the global delta thickness vector in the earth model: */
		nv=femmodel->vertices->NumberOfVertices();
		GetVectorFromInputsx(&forcingglobal,femmodel,forcingenum,VertexSIdEnum);

		forcingglobal->Set(0.);

		/*Retrieve transition vectors, used to plug from each ice cap into the global forcing:*/
		femmodel->parameters->FindParam(&transitions,&ntransitions,&transitions_m,&transitions_n,SealevelchangeTransitionsEnum);

		if(ntransitions!=earthid)_error_("TransferForcing error message: number of transition vectors is not equal to the number of icecaps!");

		/*Go through all the delta thicknesses coming from each ice cap: */
		if(my_rank==0){
			for(int i=0;i<earthid;i++){

				IssmDouble* forcingfromcap= forcings[i]; //careful, this only exists on rank 0 of the earth model!
				if(forcingfromcap){
					IssmDouble* transition=transitions[i];
					int         M=transitions_m[i];

					/*build index to plug values: */
					int* index=xNew<int>(M); for(int i=0;i<M;i++)index[i]=reCast<int>(transition[i])-1; //matlab indexing!

					/*We are going to plug this vector into the earth model, at the right vertices corresponding to this particular 
					 * ice cap: */

					forcingglobal->SetValues(M,index,forcingfromcap,ADD_VAL);
					xDelete<int>(index);
				}
			}
		}

		/*Assemble vector:*/
		forcingglobal->Assemble();

		/*Plug into elements:*/
		InputUpdateFromVectorx(femmodel,forcingglobal,forcingenum,VertexSIdEnum);
	} 

	/*Free resources:*/
	if(my_rank==0 && modelid!=earthid){
		ISSM_MPI_Wait(&send_request_1,&status);
		if(existforcing){
			ISSM_MPI_Wait(&send_request_2,&status);
			ISSM_MPI_Wait(&send_request_3,&status);
		}
	}
	if(forcings){
		for(int i=0;i<nummodels-1;i++) xDelete<IssmDouble>(forcings[i]);
		xDelete<IssmDouble*>(forcings);
	}
	if(forcing)xDelete<IssmDouble>(forcing);
	if(transfercount) xDelete<IssmDouble>(transfercount);
	if(forcingglobal)delete forcingglobal;
	if(transitions){
		for(int i=0;i<earthid;i++){
			IssmDouble* temp=transitions[i];
			xDelete<IssmDouble>(temp);
		}
		xDelete<IssmDouble*>(transitions);
		xDelete<int>(transitions_m);
		xDelete<int>(transitions_n);
	}
	if(nvs)xDelete<int>(nvs);
} /*}}}*/
void TransferSealevel(FemModel* femmodel,int forcingenum){ /*{{{*/

	/*forcing being transferred from earth to ice caps: */
	IssmDouble*  forcing=NULL; 
	IssmDouble*  forcingglobal=NULL; 

	/*transition vectors:*/
	IssmDouble** transitions=NULL;
	int          ntransitions; 
	int*         transitions_m=NULL;
	int*         transitions_n=NULL;
	int  nv;
	int  modelid,earthid,nummodels;
	int  numcoms;

	/*communicators:*/
	ISSM_MPI_Comm fromcomm;
	ISSM_MPI_Comm* tocomms=NULL;
	ISSM_MPI_Status status;
	ISSM_MPI_Request* send_requests_1=NULL;
	ISSM_MPI_Request* send_requests_2=NULL;


	/*Recover some parameters: */
	femmodel->parameters->FindParam(&modelid,ModelIdEnum);
	femmodel->parameters->FindParam(&earthid,EarthIdEnum);
	femmodel->parameters->FindParam(&nummodels,NumModelsEnum);
	int my_rank=IssmComm::GetRank();

	/*retrieve the inter communicators that will be used to send data from earth to ice caps:*/
	if(modelid==earthid){
		GenericParam<ISSM_MPI_Comm*>* parcoms = dynamic_cast<GenericParam<ISSM_MPI_Comm*>*>(femmodel->parameters->FindParamObject(IcecapToEarthCommEnum));
		if(!parcoms)_error_("TransferSealevel error message: could not find IcecapToEarthComm communicator");
		tocomms=parcoms->GetParameterValue();
		//femmodel->parameters->FindParam((int**)(&tocomms),&numcoms,IcecapToEarthCommEnum);
	}
	else{
		GenericParam<ISSM_MPI_Comm>* parcom = dynamic_cast<GenericParam<ISSM_MPI_Comm>*>(femmodel->parameters->FindParamObject(IcecapToEarthCommEnum));
		if(!parcom)_error_("TransferSealevel error message: could not find IcecapToEarthComm communicator");
		fromcomm=parcom->GetParameterValue();
		//femmodel->parameters->FindParam((int*)(&fromcomm), IcecapToEarthCommEnum);
	}

	/*Retrieve sea-level on earth model: */
	if(modelid==earthid){
		nv=femmodel->vertices->NumberOfVertices();
		GetVectorFromInputsx(&forcingglobal,femmodel,forcingenum,VertexSIdEnum);
	}

	/*Send the forcing to the ice caps:{{{*/
	if(my_rank==0){

		if(modelid==earthid){

			/*Retrieve transition vectors, used to figure out global forcing contribution to each ice cap's own elements: */
			femmodel->parameters->FindParam(&transitions,&ntransitions,&transitions_m,&transitions_n,SealevelchangeTransitionsEnum);
			if(ntransitions!=earthid) _error_("TransferSealevel error message: number of transition vectors is not equal to the number of icecaps!");

			/*Prepare requests*/
			send_requests_1 = xNew<ISSM_MPI_Request>(earthid);
			send_requests_2 = xNew<ISSM_MPI_Request>(earthid);
			for(int i=0;i<earthid;i++){
				send_requests_1[i] = ISSM_MPI_REQUEST_NULL;
				send_requests_2[i] = ISSM_MPI_REQUEST_NULL;
			}

			for(int i=0;i<earthid;i++){
				nv=transitions_m[i];
				forcing=xNew<IssmDouble>(nv);
				IssmDouble* transition=transitions[i];
				for(int j=0;j<nv;j++){
					forcing[j]=forcingglobal[reCast<int>(transition[j])-1];
				}
				ISSM_MPI_Isend(&nv, 1, ISSM_MPI_INT, 0, i, tocomms[i], &send_requests_1[i]);
				ISSM_MPI_Isend(forcing, nv, ISSM_MPI_DOUBLE, 0, i, tocomms[i], &send_requests_2[i]);
				xDelete<IssmDouble>(forcing);
			}
		}
		else{
			ISSM_MPI_Recv(&nv, 1, ISSM_MPI_INT, 0, modelid, fromcomm, &status);
			forcing=xNew<IssmDouble>(nv);
			ISSM_MPI_Recv(forcing, nv, ISSM_MPI_DOUBLE, 0, modelid, fromcomm, &status);
		}
	}
	/*}}}*/

	/*On each ice cap, spread the forcing across cpus, and update the elements dataset accordingly: {{{*/
	if(modelid!=earthid){

		ISSM_MPI_Bcast(&nv,1,ISSM_MPI_INT,0,IssmComm::GetComm());
		if(my_rank!=0)forcing=xNew<IssmDouble>(nv);
		ISSM_MPI_Bcast(forcing,nv,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

		/*Plug into elements:*/
		InputUpdateFromVectorx(femmodel,forcing,forcingenum,VertexSIdEnum);
	} 
	/*}}}*/

	/*Free resources:*/
	if(my_rank==0 && modelid==earthid){
		for(int i=0;i<earthid;i++){
			ISSM_MPI_Wait(&send_requests_1[i],&status);
			ISSM_MPI_Wait(&send_requests_2[i],&status);
		}
		xDelete<ISSM_MPI_Request>(send_requests_1);
		xDelete<ISSM_MPI_Request>(send_requests_2);
	}
	if(forcingglobal)xDelete<IssmDouble>(forcingglobal);
	if(forcing)xDelete<IssmDouble>(forcing);
	if(transitions){
		for(int i=0;i<ntransitions;i++){
			IssmDouble* temp=transitions[i];
			xDelete<IssmDouble>(temp);
		}
		xDelete<IssmDouble*>(transitions);
		xDelete<int>(transitions_m);
		xDelete<int>(transitions_n);
	}

} /*}}}*/
