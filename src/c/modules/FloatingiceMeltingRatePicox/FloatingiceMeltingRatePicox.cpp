/*!\file FloatingiceMeltingRatePicox
 * \brief: calculates Floating ice melting rate following the PICO model (Reese et al., 2017)
 */

#include "./FloatingiceMeltingRatePicox.h"
#include "../InputDuplicatex/InputDuplicatex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void FloatingiceMeltingRatePicox(FemModel* femmodel){/*{{{*/

	int maxbox;
	bool isplume;

	/*First, reset all melt to 0 */
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		int numvertices = element->GetNumberOfVertices();
		IssmDouble* values = xNewZeroInit<IssmDouble>(numvertices);
		element->AddInput(BasalforcingsFloatingiceMeltingRateEnum,values,P1Enum);
		xDelete<IssmDouble>(values);
	}

	/*PICO melt rate parameterization (Reese et al., 2018)*/
	femmodel->parameters->FindParam(&maxbox,BasalforcingsPicoMaxboxcountEnum);
	UpdateBoxIdsPico(femmodel);
	ComputeBoxAreasPico(femmodel);
	for(int i=0;i<maxbox;i++){
		UpdateBoxPico(femmodel,i);
		ComputeAverageOceanvarsPico(femmodel,i);
	}

	/*Optional buoyant plume melt rate parameterization (Lazeroms et al., 2018) */
	femmodel->parameters->FindParam(&isplume,BasalforcingsPicoIsplumeEnum);
	if(isplume) ComputeBasalMeltPlume(femmodel);
}/*}}}*/

void UpdateBoxIdsPico(FemModel* femmodel){/*{{{*/

	int         numvertices,num_basins,maxbox,basinid;
	IssmDouble  dist_max;
	IssmDouble* distances=NULL;

	femmodel->parameters->FindParam(&num_basins,BasalforcingsPicoNumBasinsEnum);
	femmodel->parameters->FindParam(&maxbox,BasalforcingsPicoMaxboxcountEnum);
	IssmDouble* dmax_basin_cpu=xNew<IssmDouble>(num_basins);

	InputDuplicatex(femmodel,MaskOceanLevelsetEnum,DistanceToGroundinglineEnum);
	femmodel->DistanceToFieldValue(MaskOceanLevelsetEnum,0.,DistanceToGroundinglineEnum);

	InputDuplicatex(femmodel,MaskIceLevelsetEnum,DistanceToCalvingfrontEnum);
	femmodel->DistanceToFieldValue(MaskIceLevelsetEnum,0.,DistanceToCalvingfrontEnum);

	/*find maximum distance to grounding line per domain and per basin*/
	IssmDouble maxdist_cpu=-1.;
	for(int i=0;i<num_basins;i++){dmax_basin_cpu[i]=-1;}
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		if(!element->IsIceInElement() || !element->IsAllFloating()) continue;
		numvertices = element->GetNumberOfVertices();
		distances=xNew<IssmDouble>(numvertices);
		element->GetInputListOnVertices(&distances[0],DistanceToGroundinglineEnum);
		element->GetInputValue(&basinid,BasalforcingsPicoBasinIdEnum);
		for(int k=0; k<numvertices; k++){
			if(fabs(distances[k])>maxdist_cpu){maxdist_cpu=fabs(distances[k]);}
			if(fabs(distances[k])>dmax_basin_cpu[basinid]){dmax_basin_cpu[basinid]=fabs(distances[k]);}
		}
		xDelete<IssmDouble>(distances);
	}

	/*Synchronize across cpus*/
	IssmDouble* dmax_basin=xNew<IssmDouble>(num_basins);
	ISSM_MPI_Allreduce((void*)&maxdist_cpu,(void*)&dist_max,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,IssmComm::GetComm());
	ISSM_MPI_Allreduce(dmax_basin_cpu,dmax_basin,num_basins,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,IssmComm::GetComm());

	/*Define maximum number of boxes per basin*/
	int* nd=xNew<int>(num_basins);
	for(int i=0; i<num_basins;i++){
		IssmDouble val=sqrt(dmax_basin[i]/dist_max)*(maxbox-1);

		#ifdef _HAVE_AD_
		_error_("Check the implementation of floor below");
		/*Do not use floor when AD is on*/
		int k=0; while(k<val+.5){k++;}
		nd[i]=k;

		#else
		nd[i]= reCast<int>(floor(val));
		#endif
	} 

	/*Assign box numbers*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->PicoUpdateBoxid(nd);
	}

	/*Cleanup and return */
	xDelete<int>(nd);
	xDelete<IssmDouble>(dmax_basin);
	xDelete<IssmDouble>(dmax_basin_cpu);

}/*}}}*/
void ComputeBoxAreasPico(FemModel* femmodel){/*{{{*/

	int num_basins,maxbox,basinid,boxid,domaintype;
	IssmDouble dist_max,area;

	femmodel->parameters->FindParam(&num_basins,BasalforcingsPicoNumBasinsEnum);
	femmodel->parameters->FindParam(&maxbox,BasalforcingsPicoMaxboxcountEnum);

	IssmDouble* boxareas=xNewZeroInit<IssmDouble>(num_basins*maxbox);

	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		if(!element->IsOnBase()) continue;
		Element* basalelement = element->SpawnBasalElement();
		if(!basalelement->IsIceInElement() || !basalelement->IsAllFloating()){
			if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
			continue;
		}
		basalelement->GetInputValue(&boxid,BasalforcingsPicoBoxIdEnum);
		basalelement->GetInputValue(&basinid,BasalforcingsPicoBasinIdEnum);
		boxareas[basinid*maxbox+boxid]+=basalelement->GetHorizontalSurfaceArea();
		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	}

	/*Synchronize across cpus*/
	IssmDouble* sumareas =xNew<IssmDouble>(num_basins*maxbox);
	ISSM_MPI_Allreduce(boxareas,sumareas,num_basins*maxbox,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
	//if(sumareas[0]==0){_error_("No elements in box 0, basal meltrates will be 0. Consider decreasing md.basalforcings.maxboxcount or refining your mesh!");}

	/*Update parameters to keep track of the new areas in future calculations*/
	femmodel->parameters->AddObject(new DoubleVecParam(BasalforcingsPicoBoxAreaEnum,sumareas,maxbox*num_basins));

	/*Cleanup and return */
	xDelete<IssmDouble>(boxareas);
	xDelete<IssmDouble>(sumareas);

}/*}}}*/
void UpdateBoxPico(FemModel* femmodel, int loopboxid){/*{{{*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->PicoUpdateBox(loopboxid);
	}
}/*}}}*/
void ComputeAverageOceanvarsPico(FemModel* femmodel, int boxid){/*{{{*/

	int num_basins, basinid, maxbox, M, domaintype;
	IssmDouble area, toc, soc, overturning;
	IssmDouble* boxareas=NULL;
	IssmDouble* overturning_weighted_avg=NULL;

	femmodel->parameters->FindParam(&num_basins,BasalforcingsPicoNumBasinsEnum);
	femmodel->parameters->FindParam(&maxbox,BasalforcingsPicoMaxboxcountEnum);
	femmodel->parameters->FindParam(&boxareas,&M, BasalforcingsPicoBoxAreaEnum);
	IssmDouble* toc_weighted_avg           = xNewZeroInit<IssmDouble>(num_basins);
	IssmDouble* soc_weighted_avg           = xNewZeroInit<IssmDouble>(num_basins);
	IssmDouble* toc_sumweightedavg         = xNewZeroInit<IssmDouble>(num_basins);
	IssmDouble* soc_sumweightedavg         = xNewZeroInit<IssmDouble>(num_basins);
	IssmDouble* overturning_sumweightedavg = xNewZeroInit<IssmDouble>(num_basins);

	/* Compute Toc and Soc weighted avg (boxes 0 to n-1) */
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);

		/*Check whether we should continue*/
		if(!element->IsOnBase()) continue;
		if(!element->IsIceInElement() || !element->IsAllFloating()) continue;
		int el_boxid;
		element->GetInputValue(&el_boxid,BasalforcingsPicoBoxIdEnum);
		if(el_boxid!=boxid) continue;

		Element* basalelement = element->SpawnBasalElement();

		Input* tocs_input=basalelement->GetInput(BasalforcingsPicoSubShelfOceanTempEnum); _assert_(tocs_input); 
		Input* socs_input=basalelement->GetInput(BasalforcingsPicoSubShelfOceanSalinityEnum); _assert_(socs_input);

		basalelement->GetInputValue(&basinid,BasalforcingsPicoBasinIdEnum);
		Gauss* gauss=basalelement->NewGauss(1); gauss->GaussPoint(0);
		tocs_input->GetInputValue(&toc,gauss);
		socs_input->GetInputValue(&soc,gauss);
		delete gauss;
		area=basalelement->GetHorizontalSurfaceArea();
		toc_weighted_avg[basinid]+=toc*area;
		soc_weighted_avg[basinid]+=soc*area;

		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	}

	/*Syncronize across cpus*/
	ISSM_MPI_Allreduce(toc_weighted_avg,toc_sumweightedavg,num_basins,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
	ISSM_MPI_Allreduce(soc_weighted_avg,soc_sumweightedavg,num_basins,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());

	for(int k=0;k<num_basins;k++){
		int p=k*maxbox+boxid; 
		if(boxareas[p]==0) continue;	
		toc_sumweightedavg[k] = toc_sumweightedavg[k]/boxareas[p];
		soc_sumweightedavg[k] = soc_sumweightedavg[k]/boxareas[p];
	}

	femmodel->parameters->AddObject(new DoubleVecParam(BasalforcingsPicoAverageTemperatureEnum,toc_sumweightedavg,num_basins));
	femmodel->parameters->AddObject(new DoubleVecParam(BasalforcingsPicoAverageSalinityEnum,soc_sumweightedavg,num_basins));

	/* Compute overturning weighted avg (box 0 only) */
	if(boxid==0){ 
		overturning_weighted_avg=xNewZeroInit<IssmDouble>(num_basins);
		for(Object* & object : femmodel->elements->objects){
			Element* element=xDynamicCast<Element*>(object);

			/*Check whether we should continue or not*/
			if(!element->IsOnBase()) continue;
			if(!element->IsIceInElement() || !element->IsAllFloating()) continue;
			int el_boxid;
			element->GetInputValue(&el_boxid,BasalforcingsPicoBoxIdEnum);
			if(el_boxid!=boxid) continue;

			Element* basalelement = element->SpawnBasalElement();
	     	Input* overturnings_input=basalelement->GetInput(BasalforcingsPicoSubShelfOceanOverturningEnum); _assert_(overturnings_input);

			basalelement->GetInputValue(&basinid,BasalforcingsPicoBasinIdEnum);
			Gauss* gauss=basalelement->NewGauss(1); gauss->GaussPoint(0);
			overturnings_input->GetInputValue(&overturning,gauss);
			delete gauss;
			area=basalelement->GetHorizontalSurfaceArea();
			overturning_weighted_avg[basinid]+=overturning*area;
			basalelement->FindParam(&domaintype,DomainTypeEnum);

			if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
		}

		/*Syncronize across cpus*/
		ISSM_MPI_Allreduce(overturning_weighted_avg,overturning_sumweightedavg,num_basins,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());

		for(int k=0;k<num_basins;k++){
			int p=k*maxbox+boxid;
			if(boxareas[p]==0.) continue;
			overturning_sumweightedavg[k] = overturning_sumweightedavg[k]/boxareas[p];
		}
		femmodel->parameters->AddObject(new DoubleVecParam(BasalforcingsPicoAverageOverturningEnum,overturning_sumweightedavg,num_basins));
	}

	/*Cleanup and return */
	xDelete<IssmDouble>(overturning_sumweightedavg);
	xDelete<IssmDouble>(toc_sumweightedavg);
	xDelete<IssmDouble>(soc_sumweightedavg);
	xDelete<IssmDouble>(overturning_weighted_avg);
	xDelete<IssmDouble>(toc_weighted_avg);
	xDelete<IssmDouble>(soc_weighted_avg);
	xDelete<IssmDouble>(boxareas);
}/*}}}*/
void ComputeBasalMeltPlume(FemModel* femmodel){/*{{{*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->PicoComputeBasalMelt();
	}
}/*}}}*/
