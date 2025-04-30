#include "./StressbalanceAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../solutionsequences/solutionsequences.h"
#include "../classes/IoModel.h"
#include "../classes/FemModel.h"
#include "../classes/Constraints/Constraints.h"
#include "../classes/Constraints/Constraint.h"
#include "../classes/Constraints/SpcStatic.h"
#include "../classes/Params/Parameters.h"
#include "../classes/Nodes.h"
#include "../classes/Node.h"
#include "../classes/Elements/Elements.h"
#include "../classes/Elements/Element.h"
#include "../modules/ModelProcessorx/ModelProcessorx.h"
#include "../modules/IoModelToConstraintsx/IoModelToConstraintsx.h"
#include "../modules/InputUpdateFromConstantx/InputUpdateFromConstantx.h"
#include "../modules/SetActiveNodesLSMx/SetActiveNodesLSMx.h"
#include "../cores/cores.h"

//#define FSANALYTICAL 10
//#define LATERALFRICTION 1
//#define DISCSLOPE 1 //testing for SSA

/*Model processing*/
void StressbalanceAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Intermediary*/
	int        i,j;
	int        finiteelement,p_fe,v_fe;
	IssmDouble g;
	IssmDouble rho_ice;
	IssmDouble FSreconditioning;
	bool       isSIA,isSSA,isL1L2,isMOLHO,isHO,isFS,iscoupling;
	bool       spcpresent = false;
	int        Mx,Nx;
	int        My,Ny;
	int        Mz,Nz;
	IssmDouble *spcvx          = NULL;
	IssmDouble *spcvy          = NULL;
	IssmDouble *spcvz          = NULL;
	IssmDouble *nodeonSSA = NULL;
	IssmDouble *nodeonHO   = NULL;
	IssmDouble *nodeonFS   = NULL;
	IssmDouble *nodeonbase      = NULL;
	IssmDouble *groundedice_ls = NULL;
	IssmDouble *vertices_type  = NULL;
	IssmDouble *surface        = NULL;
	IssmDouble *z              = NULL;
	IssmDouble *timesx=NULL;
	IssmDouble *timesy=NULL;
	IssmDouble *timesz=NULL;
   IssmDouble* values=NULL;

	/*Fetch parameters: */
	iomodel->FindConstant(&g,"md.constants.g");
	iomodel->FindConstant(&rho_ice,"md.materials.rho_ice");
	iomodel->FindConstant(&FSreconditioning,"md.stressbalance.FSreconditioning");
	iomodel->FindConstant(&isSIA,"md.flowequation.isSIA");
	iomodel->FindConstant(&isSSA,"md.flowequation.isSSA");
	iomodel->FindConstant(&isL1L2,"md.flowequation.isL1L2");
	iomodel->FindConstant(&isMOLHO,"md.flowequation.isMOLHO");
	iomodel->FindConstant(&isHO,"md.flowequation.isHO");
	iomodel->FindConstant(&isFS,"md.flowequation.isFS");

	/*Is this model only SIA??*/
	if(!isSSA && !isHO && !isFS && !isL1L2 && !isMOLHO) return;

	/*Do we have coupling*/
	if((isSIA?1.:0.) + (isSSA?1.:0.) + (isL1L2?1.:0.) + (isMOLHO?1.:0.) + (isHO?1.:0.) + (isFS?1.:0.) >1.)
	 iscoupling = true;
	else
	 iscoupling = false;

	/*If no coupling, call Regular IoModelToConstraintsx, else, use P1 elements only*/
	if(!iscoupling){

		/*Get finite element type*/
		if(isSSA)       iomodel->FindConstant(&finiteelement,"md.flowequation.fe_SSA");
		else if(isL1L2) finiteelement = P1Enum;
		else if(isMOLHO) finiteelement = P1Enum;
		else if(isHO)   iomodel->FindConstant(&finiteelement,"md.flowequation.fe_HO");
		else if(isFS){  iomodel->FindConstant(&finiteelement,"md.flowequation.fe_FS");
			/*Deduce velocity interpolation from finite element*/
			switch(finiteelement){
				case P1P1Enum              : v_fe = P1Enum;       p_fe = P1Enum;   break;
				case P1P1GLSEnum           : v_fe = P1Enum;       p_fe = P1Enum;   break;
				case MINIcondensedEnum     : v_fe = P1bubbleEnum; p_fe = P1Enum;   break;
				case MINIEnum              : v_fe = P1bubbleEnum; p_fe = P1Enum;   break;
				case TaylorHoodEnum        : v_fe = P2Enum;       p_fe = P1Enum;   break;
				case XTaylorHoodEnum       : v_fe = P2Enum;       p_fe = P1Enum;   break;
				case LATaylorHoodEnum      : v_fe = P2Enum;       p_fe = NoneEnum; break;
				case LACrouzeixRaviartEnum : v_fe = P2bubbleEnum; p_fe = NoneEnum; break;
				case OneLayerP4zEnum       : v_fe = P2xP4Enum;    p_fe = P1Enum;   break;
				case CrouzeixRaviartEnum   : v_fe = P2bubbleEnum; p_fe = P1DGEnum; break;
				default: _error_("finite element "<<EnumToStringx(finiteelement)<<" not supported");
			}
		}
		else{
			_error_("model not supported yet");
		}

		if(isFS){

			/*Constraint at the bedrock interface (v.n = vz = 0) (Coordinates will be updated according to the bed slope)*/
			iomodel->FetchData(&vertices_type,NULL,NULL,"md.flowequation.vertex_equation");
			iomodel->FetchData(&nodeonFS,NULL,NULL,"md.flowequation.borderFS");
			iomodel->FetchData(&nodeonbase,NULL,NULL,"md.mesh.vertexonbase");
			iomodel->FetchData(&groundedice_ls,NULL,NULL,"md.mask.ocean_levelset");
			if(iomodel->domaintype==Domain3DEnum){
				iomodel->FetchData(&spcvz,&Mz,&Nz,"md.stressbalance.spcvz");
			}
			else if (iomodel->domaintype==Domain2DverticalEnum){
				iomodel->FetchData(&spcvz,&Mz,&Nz,"md.stressbalance.spcvy");
			}
			else{
				_error_("not supported yet");
			}
			if(iomodel->domaintype==Domain3DEnum){
				IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvx",StressbalanceAnalysisEnum,v_fe,0);
				IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvy",StressbalanceAnalysisEnum,v_fe,1);
				IoModelToConstraintsx(constraints,iomodel,spcvz,Mz,Nz,StressbalanceAnalysisEnum,v_fe,2);
				iomodel->DeleteData(spcvz,"md.stressbalance.spcvz");
			}
			else if (iomodel->domaintype==Domain2DverticalEnum){
				IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvx",StressbalanceAnalysisEnum,v_fe,0);
				IoModelToConstraintsx(constraints,iomodel,spcvz,Mz,Nz,StressbalanceAnalysisEnum,v_fe,1);
				iomodel->DeleteData(spcvz,"md.stressbalance.spcvy");
			}
			else{
				_error_("not supported yet");
			}
			iomodel->DeleteData(vertices_type,"md.flowequation.vertex_equation");
			iomodel->DeleteData(nodeonFS,"md.flowequation.borderFS");
			iomodel->DeleteData(nodeonbase,"md.mesh.vertexonbase");
			iomodel->DeleteData(groundedice_ls,"md.mask.ocean_levelset");

			/*Pressure spc*/
			int count = constraints->Size();
			iomodel->FetchData(&vertices_type,NULL,NULL,"md.flowequation.vertex_equation");
			iomodel->FetchData(&surface,NULL,NULL,"md.geometry.surface");
			iomodel->FetchData(&z,NULL,NULL,"md.mesh.z");
			bool addpressurespc = false;
			for(i=0;i<iomodel->numberofvertices;i++){
				if(IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==NoneApproximationEnum && iomodel->my_vertices[i]){
					addpressurespc = true;
					break;
				}
			}
			if(addpressurespc){
				switch(finiteelement){
					case P1P1Enum: case P1P1GLSEnum:
						for(i=0;i<iomodel->numberofvertices;i++){
							if(iomodel->my_vertices[i]){
								if(IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==NoneApproximationEnum){
									constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+i+1,0,g*rho_ice*(surface[i]-z[i])/FSreconditioning,StressbalanceAnalysisEnum));
									count++;
								}
							}
						}
						break;
					case MINIEnum: case MINIcondensedEnum:
						for(i=0;i<iomodel->numberofvertices;i++){
							if(iomodel->my_vertices[i]){
								if(IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==NoneApproximationEnum){
									constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+iomodel->numberofelements+i+1,0,g*rho_ice*(surface[i]-z[i])/FSreconditioning,StressbalanceAnalysisEnum));
									count++;
								}
							}
						}
						break;
					case TaylorHoodEnum:
						for(i=0;i<iomodel->numberofvertices;i++){
							if(iomodel->my_vertices[i]){
								if(IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==NoneApproximationEnum){
									constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+i+1,0,g*rho_ice*(surface[i]-z[i])/FSreconditioning,StressbalanceAnalysisEnum));
									count++;
								}
							}
						}
						break;
					default:
						_error_("not implemented yet");
				}
			}
			iomodel->DeleteData(vertices_type,"md.flowequation.vertex_equation");
			iomodel->DeleteData(surface,"md.geometry.surface");
			iomodel->DeleteData(z,"md.mesh.z");
		}
		else{
			if(!isMOLHO){
				IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvx",StressbalanceAnalysisEnum,finiteelement,0);
				if(iomodel->domaintype!=Domain2DverticalEnum){
					IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvy",StressbalanceAnalysisEnum,finiteelement,1);
				}
			}
			else{//MOLHO 
				IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvx_base",StressbalanceAnalysisEnum,finiteelement,0);
				IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvx_shear",StressbalanceAnalysisEnum,finiteelement,1);
				IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvy_base",StressbalanceAnalysisEnum,finiteelement,2);
				IoModelToConstraintsx(constraints,iomodel,"md.stressbalance.spcvy_shear",StressbalanceAnalysisEnum,finiteelement,3);
			}
		}

		return;
	}

	/*Constraints: fetch data: */
	iomodel->FetchData(&spcvx,&Mx,&Nx,"md.stressbalance.spcvx");
	iomodel->FetchData(&spcvy,&My,&Ny,"md.stressbalance.spcvy");
	iomodel->FetchData(&spcvz,&Mz,&Nz,"md.stressbalance.spcvz");
	iomodel->FetchData(&nodeonSSA,NULL,NULL,"md.flowequation.borderSSA");
	if(iomodel->domaintype==Domain3DEnum)iomodel->FetchData(&nodeonHO,NULL,NULL,"md.flowequation.borderHO");
	if(iomodel->domaintype==Domain3DEnum)iomodel->FetchData(&nodeonFS,NULL,NULL,"md.flowequation.borderFS");
	if(iomodel->domaintype==Domain3DEnum)iomodel->FetchData(&nodeonbase,NULL,NULL,"md.mesh.vertexonbase");
	if(iomodel->domaintype==Domain3DEnum)iomodel->FetchData(&groundedice_ls,NULL,NULL,"md.mask.ocean_levelset");
	iomodel->FetchData(&vertices_type,NULL,NULL,"md.flowequation.vertex_equation");
	iomodel->FetchData(&surface,NULL,NULL,"md.geometry.surface");
	iomodel->FetchData(&z,NULL,NULL,"md.mesh.z");

	/*Initialize counter: */
	int count=0;

	/*figure out times: */
	timesx=xNew<IssmDouble>(Nx);
	for(j=0;j<Nx;j++){
		timesx[j]=spcvx[(Mx-1)*Nx+j];
	}
	/*figure out times: */
	timesy=xNew<IssmDouble>(Ny);
	for(j=0;j<Ny;j++){
		timesy[j]=spcvy[(My-1)*Ny+j];
	}
	/*figure out times: */
	timesz=xNew<IssmDouble>(Nz);
	for(j=0;j<Nz;j++){
		timesz[j]=spcvz[(Mz-1)*Nz+j];
	}

	/*Create spcs from x,y,z, as well as the spc values on those spcs: */
	for(i=0;i<iomodel->numberofvertices;i++){
		if(iomodel->my_vertices[i]){

			/*Start with adding spcs of coupling: zero at the border SSA/HO for the appropriate dofs*/
			if(IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==SSAHOApproximationEnum){
				/*If grionSSA, spc HO dofs: 3 & 4*/
					if (reCast<int,IssmDouble>(nodeonHO[i])){
						constraints->AddObject(new SpcStatic(count+1,i+1,0,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						constraints->AddObject(new SpcStatic(count+1,i+1,1,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						if (!xIsNan<IssmDouble>(spcvx[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,2,spcvx[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}
						if (!xIsNan<IssmDouble>(spcvy[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,3,spcvy[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}

					}
					else if (reCast<int,IssmDouble>(nodeonSSA[i])){
						constraints->AddObject(new SpcStatic(count+1,i+1,2,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						constraints->AddObject(new SpcStatic(count+1,i+1,3,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						if (!xIsNan<IssmDouble>(spcvx[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,0,spcvx[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}
						if (!xIsNan<IssmDouble>(spcvy[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,1,spcvy[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}

					}
					else _error_("if vertices_type is SSAHO, you shoud have nodeonHO or nodeonSSA");
			}
			/*Also add spcs of coupling: zero at the border HO/FS for the appropriate dofs*/
			else if (IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==HOFSApproximationEnum){
				/*If grion,HO spc FS dofs: 3 4 & 5*/
					if (reCast<int,IssmDouble>(nodeonHO[i])){
						constraints->AddObject(new SpcStatic(count+1,i+1,2,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						constraints->AddObject(new SpcStatic(count+1,i+1,3,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						constraints->AddObject(new SpcStatic(count+1,i+1,4,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						if (!xIsNan<IssmDouble>(spcvx[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,0,spcvx[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}
						if (!xIsNan<IssmDouble>(spcvy[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,1,spcvy[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}

					}
					else if (reCast<int,IssmDouble>(nodeonFS[i])){ //spc HO nodes: 1 & 2
						constraints->AddObject(new SpcStatic(count+1,i+1,0,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						constraints->AddObject(new SpcStatic(count+1,i+1,1,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						if (!xIsNan<IssmDouble>(spcvx[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,2,spcvx[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}
						if (!xIsNan<IssmDouble>(spcvy[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,3,spcvy[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}
						if (!xIsNan<IssmDouble>(spcvz[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,4,spcvz[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}
					}
					else _error_("if vertices_type is HOFS, you shoud have nodeonHO or nodeonFS");
			}
			/*Also add spcs of coupling: zero at the border HO/FS for the appropriate dofs*/
			else if (IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==SSAFSApproximationEnum){
				/*If grion,HO spc FS dofs: 3 4 & 5*/
					if (reCast<int,IssmDouble>(nodeonSSA[i])){
						constraints->AddObject(new SpcStatic(count+1,i+1,2,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						constraints->AddObject(new SpcStatic(count+1,i+1,3,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						constraints->AddObject(new SpcStatic(count+1,i+1,4,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						if (!xIsNan<IssmDouble>(spcvx[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,0,spcvx[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}
						if (!xIsNan<IssmDouble>(spcvy[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,1,spcvy[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}

					}
					else if (reCast<int,IssmDouble>(nodeonFS[i])){ //spc SSA nodes: 1 & 2
						constraints->AddObject(new SpcStatic(count+1,i+1,0,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						constraints->AddObject(new SpcStatic(count+1,i+1,1,0,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
						count++;
						if (!xIsNan<IssmDouble>(spcvx[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,2,spcvx[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}
						if (!xIsNan<IssmDouble>(spcvy[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,3,spcvy[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}
						if (!xIsNan<IssmDouble>(spcvz[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,4,spcvz[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
							count++;
						}
					}
					else _error_("if vertices_type is SSAFS, you shoud have nodeonSSA or nodeonFS");
			}
			/*Now add the regular spcs*/
			else{
				if (Mx==iomodel->numberofvertices && !xIsNan<IssmDouble>(spcvx[i])){
					constraints->AddObject(new SpcStatic(count+1,i+1,0,spcvx[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vx.
					count++;

				}
				else if (Mx==iomodel->numberofvertices+1) {
					/*figure out times and values: */
					values=xNew<IssmDouble>(Nx);
					spcpresent=false;
					for(j=0;j<Nx;j++){
						values[j]=spcvx[i*Nx+j];
						if(!xIsNan<IssmDouble>(values[j]))spcpresent=true; //NaN means no spc by default
					}

					if(spcpresent){
						constraints->AddObject(new SpcTransient(count+1,i+1,0,Nx,timesx,values,StressbalanceAnalysisEnum));
						count++;
					}
					xDelete<IssmDouble>(values);
				}
				else if (IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==SIAApproximationEnum){
					constraints->AddObject(new SpcDynamic(count+1,i+1,0,0.,StressbalanceAnalysisEnum));
					count++;
				}

				if (My==iomodel->numberofvertices && !xIsNan<IssmDouble>(spcvy[i])){
					constraints->AddObject(new SpcStatic(count+1,i+1,1,spcvy[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 1 to vy.
					count++;
				}
				else if (My==iomodel->numberofvertices+1){
					/*figure out times and values: */
					values=xNew<IssmDouble>(Ny);
					spcpresent=false;
					for(j=0;j<Ny;j++){
						values[j]=spcvy[i*Ny+j];
						if(!xIsNan<IssmDouble>(values[j]))spcpresent=true; //NaN means no spc by default
					}
					if(spcpresent){
						constraints->AddObject(new SpcTransient(count+1,i+1,1,Ny,timesy,values,StressbalanceAnalysisEnum));
						count++;
					}
					xDelete<IssmDouble>(values);
				}
				else if (IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==SIAApproximationEnum){
					constraints->AddObject(new SpcDynamic(count+1,i+1,1,0.,StressbalanceAnalysisEnum));
					count++;
				}

				if (IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==FSApproximationEnum ||  (IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==NoneApproximationEnum)){
					if (Mz==iomodel->numberofvertices && !xIsNan<IssmDouble>(spcvz[i])){
						constraints->AddObject(new SpcStatic(count+1,i+1,2,spcvz[i],StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 2 to vy
						count++;
					}
					else if (Mz==iomodel->numberofvertices+1){
						/*figure out times and values: */
						values=xNew<IssmDouble>(Nz);
						spcpresent=false;
						for(j=0;j<Nz;j++){
							values[j]=spcvz[i*Nz+j];
							if(!xIsNan<IssmDouble>(values[j]))spcpresent=true; //NaN means no spc by default
						}
						if(spcpresent){
							constraints->AddObject(new SpcTransient(count+1,i+1,2,Nz,timesz,values,StressbalanceAnalysisEnum));
							count++;
						}
						xDelete<IssmDouble>(values);
					}

				}
				if (IoCodeToEnumVertexEquation(reCast<int>(vertices_type[i]))==NoneApproximationEnum){
					constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+i+1,0,g*rho_ice*(surface[i]-z[i])/FSreconditioning,StressbalanceAnalysisEnum)); //add count'th spc, on node i+1, setting dof 2 to vy
					count++;
				}
			}
		}
	}

	/*Free data: */
	iomodel->DeleteData(spcvx,"md.stressbalance.spcvx");
	iomodel->DeleteData(spcvy,"md.stressbalance.spcvy");
	iomodel->DeleteData(spcvz,"md.stressbalance.spcvz");
	iomodel->DeleteData(nodeonSSA,"md.flowequation.borderSSA");
	if(iomodel->domaintype==Domain3DEnum)iomodel->DeleteData(nodeonHO,"md.flowequation.borderHO");
	if(iomodel->domaintype==Domain3DEnum)iomodel->DeleteData(nodeonFS,"md.flowequation.borderFS");
	if(iomodel->domaintype==Domain3DEnum)iomodel->DeleteData(nodeonbase,"md.mesh.vertexonbase");
	if(iomodel->domaintype==Domain3DEnum)iomodel->DeleteData(groundedice_ls,"md.mask.ocean_levelset");
	iomodel->DeleteData(vertices_type,"md.flowequation.vertex_equation");
	iomodel->DeleteData(surface,"md.geometry.surface");
	iomodel->DeleteData(z,"md.mesh.z");

	/*Free resources:*/
	xDelete<IssmDouble>(timesx);
	xDelete<IssmDouble>(timesy);
	xDelete<IssmDouble>(timesz);
	xDelete<IssmDouble>(values);

}/*}}}*/
void StressbalanceAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Intermediary*/
	const int   RIFTINFOSIZE = 12;
	int         i;
	int         count;
	int         penpair_ids[2];
	bool        isSSA,isL1L2,isMOLHO,isHO,isFS;
	int         numpenalties,numrifts,numriftsegments;
	IssmDouble *riftinfo       = NULL;
	IssmDouble *penalties      = NULL;
	int         assert_int;

	/*Fetch parameters: */
	iomodel->FindConstant(&isL1L2,"md.flowequation.isL1L2");
	iomodel->FindConstant(&isMOLHO,"md.flowequation.isMOLHO");
	iomodel->FindConstant(&isFS,"md.flowequation.isFS");
	iomodel->FindConstant(&isSSA,"md.flowequation.isSSA");
	iomodel->FindConstant(&isHO,"md.flowequation.isHO");
	iomodel->FindConstant(&numrifts,"md.rifts.numrifts");

	/*Is this SIA only?*/
	if(!isSSA && !isHO && !isFS && !isL1L2 && !isMOLHO) return;

	/*Initialize counter: */
	count=0;

	/*Create Penpair for penalties: */
	iomodel->FetchData(&penalties,&numpenalties,NULL,"md.stressbalance.vertex_pairing");

	for(i=0;i<numpenalties;i++){

		if(iomodel->my_vertices[reCast<int,IssmDouble>(penalties[2*i+0]-1)]){

			/*In debugging mode, check that the second node is in the same cpu*/
			assert_int=iomodel->my_vertices[reCast<int,IssmDouble>(penalties[2*i+1]-1)]; _assert_(assert_int);

			/*Get node ids*/
			penpair_ids[0]=reCast<int,IssmDouble>(penalties[2*i+0]);
			penpair_ids[1]=reCast<int,IssmDouble>(penalties[2*i+1]);

			/*Create Load*/
			loads->AddObject(new Penpair(count+1,&penpair_ids[0]));
			count++;
		}
	}

	/*Free resources: */
	iomodel->DeleteData(penalties,"md.stressbalance.vertex_pairing");

	/*Create Riffront loads for rifts: */
	if(numrifts){
		iomodel->FetchData(&riftinfo,&numriftsegments,NULL,"md.rifts.riftstruct");
		iomodel->FetchData(5,"md.rifts.riftstruct","md.geometry.thickness","md.geometry.base","md.geometry.surface","md.mask.ocean_levelset");
		for(i=0;i<numriftsegments;i++){
			if(iomodel->my_elements[reCast<int,IssmDouble>(*(riftinfo+RIFTINFOSIZE*i+2))-1]){
				loads->AddObject(new Riftfront(count+1,i,iomodel));
				count++;
			}
		}
		iomodel->DeleteData(5,"md.rifts.riftstruct","md.geometry.thickness","md.geometry.base","md.geometry.surface","md.mask.ocean_levelset");
		xDelete<IssmDouble>(riftinfo);
	}
}/*}}}*/
void StressbalanceAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Intermediary*/
	bool isSSA,isL1L2,isMOLHO,isHO,isFS,iscoupling;
	int  finiteelement=-1,approximation=-1;

	/*Fetch parameters: */
	iomodel->FindConstant(&isSSA,"md.flowequation.isSSA");
	iomodel->FindConstant(&isL1L2,"md.flowequation.isL1L2");
	iomodel->FindConstant(&isMOLHO,"md.flowequation.isMOLHO");
	iomodel->FindConstant(&isHO,"md.flowequation.isHO");
	iomodel->FindConstant(&isFS,"md.flowequation.isFS");

	/*Now, check that we have non SIA elements */
	if(!isSSA && !isL1L2 && !isMOLHO && !isHO && !isFS) return;

	/*Do we have coupling*/
	if( (isSSA?1.:0.) + (isL1L2?1.:0.) + (isMOLHO?1.:0.) + (isHO?1.:0.) + (isFS?1.:0.) >1.)
	 iscoupling = true;
	else
	 iscoupling = false;

	/*If no coupling, call Regular CreateNodes, else, use P1 elements only*/
	if(!iscoupling){

		/*Get finite element type*/
		if(isSSA){
			approximation=SSAApproximationEnum;
			iomodel->FindConstant(&finiteelement,"md.flowequation.fe_SSA");
		}
		else if(isL1L2){
			approximation = L1L2ApproximationEnum;
			finiteelement = P1Enum;
		}
		else if(isMOLHO){
			approximation = MOLHOApproximationEnum;
			finiteelement = P1Enum;
		}
		else if(isHO){
			approximation = HOApproximationEnum;
			iomodel->FindConstant(&finiteelement,"md.flowequation.fe_HO");
		}
		else if(isFS){
			approximation = FSApproximationEnum;
			iomodel->FindConstant(&finiteelement,"md.flowequation.fe_FS");
		}
		if(!isamr){
			iomodel->FetchData(3,"md.flowequation.borderSSA","md.flowequation.vertex_equation","md.stressbalance.referential");
			if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(3,"md.mesh.vertexonbase","md.mesh.vertexonsurface","md.flowequation.borderFS");
		}
		::CreateNodes(nodes,iomodel,StressbalanceAnalysisEnum,finiteelement,isamr,approximation);
		if(!isamr){
			iomodel->DeleteData(6,"md.mesh.vertexonbase","md.mesh.vertexonsurface","md.flowequation.borderSSA","md.flowequation.vertex_equation",
					"md.stressbalance.referential","md.flowequation.borderFS");
		}
	}
	else{
		/*Coupling: we are going to create P1 Elements only*/
		iomodel->FetchData(6,"md.mesh.vertexonbase","md.mesh.vertexonsurface","md.flowequation.borderSSA","md.flowequation.borderFS",
					"md.flowequation.vertex_equation","md.stressbalance.referential");
		if(isFS){
			int* approximations = xNew<int>(2*iomodel->numberofvertices+iomodel->numberofelements);
			for(int i=0;i<iomodel->numberofvertices;i++){
				approximation = IoCodeToEnumVertexEquation(reCast<int>(iomodel->Data("md.flowequation.vertex_equation")[i]));
				if(approximation==FSApproximationEnum)  approximation=FSvelocityEnum;
				approximations[i] = approximation;
			}
			for(int i=0;i<iomodel->numberofelements;i++) approximations[iomodel->numberofvertices+i] = FSvelocityEnum;
			for(int i=0;i<iomodel->numberofvertices;i++) approximations[iomodel->numberofvertices+iomodel->numberofelements+i] = FSpressureEnum;
			::CreateNodes(nodes,iomodel,StressbalanceAnalysisEnum,MINIcondensedEnum,isamr,0,approximations);
			xDelete<int>(approximations);

			for(Object* & object: nodes->objects){
				Node* node=xDynamicCast<Node*>(object);
				int   sid = node->Sid();
				if(sid>=iomodel->numberofvertices+iomodel->numberofelements){
					/*Constrain pressure if not FS*/
					int id = sid - (iomodel->numberofvertices+iomodel->numberofelements);
					approximation=IoCodeToEnumVertexEquation(reCast<int>(iomodel->Data("md.flowequation.vertex_equation")[id]));
					if(approximation==HOApproximationEnum || approximation==SSAApproximationEnum){
						node->Deactivate();
					}
				}
			}
		}
		else{
			int* approximations = xNew<int>(iomodel->numberofvertices);
			for(int i=0;i<iomodel->numberofvertices;i++) approximations[i] = IoCodeToEnumVertexEquation(reCast<int>(iomodel->Data("md.flowequation.vertex_equation")[i]));
			::CreateNodes(nodes,iomodel,StressbalanceAnalysisEnum,P1Enum,isamr,0,approximations);
			xDelete<int>(approximations);
		}
		iomodel->DeleteData(6,"md.mesh.vertexonbase","md.mesh.vertexonsurface","md.flowequation.borderSSA","md.flowequation.borderFS",
					"md.flowequation.vertex_equation","md.stressbalance.referential");
	}
}/*}}}*/
int  StressbalanceAnalysis::DofsPerNode(int** pdoftype,int domaintype,int approximation){/*{{{*/

	/*output*/
	int *doftype = NULL;
	int  numdofs;

	switch(approximation){
		case SSAApproximationEnum:
			 switch(domaintype){
				 case Domain3DEnum:           numdofs=2; break;
				 case Domain2DhorizontalEnum: numdofs=2; break;
				 case Domain2DverticalEnum:   numdofs=1; break;
				 default: _error_("mesh type not supported yet");
			 }
			 break;
		case L1L2ApproximationEnum: numdofs = 2; break;
		case MOLHOApproximationEnum: numdofs = 4; break;
		case HOApproximationEnum:
			 switch(domaintype){
				 case Domain3DEnum:         numdofs=2; break;
				 case Domain2DverticalEnum: numdofs=1; break;
				 default: _error_("mesh type not supported yet");
			 }
			 break;
		case SIAApproximationEnum:  numdofs =2; break;
		case FSvelocityEnum:
			 switch(domaintype){
				 case Domain3DEnum:         numdofs=3; break;
				 case Domain2DverticalEnum: numdofs=2; break;
				 default: _error_("mesh type not supported yet");
			}
			break;
		case FSpressureEnum: numdofs=1; break;
		case NoneApproximationEnum:
			 switch(domaintype){
				 case Domain3DEnum:         numdofs=4; break;
				 case Domain2DverticalEnum: numdofs=3; break;
				 default: _error_("mesh type not supported yet");
			}
			break;
		case SSAHOApproximationEnum:
			numdofs=4;
			doftype=xNew<int>(numdofs);
			doftype[0]=SSAApproximationEnum;
			doftype[1]=SSAApproximationEnum;
			doftype[2]=HOApproximationEnum;
			doftype[3]=HOApproximationEnum;
			break;
		case HOFSApproximationEnum:
			numdofs=5;
			doftype=xNew<int>(numdofs);
			doftype[0]=HOApproximationEnum;
			doftype[1]=HOApproximationEnum;
			doftype[2]=FSvelocityEnum;
			doftype[3]=FSvelocityEnum;
			doftype[4]=FSvelocityEnum;
			break;
		case SSAFSApproximationEnum:
			numdofs=5;
			doftype=xNew<int>(numdofs);
			doftype[0]=SSAApproximationEnum;
			doftype[1]=SSAApproximationEnum;
			doftype[2]=FSvelocityEnum;
			doftype[3]=FSvelocityEnum;
			doftype[4]=FSvelocityEnum;
			break;
		default:
			_error_("Approximation " << EnumToStringx(approximation) << " not implemented yet");
	}

	/*Assign output pointer and return*/
	*pdoftype = doftype;
	return numdofs;
}/*}}}*/
void StressbalanceAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Intermediaries*/
	int    materials_type,finiteelement,fe_FS;
	int    approximation;
	int*   finiteelement_list=NULL;
	bool   isSSA,isL1L2,isMOLHO,isHO,isFS,iscoupling,ishydrologylayer;
	bool   control_analysis;
	bool   dakota_analysis;
	bool   ismovingfront;

	/*Fetch constants needed: */
	iomodel->FindConstant(&isSSA,"md.flowequation.isSSA");
	iomodel->FindConstant(&isL1L2,"md.flowequation.isL1L2");
	iomodel->FindConstant(&isMOLHO,"md.flowequation.isMOLHO");
	iomodel->FindConstant(&isHO,"md.flowequation.isHO");
	iomodel->FindConstant(&isFS,"md.flowequation.isFS");
	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");
	iomodel->FindConstant(&materials_type,"md.materials.type");
	iomodel->FindConstant(&ismovingfront,"md.transient.ismovingfront");
	iomodel->FindConstant(&ishydrologylayer,"md.stressbalance.ishydrologylayer");

	/*return if no processing required*/
	if(!isSSA && !isL1L2 && !isMOLHO && !isHO && !isFS) return;

	/*Fetch data needed and allocate vectors: */
	iomodel->FetchData(1,"md.flowequation.element_equation");
	finiteelement_list=xNewZeroInit<int>(iomodel->numberofelements);

	/*Do we have coupling*/
	if( (isSSA?1.:0.) + (isL1L2?1.:0.) + (isMOLHO?1.:0.) + (isHO?1.:0.) + (isFS?1.:0.) >1.)
	 iscoupling = true;
	else
	 iscoupling = false;

	/*Get finite element type*/
	if(!iscoupling){
		if(isSSA)       iomodel->FindConstant(&finiteelement,"md.flowequation.fe_SSA");
		else if(isL1L2) finiteelement = P1Enum;
		else if(isMOLHO) finiteelement = P1Enum;
		else if(isHO)   iomodel->FindConstant(&finiteelement,"md.flowequation.fe_HO");
		else if(isFS)   iomodel->FindConstant(&finiteelement,"md.flowequation.fe_FS");
		for(int i=0;i<iomodel->numberofelements;i++){
			finiteelement_list[i]=finiteelement;
		}
	}
	else{
		if(isFS){
			for(int i=0;i<iomodel->numberofelements;i++){
				approximation=IoCodeToEnumElementEquation(reCast<int>(iomodel->Data("md.flowequation.element_equation")[i]));
				if(approximation==FSApproximationEnum || approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
					finiteelement_list[i]=MINIcondensedEnum;
				}
				else{
					finiteelement_list[i]=P1Enum;
				}
			}
		}
		else{
			finiteelement = P1Enum;
			for(int i=0;i<iomodel->numberofelements;i++){
				finiteelement_list[i]=finiteelement;
			}
		}
	}

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement_list[i]);

			/*Need to know the type of approximation for this element*/
			if(iomodel->Data("md.flowequation.element_equation")){
				inputs->SetInput(ApproximationEnum,counter,IoCodeToEnumElementEquation(reCast<int>(iomodel->Data("md.flowequation.element_equation")[i])));
			}

			counter++;
		}
	}

	/*Create inputs: */
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum,0.);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum,0.);
	/*Hydrology layer*/
	if(ishydrologylayer){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.watercolumn",HydrologySheetThicknessEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_height",HydrologyBumpHeightEnum);
	}
	/*MOLHO*/
	if(isMOLHO){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxShearEnum,0.);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyShearEnum,0.);
		/*3D MOLHO also need to have VxBase and VyBase for reconstruting Vx and Vy*/
		if (iomodel->domaintype==Domain3DEnum) {
			iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxBaseEnum,0.);
			iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyBaseEnum,0.);
		}
	}
   if(iomodel->domaintype==Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxBaseEnum,0.);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyBaseEnum,0.);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxSurfaceEnum,0.);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VySurfaceEnum,0.);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.stressbalance.loadingforcex",LoadingforceXEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.stressbalance.loadingforcey",LoadingforceYEnum);
	#ifdef LATERALFRICTION
	iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonboundary",MeshVertexonboundaryEnum);
	#endif

	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	if(iomodel->domaintype==Domain3DEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.flowequation.borderFS",FlowequationBorderFSEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.stressbalance.loadingforcez",LoadingforceZEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vz",VzEnum,0.);
	}
	if(isFS){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.pressure",PressureEnum,0.);

		/*Add basal forcings to compute melt rate*/
		bool isstochastic;
	   int basalforcing_model;
	   int melt_parameterization;
	   iomodel->FindConstant(&basalforcing_model,"md.basalforcings.model");
	   iomodel->FindConstant(&isstochastic,"md.stochasticforcing.isstochasticforcing");
	   iomodel->FindConstant(&melt_parameterization,"md.frontalforcings.parameterization");
		switch(basalforcing_model){
			case FloatingMeltRateEnum:
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.floatingice_melting_rate",BasalforcingsFloatingiceMeltingRateEnum);
				if(isstochastic){
					iomodel->FetchDataToInput(inputs,elements,"md.stochasticforcing.default_id",StochasticForcingDefaultIdEnum);
					iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.floatingice_melting_rate",BaselineBasalforcingsFloatingiceMeltingRateEnum);
				}
				break;
			case LinearFloatingMeltRateEnum:
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.perturbation_melting_rate",BasalforcingsPerturbationMeltingRateEnum,0.);
				break;
			case MismipFloatingMeltRateEnum:
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.meltrate_factor",BasalforcingsMeltrateFactorEnum);
				break;
			case MantlePlumeGeothermalFluxEnum:
				break;
			case SpatialLinearFloatingMeltRateEnum:
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.deepwater_melting_rate",BasalforcingsSpatialDeepwaterMeltingRateEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.deepwater_elevation",BasalforcingsSpatialDeepwaterElevationEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.upperwater_melting_rate",BasalforcingsSpatialUpperwaterMeltingRateEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.upperwater_elevation",BasalforcingsSpatialUpperwaterElevationEnum);
				if(isstochastic){
					iomodel->FetchDataToInput(inputs,elements,"md.stochasticforcing.default_id",StochasticForcingDefaultIdEnum);
					iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.deepwater_melting_rate",BaselineBasalforcingsSpatialDeepwaterMeltingRateEnum);
				}
				break;
			case BasalforcingsPicoEnum:
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.basin_id",BasalforcingsPicoBasinIdEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.overturning_coeff",BasalforcingsPicoOverturningCoeffEnum);
				break;
			case BasalforcingsIsmip6Enum:
				iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.basin_id",BasalforcingsIsmip6BasinIdEnum);
				break;
			case BeckmannGoosseFloatingMeltRateEnum:
				bool isthermalforcing;
				iomodel->FindConstant(&isthermalforcing,"md.basalforcings.isthermalforcing");
      	   iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.meltrate_factor",BasalforcingsMeltrateFactorEnum);
      	   if(isthermalforcing==0){
      	      iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.ocean_salinity",BasalforcingsOceanSalinityEnum);
      	      iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.ocean_temp",BasalforcingsOceanTempEnum);
      	   }
      	   else if(melt_parameterization!=FrontalForcingsRignotarmaEnum){
      	      iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.ocean_thermalforcing",ThermalForcingEnum);
      	   }
				break;
			default:
				_error_("Basal forcing model "<<EnumToStringx(basalforcing_model)<<" not supported yet");
		}
	}
	
	/*LATH parameters*/
	iomodel->FindConstant(&fe_FS,"md.flowequation.fe_FS");
	if(fe_FS==LATaylorHoodEnum || fe_FS==LACrouzeixRaviartEnum){
		InputUpdateFromConstantx(inputs,elements,0.,SigmaNNEnum);
	}

	/*Friction*/
	FrictionUpdateInputs(elements, inputs, iomodel);

	/*Free data: */
	iomodel->DeleteData(1,"md.flowequation.element_equation");
	xDelete<int>(finiteelement_list);
}/*}}}*/
void StressbalanceAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*Intermediaries*/
	int     fe_FS;
	int     numoutputs;
	char**  requestedoutputs = NULL;
	int     materials_type;

	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isSIA",FlowequationIsSIAEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isSSA",FlowequationIsSSAEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isL1L2",FlowequationIsL1L2Enum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isMOLHO",FlowequationIsMOLHOEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isHO",FlowequationIsHOEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isFS",FlowequationIsFSEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.fe_FS",FlowequationFeFSEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isNitscheBC",FlowequationIsNitscheEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.FSNitscheGamma",FeFSNitscheGammaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.restol",StressbalanceRestolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.reltol",StressbalanceReltolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.abstol",StressbalanceAbstolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.ishydrologylayer",StressbalanceIsHydrologyLayerEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.isnewton",StressbalanceIsnewtonEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.maxiter",StressbalanceMaxiterEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.penalty_factor",StressbalancePenaltyFactorEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.rift_penalty_threshold",StressbalanceRiftPenaltyThresholdEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.FSreconditioning",StressbalanceFSreconditioningEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.shelf_dampening",StressbalanceShelfDampeningEnum));

	/*XTH LATH parameters*/
	iomodel->FindConstant(&fe_FS,"md.flowequation.fe_FS");
	if(fe_FS==XTaylorHoodEnum || fe_FS==LATaylorHoodEnum || fe_FS==LACrouzeixRaviartEnum){
		parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.augmented_lagrangian_r",AugmentedLagrangianREnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.augmented_lagrangian_rlambda",AugmentedLagrangianRlambdaEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.XTH_theta",AugmentedLagrangianThetaEnum));
	}

	iomodel->FindConstant(&materials_type,"md.materials.type");
	if(materials_type==MatdamageiceEnum){
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.law",DamageLawEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.kappa",DamageKappaEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.damage.stress_threshold",DamageStressThresholdEnum));
	}

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.stressbalance.requested_outputs");
	parameters->AddObject(new IntParam(StressbalanceNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(StressbalanceRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.stressbalance.requested_outputs");

	/*Friction*/
	FrictionUpdateParameters(parameters, iomodel);

}/*}}}*/

/*Finite Element Analysis*/
void           StressbalanceAnalysis::Core(FemModel* femmodel){/*{{{*/

	/*Intermediaries*/
	bool isSSA,isL1L2,isMOLHO,isHO,isFS;
	bool conserve_loads = true;
	int  newton,domaintype,fe_FS;

	/* recover parameters:*/
	femmodel->parameters->FindParam(&isSSA,FlowequationIsSSAEnum);
	femmodel->parameters->FindParam(&isL1L2,FlowequationIsL1L2Enum);
	femmodel->parameters->FindParam(&isMOLHO,FlowequationIsMOLHOEnum);
	femmodel->parameters->FindParam(&isHO,FlowequationIsHOEnum);
	femmodel->parameters->FindParam(&isFS,FlowequationIsFSEnum);
	femmodel->parameters->FindParam(&fe_FS,FlowequationFeFSEnum);
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&newton,StressbalanceIsnewtonEnum);

	if(isFS && !(isSSA || isHO || isL1L2 || isMOLHO)){
		femmodel->SetCurrentConfiguration(StressbalanceAnalysisEnum);

		bool is_schur_cg_solver = false;
		#ifdef _HAVE_PETSC_
		int solver_type;
		PetscOptionsDetermineSolverType(&solver_type);

		if(solver_type==FSSolverEnum) is_schur_cg_solver = true;
		#endif

		if(is_schur_cg_solver)
		 solutionsequence_schurcg(femmodel);
		else if (fe_FS==XTaylorHoodEnum)
		 solutionsequence_la_theta(femmodel);
		else if (fe_FS==LATaylorHoodEnum || fe_FS==LACrouzeixRaviartEnum)
		 solutionsequence_la(femmodel);
		else if(newton>0)
		 solutionsequence_newton(femmodel);
		else
		 solutionsequence_nonlinear(femmodel,conserve_loads);
	}
	else if(!isFS && (isSSA || isHO || isL1L2 || isMOLHO)){
		femmodel->SetCurrentConfiguration(StressbalanceAnalysisEnum);
		if(newton>0)
		 solutionsequence_newton(femmodel);
		else
		 solutionsequence_nonlinear(femmodel,conserve_loads);

		if(domaintype==Domain2DverticalEnum && isSSA){
			femmodel->parameters->SetParam(VxEnum,InputToExtrudeEnum);
			extrudefrombase_core(femmodel);
			femmodel->parameters->SetParam(VelEnum,InputToExtrudeEnum);
			extrudefrombase_core(femmodel);
		}
	}
	else if ((isSSA || isL1L2 || isMOLHO || isHO) && isFS){
		if(VerboseSolution()) _printf0_("   computing coupling between lower order models and FS\n");
		solutionsequence_FScoupling_nonlinear(femmodel,conserve_loads);
	}
	else{
		_error_("not supported");
	}

}/*}}}*/
void           StressbalanceAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreateDVector(Element* element){/*{{{*/

	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	switch(approximation){
		case FSApproximationEnum:
			return CreateDVectorFS(element);
		default:
			return NULL; //no need for doftypes outside of FS approximation
	}
	return NULL;

}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/

	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	switch(approximation){
		case SSAApproximationEnum:
			return CreateJacobianMatrixSSA(element);
		case HOApproximationEnum:
			return CreateJacobianMatrixHO(element);
		case FSApproximationEnum:
			return CreateJacobianMatrixFS(element);
		case NoneApproximationEnum:
			return NULL;
		default:
			_error_("Approximation "<<EnumToStringx(approximation)<<" not supported");
	}
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrix(Element* element){/*{{{*/
	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	switch(approximation){
		case SIAApproximationEnum:
			return NULL;
		case SSAApproximationEnum:
			return CreateKMatrixSSA(element);
		case L1L2ApproximationEnum:
			return CreateKMatrixL1L2(element);
		case MOLHOApproximationEnum:
			return CreateKMatrixMOLHO(element);
		case HOApproximationEnum:
			return CreateKMatrixHO(element);
		case FSApproximationEnum:
			return CreateKMatrixFS(element);
		case SSAHOApproximationEnum:
			return CreateKMatrixSSAHO(element);
		case HOFSApproximationEnum:
			return CreateKMatrixHOFS(element);
		case SSAFSApproximationEnum:
			return CreateKMatrixSSAFS(element);
		case NoneApproximationEnum:
			return NULL;
		default:
			_error_("Approximation "<<EnumToStringx(approximation)<<" not supported");
	}
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVector(Element* element){/*{{{*/

	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	switch(approximation){
		case SIAApproximationEnum:
			return NULL;
		case SSAApproximationEnum:
			return CreatePVectorSSA(element);
		case L1L2ApproximationEnum:
			return CreatePVectorL1L2(element);
		case MOLHOApproximationEnum:
			return CreatePVectorMOLHO(element);
		case HOApproximationEnum:
			return CreatePVectorHO(element);
		case FSApproximationEnum:
			return CreatePVectorFS(element);
		case SSAHOApproximationEnum:
			return CreatePVectorSSAHO(element);
		case HOFSApproximationEnum:
			return CreatePVectorHOFS(element);
		case SSAFSApproximationEnum:
			return CreatePVectorSSAFS(element);
		case NoneApproximationEnum:
			return NULL;
		default:
			_error_("Approximation "<<EnumToStringx(approximation)<<" not supported");
	}
}/*}}}*/
void           StressbalanceAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	switch(approximation){
		case FSApproximationEnum: case NoneApproximationEnum:
			GetSolutionFromInputsFS(solution,element);
			return;
		case SSAApproximationEnum: case HOApproximationEnum: case L1L2ApproximationEnum: case SIAApproximationEnum:
			GetSolutionFromInputsHoriz(solution,element);
			return;
		case MOLHOApproximationEnum:
			GetSolutionFromInputsMOLHO(solution,element);
			return;
		case SSAHOApproximationEnum: case HOFSApproximationEnum: case SSAFSApproximationEnum:
			/*the elements around will create the solution*/
			return;
		default:
			_error_("Approximation "<<EnumToStringx(approximation)<<"("<<approximation<<") not supported");
	}
}/*}}}*/
void           StressbalanceAnalysis::GetSolutionFromInputsHoriz(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	IssmDouble   vx,vy;
	int          domaintype,dim,approximation,dofpernode;
	int*         doflist = NULL;

	/*Get some parameters*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; dofpernode = 2; break;
		case Domain2DverticalEnum:   dim = 2; dofpernode = 1; break;
		case Domain3DEnum:           dim = 3; dofpernode = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*dofpernode;
	element->GetInputValue(&approximation,ApproximationEnum);

	/*Fetch dof list and allocate solution vector*/
	element->GetDofList(&doflist,approximation,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numdof);

	/*Get inputs*/
	Input* vx_input=element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=NULL;
	if(domaintype!=Domain2DverticalEnum){vy_input=element->GetInput(VyEnum); _assert_(vy_input);}

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	Gauss* gauss=element->NewGauss();
	for(int i=0;i<numnodes;i++){
		gauss->GaussNode(element->FiniteElement(),i);

		/*Recover vx and vy*/
		vx_input->GetInputValue(&vx,gauss);
		values[i*dofpernode+0]=vx;
		if(dofpernode==2){
			vy_input->GetInputValue(&vy,gauss);
			values[i*dofpernode+1]=vy;
		}
	}

	solution->SetValues(numdof,doflist,values,INS_VAL);

	/*Free resources:*/
	delete gauss;
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);
}/*}}}*/
void           StressbalanceAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	switch(approximation){
		case FSApproximationEnum: case NoneApproximationEnum:
			InputUpdateFromSolutionFS(solution,element);
			return;
		case SIAApproximationEnum:
			return;
		case SSAApproximationEnum:
			InputUpdateFromSolutionSSA(solution,element);
			return;
		case HOApproximationEnum:
			InputUpdateFromSolutionHO(solution,element);
			return;
		case L1L2ApproximationEnum:
			InputUpdateFromSolutionL1L2(solution,element);
			return;
		case MOLHOApproximationEnum:
			InputUpdateFromSolutionMOLHO(solution,element);
			return;
		case SSAHOApproximationEnum:
			InputUpdateFromSolutionSSAHO(solution,element);
			return;
		case HOFSApproximationEnum:
			InputUpdateFromSolutionHOFS(solution,element);
			return;
		case SSAFSApproximationEnum:
			InputUpdateFromSolutionSSAFS(solution,element);
			return;
		default:
			_error_("Approximation "<<EnumToStringx(approximation)<<" not supported");
	}
}/*}}}*/
void           StressbalanceAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
}/*}}}*/

/*SSA*/
ElementMatrix* StressbalanceAnalysis::CreateJacobianMatrixSSA(Element* element){/*{{{*/

	/*Intermediaries*/
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement(true);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Intermediaries */
	IssmDouble Jdet,thickness;
	IssmDouble eps1dotdphii,eps1dotdphij;
	IssmDouble eps2dotdphii,eps2dotdphij;
	IssmDouble mu_prime;
	IssmDouble epsilon[3];/* epsilon=[exx,eyy,exy];*/
	IssmDouble eps1[2],eps2[2];
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element matrix, vectors and Gaussian points*/
	ElementMatrix* Ke=this->CreateKMatrixSSA(element); //Initialize Jacobian with regular SSA (first part of the Gateau derivative)
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* thickness_input = basalelement->GetInput(ThicknessEnum);_assert_(thickness_input);
	Input* vx_input        = basalelement->GetInput(VxEnum);       _assert_(vx_input);
	Input* vy_input        = basalelement->GetInput(VyEnum);       _assert_(vy_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = basalelement->NewGauss(2);
	while(gauss->next()){

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		thickness_input->GetInputValue(&thickness, gauss);
		basalelement->StrainRateSSA(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
		basalelement->material->ViscositySSADerivativeEpsSquare(&mu_prime,&epsilon[0],gauss);
		eps1[0]=2*epsilon[0]+epsilon[1];   eps2[0]=epsilon[2];
		eps1[1]=epsilon[2];                eps2[1]=epsilon[0]+2*epsilon[1];

		IssmDouble factor = gauss->weight*Jdet*2.*mu_prime*thickness;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				eps1dotdphii=eps1[0]*dbasis[0*numnodes+i]+eps1[1]*dbasis[1*numnodes+i];
				eps1dotdphij=eps1[0]*dbasis[0*numnodes+j]+eps1[1]*dbasis[1*numnodes+j];
				eps2dotdphii=eps2[0]*dbasis[0*numnodes+i]+eps2[1]*dbasis[1*numnodes+i];
				eps2dotdphij=eps2[0]*dbasis[0*numnodes+j]+eps2[1]*dbasis[1*numnodes+j];

				Ke->values[2*numnodes*(2*i+0)+2*j+0]+=factor*eps1dotdphij*eps1dotdphii;
				Ke->values[2*numnodes*(2*i+0)+2*j+1]+=factor*eps2dotdphij*eps1dotdphii;
				Ke->values[2*numnodes*(2*i+1)+2*j+0]+=factor*eps1dotdphij*eps2dotdphii;
				Ke->values[2*numnodes*(2*i+1)+2*j+1]+=factor*eps2dotdphij*eps2dotdphii;
			}
		}
	}

	/*Transform Coordinate System*/
	basalelement->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);

	/*clean-up and return*/
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;

}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixSSA(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum: case Domain2DverticalEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement(true);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixSSAViscous(basalelement);
	ElementMatrix* Ke2=CreateKMatrixSSAFriction(basalelement);
	#ifdef LATERALFRICTION
	ElementMatrix* Ke3=CreateKMatrixSSALateralFriction(basalelement);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2,Ke3);
	delete Ke3;
	#else
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);
	#endif

	/*clean-up and return*/
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	delete Ke1;
	delete Ke2;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixSSAFriction(Element* element){/*{{{*/

	/*Return if element is inactive*/
	if(element->IsAllFloating() || !element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int         dim,domaintype;
	bool        mainlyfloating;
	int         friction_style,point1;
	IssmDouble  alpha2,Jdet,fraction1,fraction2;
	IssmDouble  gllevelset,phi=1.;
	IssmDouble *xyz_list  = NULL;
	Gauss*      gauss     = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1;break;
		case Domain2DhorizontalEnum: dim = 2;break;
		case Domain3DEnum:           dim = 2;break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*dim;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke = element->NewElementMatrix(SSAApproximationEnum);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&friction_style,GroundinglineFrictionInterpolationEnum);
	Input* surface_input    = element->GetInput(SurfaceEnum); _assert_(surface_input);
	Input* gllevelset_input = NULL;

	/*build friction object, used later on: */
	Friction* friction=new Friction(element,dim);

	/*Recover portion of element that is grounded*/
	if(!(friction_style==SubelementFriction2Enum)) phi=element->GetGroundedPortion(xyz_list);
	if(friction_style==SubelementFriction2Enum){
		gllevelset_input=element->GetInput(MaskOceanLevelsetEnum); _assert_(gllevelset_input);
		element->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating, MaskOceanLevelsetEnum,0);
	   gauss = element->NewGauss(point1,fraction1,fraction2,mainlyfloating,2);
	}
	else{
		gauss = element->NewGauss(2);
	}

	/* Start  looping on the number of gaussian points: */
	while(gauss->next()){

		friction->GetAlpha2(&alpha2,gauss);
		if(friction_style==SubelementFriction1Enum) alpha2=phi*alpha2;
		else if(friction_style==SubelementFriction2Enum){
			gllevelset_input->GetInputValue(&gllevelset, gauss);
			if(gllevelset<0.) alpha2=0.;
		}
		else if(friction_style==NoFrictionOnPartiallyFloatingEnum){
			if (phi<0.99999999) alpha2=0.;
		}
		else  _error_("friction interpolation "<<EnumToStringx(friction_style)<<" not implemented yet");

		element->NodalFunctions(basis,gauss);
		element->JacobianDeterminant(&Jdet, xyz_list,gauss);

		IssmDouble factor = alpha2*gauss->weight*Jdet;
		if(dim==2){
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[2*i*2*numnodes+2*j]       += factor*basis[i]*basis[j];
					Ke->values[(2*i+1)*2*numnodes+2*j+1] += factor*basis[i]*basis[j];
				}
			}
		}
		else{
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += factor*basis[i]*basis[j];
				}
			}
		}
	}

	/*Transform Coordinate System*/
	if(dim==2) element->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	delete gauss;
	delete friction;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixSSALateralFriction(Element* element){/*{{{*/

	/*Return if element is inactive*/
	if(!element->IsIceInElement()) return NULL;

	/*If no boundary, return NULL*/
	if(!element->IsFaceOnBoundary()) return NULL;

	/*Intermediaries*/
	IssmDouble  alpha2;
	IssmDouble  Jdet;
	int         domaintype;
	IssmDouble  icefront;
	IssmDouble *xyz_list          = NULL;
	IssmDouble *xyz_list_boundary = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype==Domain2DverticalEnum) return NULL; //not supported yet

	/*Fetch number of nodes and dof for this finite element*/
	int dim      = 2;
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*dim;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke = element->NewElementMatrix(SSAApproximationEnum);
	IssmDouble*    B  = xNew<IssmDouble>(dim*numdof);
	IssmDouble*    D  = xNewZeroInit<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetLevelCoordinates(&xyz_list_boundary,xyz_list,MeshVertexonboundaryEnum,1.);
	Input* icelevelset_input = element->GetInput(MaskIceLevelsetEnum); _assert_(icelevelset_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(xyz_list,xyz_list_boundary,3);
	while(gauss->next()){

		this->GetBSSAFriction(B,element,dim,xyz_list,gauss);
		element->JacobianDeterminantSurface(&Jdet,xyz_list_boundary,gauss);
		icelevelset_input->GetInputValue(&icefront, gauss);
		if(icefront==0.)
		 alpha2=0.;
		else
		 alpha2=2.e+12;
		IssmDouble factor = alpha2*gauss->weight*Jdet;
		for(int i=0;i<dim;i++) D[i*dim+i]=factor;

		TripleMultiply(B,dim,numdof,1,
					D,dim,dim,0,
					B,dim,numdof,0,
					&Ke->values[0],1);
	}

	/*Transform Coordinate System*/
	if(dim==2) element->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_boundary);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(D);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixSSAViscous(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int         dim,domaintype;
	IssmDouble  viscosity,thickness,Jdet;
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*dim;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix(SSAApproximationEnum);
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* thickness_input=element->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* vx_input=element->GetInput(VxEnum);               _assert_(vx_input);
	Input* vy_input    = NULL;
	if(dim==2){
		vy_input    = element->GetInput(VyEnum);       _assert_(vy_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		thickness_input->GetInputValue(&thickness, gauss);
		element->material->ViscositySSA(&viscosity,dim,xyz_list,gauss,vx_input,vy_input);

		IssmDouble factor = gauss->weight*Jdet*viscosity*thickness;
		if(dim==2){
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[2*i*2*numnodes+2*j] += factor*(
								4.*dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]
								);
					Ke->values[2*i*2*numnodes+2*j+1] += factor*(
								2.*dbasis[1*numnodes+j]*dbasis[0*numnodes+i] + dbasis[0*numnodes+j]*dbasis[1*numnodes+i]
								);
					Ke->values[(2*i+1)*2*numnodes+2*j] += factor*(
								2.*dbasis[0*numnodes+j]*dbasis[1*numnodes+i] + dbasis[1*numnodes+j]*dbasis[0*numnodes+i]
								);
					Ke->values[(2*i+1)*2*numnodes+2*j+1] += factor*(
								dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + 4.*dbasis[1*numnodes+j]*dbasis[1*numnodes+i]
								);
				}
			}
		}
		else{
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += factor*(
								4.*dbasis[0*numnodes+j]*dbasis[0*numnodes+i]
								);
				}
			}
		}
	}

	/*Transform Coordinate System*/
	if(dim==2) element->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	return Ke;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorSSA(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum: case Domain2DverticalEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*compute all load vectors for this element*/
	ElementVector* pe1=CreatePVectorSSADrivingStress(basalelement);
	ElementVector* pe2=CreatePVectorSSAFront(basalelement);
	ElementVector* pe =new ElementVector(pe1,pe2);

	/*clean-up and return*/
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	delete pe1;
	delete pe2;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorSSADrivingStress(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         dim,domaintype;
	IssmDouble  thickness,Jdet,h,h_r,slope[2],hydrologyslope[2];
	IssmDouble* xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1;break;
		case Domain2DhorizontalEnum: dim = 2;break;
		case Domain3DEnum:           dim = 2;break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and vectors*/
	ElementVector* pe    = element->NewElementVector(SSAApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input*     thickness_input=element->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input*     surface_input  =element->GetInput(SurfaceEnum);   _assert_(surface_input);
	Input*     hydrologysheetthickness_input;
	Input*     hr_input;
	Input*     h_input;
	IssmDouble rhog = element->FindParam(MaterialsRhoIceEnum)*element->FindParam(ConstantsGEnum);
	bool       ishydrologylayer;
	element->FindParam(&ishydrologylayer,StressbalanceIsHydrologyLayerEnum);
	if(ishydrologylayer){
		hr_input  = element->GetInput(HydrologyBumpHeightEnum);       _assert_(hr_input);
		h_input   = element->GetInput(HydrologySheetThicknessEnum);   _assert_(h_input);
	}

	/* Start  looping on the number of gaussian points: */
	#ifndef DISCSLOPE
	Gauss* gauss=element->NewGauss(2);
	#else
	Gauss* gauss=NULL;
	Gauss* gauss_subelem=NULL;
	bool partly_floating=false;
	bool mainlyfloating=false;
	int point1;
	IssmDouble fraction1,fraction2;
	IssmDouble phi=element->GetGroundedPortion(xyz_list);
	if(phi>0.00000001 && phi<0.99999999) partly_floating=true;

	int ig=-1;// needed for driving stress parameterization
	if(partly_floating){
		element->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating, MaskOceanLevelsetEnum,0);
	   gauss=element->NewGauss(point1,fraction1,fraction2,3); //considering the entire element
		gauss_subelem=element->NewGauss(fraction1,fraction2,3);//gauss on each subelement
	}
	else{
		gauss=element->NewGauss(2);//original
	}
	#endif

	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis, gauss);

		thickness_input->GetInputValue(&thickness,gauss); _assert_(thickness>0);
		surface_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);

		#ifdef DISCSLOPE
		if(gauss_subelem && partly_floating){
			ig++;
			gauss_subelem->next();
			_assert_(std::abs(gauss_subelem->weight-gauss->weight)<0.0000001);
			/*Compute the discontinuous surface slope for this gauss point/subelement*/
			this->ComputeSurfaceSlope(&slope[0],gauss_subelem,gauss,point1,fraction1,fraction2,ig,dim,element);
			this->ComputeHydrologySlope(&hydrologyslope[0],gauss_subelem,gauss,point1,fraction1,fraction2,ig,dim,element);
		}
		#endif

		/*Change slope based on hydrology model if need be*/
		if(ishydrologylayer){
			hr_input->GetInputValue(&h_r,gauss);
			h_input->GetInputValue(&h,gauss);
			if(h>h_r){
				h_input->GetInputDerivativeValue(&hydrologyslope[0],xyz_list,gauss);
				slope[0] += hydrologyslope[0];
				if(dim==2) slope[1] += hydrologyslope[1];
			}
		}

		IssmDouble factor = rhog*thickness*Jdet*gauss->weight;
		for(int i=0;i<numnodes;i++){
			pe->values[i*dim+0]+=-factor*slope[0]*basis[i];
			if(dim==2) pe->values[i*dim+1]+=-factor*slope[1]*basis[i];
		}
	}

	/*Transform coordinate system*/
	if(dim==2) element->TransformLoadVectorCoord(pe,XYEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	#ifdef DISCSLOPE
	if(gauss_subelem) delete gauss_subelem;
	#endif
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorSSAFront(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*If no front, return NULL*/
	if(!element->IsIcefront()) return NULL;

	/*Intermediaries*/
	int         dim,domaintype;
	IssmDouble  Jdet,thickness,base,sealevel,water_pressure,ice_pressure;
	IssmDouble  surface_under_water,base_under_water,pressure;
	IssmDouble *xyz_list = NULL;
	IssmDouble *xyz_list_front = NULL;
	IssmDouble  normal[2];

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1;break;
		case Domain2DhorizontalEnum: dim = 2;break;
		case Domain3DEnum:           dim = 2;break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector(SSAApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	Input* thickness_input = element->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* base_input       = element->GetInput(BaseEnum);       _assert_(base_input);
	Input* sealevel_input       = element->GetInput(SealevelEnum);       _assert_(sealevel_input);
	IssmDouble rho_water   = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice     = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble gravity     = element->FindParam(ConstantsGEnum);
	element->GetVerticesCoordinates(&xyz_list);
	element->GetIcefrontCoordinates(&xyz_list_front,xyz_list,MaskIceLevelsetEnum);
	element->NormalSection(&normal[0],xyz_list_front);

	/*Start looping on Gaussian points*/
	Gauss* gauss=element->NewGauss(xyz_list,xyz_list_front,3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		sealevel_input->GetInputValue(&sealevel,gauss);
		base_input->GetInputValue(&base,gauss);
		element->JacobianDeterminantSurface(&Jdet,xyz_list_front,gauss);
		element->NodalFunctions(basis,gauss);

		surface_under_water = min(0.,thickness+base-sealevel); // 0 if the top of the glacier is above water level
		base_under_water    = min(0.,base-sealevel);           // 0 if the bottom of the glacier is above water level
		water_pressure = 1.0/2.0*gravity*rho_water*(surface_under_water*surface_under_water - base_under_water*base_under_water);
		ice_pressure   = 1.0/2.0*gravity*rho_ice*thickness*thickness;
		pressure = ice_pressure + water_pressure;

		IssmDouble factor = Jdet*gauss->weight*pressure;
		for(int i=0;i<numnodes;i++){
			pe->values[dim*i+0]+= factor*normal[0]*basis[i];
			if(dim==2) pe->values[dim*i+1]+= factor*normal[1]*basis[i];
		}
	}

	/*Transform coordinate system*/
	if(dim==2) element->TransformLoadVectorCoord(pe,XYEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_front);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           StressbalanceAnalysis::GetBSSA(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3] where Bi is of size 3*2.
	 * For node i, Bi can be expressed in the actual coordinate system
	 * by:
	 *                   2D                      1D
	 *       Bi=[ dN/dx           0    ]   Bi=[ dN/dx ]
	 *          [   0           dN/dy  ]
	 *          [ 1/2*dN/dy  1/2*dN/dx ]
	 * where N is the finiteelement function for node i.
	 *
	 * We assume B has been allocated already, of size: 3x(2*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* dbasis=xNew<IssmDouble>(dim*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B: */
	if(dim==2){
		for(int i=0;i<numnodes;i++){
			B[2*numnodes*0+2*i+0] = dbasis[0*numnodes+i];
			B[2*numnodes*0+2*i+1] = 0.;
			B[2*numnodes*1+2*i+0] = 0.;
			B[2*numnodes*1+2*i+1] = dbasis[1*numnodes+i];
			B[2*numnodes*2+2*i+0] = .5*dbasis[1*numnodes+i];
			B[2*numnodes*2+2*i+1] = .5*dbasis[0*numnodes+i];
		}
	}
	else{
		for(int i=0;i<numnodes;i++){
			B[i] = dbasis[i];
		}
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetBSSAFriction(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3] where Bi is square and of size 2.
	 * For node i, Bi can be expressed in the actual coordinate system
	 * by:
	 *                       2D             1D
	 *                 Bi=[ N   0 ]    Bi = N
	 *                    [ 0   N ]
	 * where N is the finiteelement function for node i.
	 *
	 * We assume B has been allocated already, of size: 2 x (numdof*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* basis=xNew<IssmDouble>(numnodes);
	element->NodalFunctions(basis,gauss);

	/*Build L: */
	if(dim==2){
		for(int i=0;i<numnodes;i++){
			B[2*numnodes*0+2*i+0] = basis[i];
			B[2*numnodes*0+2*i+1] = 0.;
			B[2*numnodes*1+2*i+0] = 0.;
			B[2*numnodes*1+2*i+1] = basis[i];
		}
	}
	else{
		for(int i=0;i<numnodes;i++){
			B[i] = basis[i];
		}
	}

	/*Clean-up*/
	xDelete<IssmDouble>(basis);
}/*}}}*/
void           StressbalanceAnalysis::GetBSSAprime(IssmDouble* Bprime,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B'  matrix. B'=[B1' B2' B3'] where Bi' is of size 3*2.
	 * For node i, Bi' can be expressed in the actual coordinate system
	 * by:
	 *                         2D                        1D
	 *       Bi_prime=[ 2*dN/dx    dN/dy ]     Bi_prime=[ 2*dN/dx ]
	 *                [   dN/dx  2*dN/dy ]
	 *                [   dN/dy    dN/dx ]
	 * where hNis the finiteelement function for node i.
	 *
	 * We assume B' has been allocated already, of size: 3x(2*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* dbasis=xNew<IssmDouble>(dim*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B': */
	if(dim==2){
		for(int i=0;i<numnodes;i++){
			Bprime[2*numnodes*0+2*i+0] = 2.*dbasis[0*numnodes+i];
			Bprime[2*numnodes*0+2*i+1] =    dbasis[1*numnodes+i];
			Bprime[2*numnodes*1+2*i+0] =    dbasis[0*numnodes+i];
			Bprime[2*numnodes*1+2*i+1] = 2.*dbasis[1*numnodes+i];
			Bprime[2*numnodes*2+2*i+0] =    dbasis[1*numnodes+i];
			Bprime[2*numnodes*2+2*i+1] =    dbasis[0*numnodes+i];
		}
	}
	else{
		for(int i=0;i<numnodes;i++){
			Bprime[i] = 2.*dbasis[i];
		}
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);
}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolutionSSA(IssmDouble* solution,Element* element){/*{{{*/

	int         i,dim,domaintype;
	IssmDouble  rho_ice,g;
	int*        doflist=NULL;
	IssmDouble* xyz_list=NULL;
	Element*    basalelement=NULL;

	/*Deal with pressure first*/
	int numvertices = element->GetNumberOfVertices();
	IssmDouble* pressure  = xNew<IssmDouble>(numvertices);
	IssmDouble* thickness = xNew<IssmDouble>(numvertices);
	IssmDouble* surface   = xNew<IssmDouble>(numvertices);

	element->FindParam(&domaintype,DomainTypeEnum);
	rho_ice =element->FindParam(MaterialsRhoIceEnum);
	g       =element->FindParam(ConstantsGEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->GetInputListOnVertices(thickness,ThicknessEnum);
			for(i=0;i<numvertices;i++) pressure[i]=rho_ice*g*thickness[i];
			dim=2;
			break;
		case Domain3DEnum:
			element->GetVerticesCoordinates(&xyz_list);
			element->GetInputListOnVertices(surface,SurfaceEnum);
			for(i=0;i<numvertices;i++) pressure[i]=rho_ice*g*(surface[i]-xyz_list[i*3+2]);
			dim=2;
			break;
		case Domain2DverticalEnum:
			element->GetVerticesCoordinates(&xyz_list);
			element->GetInputListOnVertices(surface,SurfaceEnum);
			for(i=0;i<numvertices;i++) pressure[i]=rho_ice*g*(surface[i]-xyz_list[i*3+1]);
			dim=1;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	element->AddInput(PressureEnum,pressure,P1Enum);
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(thickness);
	xDelete<IssmDouble>(surface);

	/*Get basal element*/
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum: case Domain2DverticalEnum:
			if(!element->IsOnBase()){xDelete<IssmDouble>(xyz_list); return;}
			basalelement=element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();
	int numdof   = numnodes*dim;

	/*Fetch dof list and allocate solution vectors*/
	basalelement->GetDofListLocal(&doflist,SSAApproximationEnum,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numdof);
	IssmDouble* vx     = xNew<IssmDouble>(numnodes);
	IssmDouble* vy     = xNew<IssmDouble>(numnodes);
	IssmDouble* vz     = xNew<IssmDouble>(numnodes);
	IssmDouble* vel    = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numdof;i++) values[i]=solution[doflist[i]];

	/*Transform solution in Cartesian Space*/
	if(dim==2) basalelement->TransformSolutionCoord(&values[0],XYEnum);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	for(i=0;i<numnodes;i++){
		vx[i]=values[i*dim+0];
		if(xIsNan<IssmDouble>(vx[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vx[i])) _error_("Inf found in solution vector");

		if(dim==2){
			vy[i]=values[i*dim+1];
			if(xIsNan<IssmDouble>(vy[i])) _error_("NaN found in solution vector");
			if(xIsInf<IssmDouble>(vy[i])) _error_("Inf found in solution vector");
		}
	}

	/*Get Vz and compute vel*/
	if(dim==2){
		basalelement->GetInputListOnNodes(&vz[0],VzEnum,0.);
		for(i=0;i<numnodes;i++) vel[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	}
	else{
		basalelement->GetInputListOnNodes(&vy[0],VyEnum,0.);
		for(i=0;i<numnodes;i++) vel[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
	}

	/*Add vx and vy as inputs to the tria element: */
	/*Also add surface vx and vy for the misfit, and base vx and vy for friction*/
	element->AddBasalInput(VxEnum,vx,element->GetElementType());
	element->AddBasalInput(VxSurfaceEnum,vx,element->GetElementType());
	element->AddBasalInput(VxBaseEnum,vx,element->GetElementType());
	if(dim==2) {
		element->AddBasalInput(VyEnum,vy,element->GetElementType());
		element->AddBasalInput(VySurfaceEnum,vy,element->GetElementType());
		element->AddBasalInput(VyBaseEnum,vy,element->GetElementType());
	}
	element->AddBasalInput(VelEnum,vel,element->GetElementType());

	/*Free resources:*/
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vz);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
#ifdef DISCSLOPE
void StressbalanceAnalysis::ComputeSurfaceSlope(IssmDouble* slope,Gauss* gauss_DG,Gauss* gauss_CG,int point1,IssmDouble fraction1,IssmDouble fraction2,int ig,int dim,Element* element){/*{{{*/

	/*Compute the surface slope for each subelement, for a given integration point (gauss)*/
	int numnodes=element->GetNumberOfNodes();
	IssmDouble rho_ice=element->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water=element->FindParam(MaterialsRhoSeawaterEnum); //ocean
	IssmDouble* H=xNew<IssmDouble>(numnodes);
	IssmDouble* S=xNew<IssmDouble>(numnodes);
	IssmDouble* S_subelem=xNew<IssmDouble>(numnodes);
	IssmDouble H_f1,H_f2,S_f1,S_f2;

	/*Get nodal deriviatives of the subelements related to the Gauss point*/
	IssmDouble* dbasis_subelem=xNew<IssmDouble>(dim*numnodes);//CG basis for each subelement
	this->NodalFunctionsDerivativesRGB(dbasis_subelem,gauss_DG,gauss_CG,point1,fraction1,fraction2,ig,dim,element);

	/*Define thickness at the grounding line (on element edges)*/
	element->GetInputListOnNodes(H,ThicknessEnum); // thickness on vertices
	switch(point1){//{{{
		case 0:
			H_f1=H[0]+(H[1]-H[0])*fraction1;
			H_f2=H[0]+(H[2]-H[0])*fraction2;
			break;
		case 1:
			H_f1=H[1]+(H[2]-H[1])*fraction1;
			H_f2=H[1]+(H[0]-H[1])*fraction2;
			break;
		case 2:
			H_f1=H[2]+(H[0]-H[2])*fraction1;
			H_f2=H[2]+(H[1]-H[2])*fraction2;
			break;
		default:
			_error_("index "<<point1<<" not supported yet");
	}//}}}

	/*Define surface at the grounding (on element edges)*/
	S_f1=H_f1*(1-rho_ice/rho_water);
	S_f2=H_f2*(1-rho_ice/rho_water);
	element->GetInputListOnNodes(S,SurfaceEnum); // surface on vertices

	/*Define surface on the subelement vertices*/
   if(ig<4){ // BLUE element itapopo only if order is = 3
		switch(point1){ //{{{
			case 0:
				S_subelem[0]=S[0];
				S_subelem[1]=S_f1;
				S_subelem[2]=S_f2;
				break;
			case 1:
				S_subelem[0]=S[1];
				S_subelem[1]=S_f1;
				S_subelem[2]=S_f2;
				break;
			case 2:
				S_subelem[0]=S[2];
				S_subelem[1]=S_f1;
				S_subelem[2]=S_f2;
				break;
			default:
				_error_("index "<<point1<<" not supported yet");
		}//}}}
	}
	if(ig>3 && ig<8){ // GREEN element
		switch(point1){ //{{{
			case 0:
				S_subelem[0]=S_f1;
				S_subelem[1]=S[2];
				S_subelem[2]=S_f2;
				break;
			case 1:
				S_subelem[0]=S_f1;
				S_subelem[1]=S[0];
				S_subelem[2]=S_f2;
				break;
			case 2:
				S_subelem[0]=S_f1;
				S_subelem[1]=S[1];
				S_subelem[2]=S_f2;
				break;
			default:
				_error_("index "<<point1<<" not supported yet");
		}//}}}
	}
	if(ig>7){ // RED element
		switch(point1){ //{{{
			case 0:
				S_subelem[0]=S_f1;
				S_subelem[1]=S[1];
				S_subelem[2]=S[2];
				break;
			case 1:
				S_subelem[0]=S_f1;
				S_subelem[1]=S[2];
				S_subelem[2]=S[0];
				break;
			case 2:
				S_subelem[0]=S_f1;
				S_subelem[1]=S[0];
				S_subelem[2]=S[1];
				break;
			default:
				_error_("index "<<point1<<" not supported yet");
		}//}}}
	}

	//_printf_("\t"<<"H[1]\t"<<"H[2]\t"<<"H[3]\t"<<"S[1]\t"<<"S[2]\t"<<"S[3]\n");
	//_printf_("\t"<<H[0]<<"\t"<<H[1]<<"\t"<<H[2]<<"\t"<<S[0]<<"\t"<<S[1]<<"\t"<<S[2]<<"\n");
	//_printf_("\t H_f1 \t H_f2 \t S_f1 \t S_f2 \n");
	//_printf_("\t"<< H_f1 <<"\t"<< H_f2 << "\t" << S_f1 << "\t" << S_f2<< "\n");
	//_printf_("\t"<<"S_subelem[1]\t"<<"S_subelem[2]\t"<<"S_subelem[3]\n");
	//_printf_("\t"<<S_subelem[0]<<"\t"<<S_subelem[1]<<"\t"<<S_subelem[2]<<"\n");

	/*Compute slope*/
	slope[0]=0;
	slope[1]=0;
	for(int i=0;i<numnodes;i++){
		slope[0]+=S_subelem[i]*dbasis_subelem[numnodes*0+i]; //dSdx
		slope[1]+=S_subelem[i]*dbasis_subelem[numnodes*1+i]; //dSdy
	}

	/*Clean up*/
	xDelete<IssmDouble>(dbasis_subelem);
	xDelete<IssmDouble>(H);
	xDelete<IssmDouble>(S);
	xDelete<IssmDouble>(S_subelem);

}/*}}}*/
void StressbalanceAnalysis::ComputeHydrologySlope(IssmDouble* hydrologyslope,Gauss* gauss_DG,Gauss* gauss_CG,int point1,IssmDouble fraction1,IssmDouble fraction2,int ig,int dim,Element* element){/*{{{*/

	/*Compute the surface slope for each subelement, for a given integration point (gauss)*/
	int numnodes=element->GetNumberOfNodes();
	IssmDouble* h=xNew<IssmDouble>(numnodes);
	IssmDouble* Hl=xNew<IssmDouble>(numnodes);
	IssmDouble* h_r=xNew<IssmDouble>(numnodes);
	IssmDouble* Hl_subelem=xNew<IssmDouble>(numnodes);
	IssmDouble Hl_f1,Hl_f2;

	/*Get nodal deriviatives of the subelements related to the Gauss point*/
	IssmDouble* dbasis_subelem=xNew<IssmDouble>(dim*numnodes);//CG basis for each subelement
	this->NodalFunctionsDerivativesRGB(dbasis_subelem,gauss_DG,gauss_CG,point1,fraction1,fraction2,ig,dim,element);

	/*Define hydrology layer at the grounding line (on element edges)*/
	element->GetInputListOnNodes(h,HydrologySheetThicknessEnum); // hydrology sheet thickness on vertices
	element->GetInputListOnNodes(h_r,HydrologySheetThicknessEnum); // bedrock bump height on vertices
	for (int i=0;i<numnodes;i++){
		if(h[i]<h_r[i]){
			Hl[i] = 0;
		}
		else{
			Hl[i] = h[i] - h_r[i];
		}
	}
	switch(point1){//{{{
		case 0:
			Hl_f1=H[0]+(Hl[1]-Hl[0])*fraction1;
			Hl_f2=H[0]+(Hl[2]-Hl[0])*fraction2;
			break;
		case 1:
			Hl_f1=Hl[1]+(Hl[2]-Hl[1])*fraction1;
			Hl_f2=Hl[1]+(Hl[0]-Hl[1])*fraction2;
			break;
		case 2:
			Hl_f1=Hl[2]+(Hl[0]-Hl[2])*fraction1;
			Hl_f2=Hl[2]+(Hl[1]-Hl[2])*fraction2;
			break;
		default:
			_error_("index "<<point1<<" not supported yet");
	}//}}}


	/*Define surface on the subelement vertices*/
   if(ig<4){ // BLUE element itapopo only if order is = 3
		switch(point1){ //{{{
			case 0:
				Hl_subelem[0]=Hl[0];
				Hl_subelem[1]=Hl_f1;
				Hl_subelem[2]=Hl_f2;
				break;
			case 1:
				Hl_subelem[0]=Hl[1];
				Hl_subelem[1]=Hl_f1;
				Hl_subelem[2]=Hl_f2;
				break;
			case 2:
				Hl_subelem[0]=Hl[2];
				Hl_subelem[1]=Hl_f1;
				Hl_subelem[2]=Hl_f2;
				break;
			default:
				_error_("index "<<point1<<" not supported yet");
		}//}}}
	}
	if(ig>3 && ig<8){ // GREEN element
		switch(point1){ //{{{
			case 0:
				Hl_subelem[0]=Hl_f1;
				Hl_subelem[1]=Hl[2];
				Hl_subelem[2]=Hl_f2;
				break;
			case 1:
				Hl_subelem[0]=Hl_f1;
				Hl_subelem[1]=Hl[0];
				Hl_subelem[2]=Hl_f2;
				break;
			case 2:
				Hl_subelem[0]=Hl_f1;
				Hl_subelem[1]=Hl[1];
				Hl_subelem[2]=Hl_f2;
				break;
			default:
				_error_("index "<<point1<<" not supported yet");
		}//}}}
	}
	if(ig>7){ // RED element
		switch(point1){ //{{{
			case 0:
				Hl_subelem[0]=Hl_f1;
				Hl_subelem[1]=Hl[1];
				Hl_subelem[2]=Hl[2];
				break;
			case 1:
				Hl_subelem[0]=Hl_f1;
				Hl_subelem[1]=Hl[2];
				Hl_subelem[2]=Hl[0];
				break;
			case 2:
				Hl_subelem[0]=Hl_f1;
				Hl_subelem[1]=Hl[0];
				Hl_subelem[2]=Hl[1];
				break;
			default:
				_error_("index "<<point1<<" not supported yet");
		}//}}}
	}

	/*Compute slope*/
	hydrologyslope[0]=0;
	hydrologyslope[1]=0;
	for(int i=0;i<numnodes;i++){
		hydrologyslope[0]+=Hl_subelem[i]*dbasis_subelem[numnodes*0+i]; //dHldx
		hydrologyslope[1]+=Hl_subelem[i]*dbasis_subelem[numnodes*1+i]; //dHdy
	}

	/*Clean up*/
	xDelete<IssmDouble>(dbasis_subelem);
	xDelete<IssmDouble>(h);
	xDelete<IssmDouble>(Hl);
	xDelete<IssmDouble>(h_r);
	xDelete<IssmDouble>(Sl_subelem);

}/*}}}*/
void StressbalanceAnalysis::NodalFunctionsDerivativesRGB(IssmDouble* dbasis_subelem,Gauss* gauss_DG,Gauss* gauss_CG,int point1,IssmDouble fraction1,IssmDouble fraction2,int ig,int dim,Element* element){/*{{{*/

	/*Fetch number of nodes for this finite element*/
   int numnodes = element->GetNumberOfNodes();
   IssmDouble dbasis_ref[dim*numnodes];
   IssmDouble Jinv[2][2];
   //IssmDouble* dbasis_subelem = xNew<IssmDouble>(dim*numnodes);//CG basis for each subelement
   IssmDouble* xyz_list_RGB = xNew<IssmDouble>(3*numnodes); // x,y,z per node
   IssmDouble* xyz_list = NULL;

   element->GetVerticesCoordinates(&xyz_list);
   Tria* tria = xDynamicCast<Tria*>(element);

   /*Get nodal functions derivatives*/
   tria->GetNodalFunctionsDerivativesReference(dbasis_ref,gauss_DG,P1Enum); //gauss is not used here

   if(ig<4){ // BLUE element itapopo only if order is = 3
      switch(point1){//{{{
         case 0:
            /*Vertex 0 - local in the subelement*/
            xyz_list_RGB[3*0+0]=xyz_list[3*0+0]; //x
            xyz_list_RGB[3*0+1]=xyz_list[3*0+1]; //y
            /*Vertex 1*/
            xyz_list_RGB[3*1+0]=xyz_list[3*0+0]+(xyz_list[3*1+0]-xyz_list[3*0+0])*fraction1; //x
            xyz_list_RGB[3*1+1]=xyz_list[3*0+1]+(xyz_list[3*1+1]-xyz_list[3*0+1])*fraction1; //y
            /*Vertex 2*/
            xyz_list_RGB[3*2+0]=xyz_list[3*0+0]+(xyz_list[3*2+0]-xyz_list[3*0+0])*fraction2; //x;
            xyz_list_RGB[3*2+1]=xyz_list[3*0+1]+(xyz_list[3*2+1]-xyz_list[3*0+1])*fraction2; //y;
            break;
			 case 1:
            /*Vertex 0 - local in the subelement*/
            xyz_list_RGB[3*0+0]=xyz_list[3*1+0]; //x
            xyz_list_RGB[3*0+1]=xyz_list[3*1+1]; //y
            /*Vertex 1*/
            xyz_list_RGB[3*1+0]=xyz_list[3*1+0]+(xyz_list[3*2+0]-xyz_list[3*1+0])*fraction1; //x
            xyz_list_RGB[3*1+1]=xyz_list[3*1+1]+(xyz_list[3*2+1]-xyz_list[3*1+1])*fraction1; //y
            /*Vertex 2*/
            xyz_list_RGB[3*2+0]=xyz_list[3*1+0]+(xyz_list[3*0+0]-xyz_list[3*1+0])*fraction2; //x;
            xyz_list_RGB[3*2+1]=xyz_list[3*1+1]+(xyz_list[3*0+1]-xyz_list[3*1+1])*fraction2; //y;
            break;
			case 2:
            /*Vertex 0 - local in the subelement*/
            xyz_list_RGB[3*0+0]=xyz_list[3*2+0]; //x
            xyz_list_RGB[3*0+1]=xyz_list[3*2+1]; //y
            /*Vertex 1*/
            xyz_list_RGB[3*1+0]=xyz_list[3*2+0]+(xyz_list[3*0+0]-xyz_list[3*2+0])*fraction1; //x
            xyz_list_RGB[3*1+1]=xyz_list[3*2+1]+(xyz_list[3*0+1]-xyz_list[3*2+1])*fraction1; //y
            /*Vertex 2*/
            xyz_list_RGB[3*2+0]=xyz_list[3*2+0]+(xyz_list[3*1+0]-xyz_list[3*2+0])*fraction2; //x;
            xyz_list_RGB[3*2+1]=xyz_list[3*2+1]+(xyz_list[3*1+1]-xyz_list[3*2+1])*fraction2; //y;
            break;
         default:
            _error_("index "<<point1<<" not supported yet");
      }//}}}
   }
	if(ig>3 && ig<8){ // GREEN element
      switch(point1){//{{{
         case 0:
            /*Vertex 0 - local in the subelement*/
            xyz_list_RGB[3*0+0]=xyz_list[3*0+0]+(xyz_list[3*1+0]-xyz_list[3*0+0])*fraction1; //x
            xyz_list_RGB[3*0+1]=xyz_list[3*0+1]+(xyz_list[3*1+1]-xyz_list[3*0+1])*fraction1; //y
            /*Vertex 1*/
            xyz_list_RGB[3*1+0]=xyz_list[3*2+0]; //x
            xyz_list_RGB[3*1+1]=xyz_list[3*2+1]; //y
            /*Vertex 2*/
            xyz_list_RGB[3*2+0]=xyz_list[3*0+0]+(xyz_list[3*2+0]-xyz_list[3*0+0])*fraction2; //x;
            xyz_list_RGB[3*2+1]=xyz_list[3*0+1]+(xyz_list[3*2+1]-xyz_list[3*0+1])*fraction2; //y;
            break;
			 case 1:
            /*Vertex 0*/
            xyz_list_RGB[3*0+0]=xyz_list[3*1+0]+(xyz_list[3*2+0]-xyz_list[3*1+0])*fraction1; //x
            xyz_list_RGB[3*0+1]=xyz_list[3*1+1]+(xyz_list[3*2+1]-xyz_list[3*1+1])*fraction1; //y
            /*Vertex 1*/
            xyz_list_RGB[3*1+0]=xyz_list[3*0+0]; //x
            xyz_list_RGB[3*1+1]=xyz_list[3*0+1]; //y
            /*Vertex 2*/
            xyz_list_RGB[3*2+0]=xyz_list[3*1+0]+(xyz_list[3*0+0]-xyz_list[3*1+0])*fraction2; //x;
            xyz_list_RGB[3*2+1]=xyz_list[3*1+1]+(xyz_list[3*0+1]-xyz_list[3*1+1])*fraction2; //y;
            break;
			case 2:
            /*Vertex 0*/
            xyz_list_RGB[3*0+0]=xyz_list[3*2+0]+(xyz_list[3*0+0]-xyz_list[3*2+0])*fraction1; //x
            xyz_list_RGB[3*0+1]=xyz_list[3*2+1]+(xyz_list[3*0+1]-xyz_list[3*2+1])*fraction1; //y
            /*Vertex 1*/
            xyz_list_RGB[3*1+0]=xyz_list[3*1+0]; //x
            xyz_list_RGB[3*1+1]=xyz_list[3*1+1]; //y
            /*Vertex 2*/
            xyz_list_RGB[3*2+0]=xyz_list[3*2+0]+(xyz_list[3*1+0]-xyz_list[3*2+0])*fraction2; //x;
            xyz_list_RGB[3*2+1]=xyz_list[3*2+1]+(xyz_list[3*1+1]-xyz_list[3*2+1])*fraction2; //y;
            break;
         default:
            _error_("index "<<point1<<" not supported yet");
      }//}}}
   }
	if(ig>7){ // RED element
      switch(point1){//{{{
         case 0:
            /*Vertex 0*/
            xyz_list_RGB[3*0+0]=xyz_list[3*0+0]+(xyz_list[3*1+0]-xyz_list[3*0+0])*fraction1; //x
            xyz_list_RGB[3*0+1]=xyz_list[3*0+1]+(xyz_list[3*1+1]-xyz_list[3*0+1])*fraction1; //y
            /*Vertex 1*/
            xyz_list_RGB[3*1+0]=xyz_list[3*1+0]; //x
            xyz_list_RGB[3*1+1]=xyz_list[3*1+1]; //y
            /*Vertex 2*/
            xyz_list_RGB[3*2+0]=xyz_list[3*2+0]; //x;
            xyz_list_RGB[3*2+1]=xyz_list[3*2+1]; //y;
            break;
         case 1:
				/*Vertex 0*/
            xyz_list_RGB[3*0+0]=xyz_list[3*1+0]+(xyz_list[3*2+0]-xyz_list[3*1+0])*fraction1; //x
            xyz_list_RGB[3*0+1]=xyz_list[3*1+1]+(xyz_list[3*2+1]-xyz_list[3*1+1])*fraction1; //y
            /*Vertex 1*/
            xyz_list_RGB[3*1+0]=xyz_list[3*2+0]; //x
            xyz_list_RGB[3*1+1]=xyz_list[3*2+1]; //y
            /*Vertex 2*/
            xyz_list_RGB[3*2+0]=xyz_list[3*0+0]; //x;
            xyz_list_RGB[3*2+1]=xyz_list[3*0+1]; //y;
            break;
			case 2:
            /*Vertex 0*/
            xyz_list_RGB[3*0+0]=xyz_list[3*2+0]+(xyz_list[3*0+0]-xyz_list[3*2+0])*fraction1; //x
            xyz_list_RGB[3*0+1]=xyz_list[3*2+1]+(xyz_list[3*0+1]-xyz_list[3*2+1])*fraction1; //y
            /*Vertex 1*/
            xyz_list_RGB[3*1+0]=xyz_list[3*0+0]; //x
            xyz_list_RGB[3*1+1]=xyz_list[3*0+1]; //y
            /*Vertex 2*/
            xyz_list_RGB[3*2+0]=xyz_list[3*1+0]; //x;
            xyz_list_RGB[3*2+1]=xyz_list[3*1+1]; //y;
            break;
         default:
            _error_("index "<<point1<<" not supported yet");
      }//}}}
   }

	//std::cout<<"	ig="<<ig<<std::endl;
	//std::cout<<"	"<<xyz_list_RGB[0*numnodes+0]<<"\t"<<xyz_list_RGB[1*numnodes+0]<<"\t"<<xyz_list_RGB[2*numnodes+0]<<std::endl; //x
	//std::cout<<"	"<<xyz_list_RGB[0*numnodes+1]<<"\t"<<xyz_list_RGB[1*numnodes+1]<<"\t"<<xyz_list_RGB[2*numnodes+1]<<std::endl; //y

   /*Get Jacobian*/
	tria->GetJacobianInvert(&Jinv[0][0],xyz_list_RGB,gauss_CG); // gauss is not used here

	/*Build dbasis CG in the subelement*/
	for(int i=0;i<numnodes;i++){
		dbasis_subelem[numnodes*0+i]=Jinv[0][0]*dbasis_ref[0*numnodes+i]+Jinv[0][1]*dbasis_ref[1*numnodes+i];//dPhi_i/dx
		dbasis_subelem[numnodes*1+i]=Jinv[1][0]*dbasis_ref[0*numnodes+i]+Jinv[1][1]*dbasis_ref[1*numnodes+i];//dPhi_i/dy
	}

   xDelete<IssmDouble>(xyz_list);
   xDelete<IssmDouble>(xyz_list_RGB);
   //xDelete<IssmDouble>(dbasis_subelem);

}/*}}}*/
#endif

/*L1L2*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixL1L2(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixL1L2Viscous(element);
	ElementMatrix* Ke2=CreateKMatrixL1L2Friction(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixL1L2Friction(Element* element){/*{{{*/

	if(!element->IsOnBase() || element->IsAllFloating()) return NULL;
	Element* basalelement = element->SpawnBasalElement();
	ElementMatrix* Ke = CreateKMatrixSSAFriction(basalelement);

	/*clean-up and return*/
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixL1L2Viscous(Element* element){/*{{{*/

	/*Intermediaries*/
	IssmDouble  viscosity,Jdet;
	IssmDouble *xyz_list = NULL;

	/*Get element on base*/
	Element* basalelement = element->GetBasalElement()->SpawnBasalElement(true);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke     = basalelement->NewElementMatrix(L1L2ApproximationEnum);
	IssmDouble*    dbasis = xNew<IssmDouble>(3*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* surface_input = element->GetInput(SurfaceEnum); _assert_(surface_input);
	Input* vx_input      = element->GetInput(VxEnum);      _assert_(vx_input);
	Input* vy_input      = element->GetInput(VyEnum);      _assert_(vy_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss      = element->NewGauss(5);
	Gauss* gauss_base = basalelement->NewGauss();
	while(gauss->next()){
		gauss->SynchronizeGaussBase(gauss_base);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctionsDerivatives(dbasis,xyz_list,gauss_base);

		element->material->ViscosityL1L2(&viscosity,xyz_list,gauss,vx_input,vy_input,surface_input);

		IssmDouble factor=gauss->weight*Jdet*viscosity;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[2*i*2*numnodes+2*j] += factor*(
							4.*dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]
							);
				Ke->values[2*i*2*numnodes+2*j+1] += factor*(
							2.*dbasis[1*numnodes+j]*dbasis[0*numnodes+i] + dbasis[0*numnodes+j]*dbasis[1*numnodes+i]
							);
				Ke->values[(2*i+1)*2*numnodes+2*j] += factor*(
							2.*dbasis[0*numnodes+j]*dbasis[1*numnodes+i] + dbasis[1*numnodes+j]*dbasis[0*numnodes+i]
							);
				Ke->values[(2*i+1)*2*numnodes+2*j+1] += factor*(
							dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + 4.*dbasis[1*numnodes+j]*dbasis[1*numnodes+i]
							);
			}
		}
	}

	/*Transform Coordinate System*/
	basalelement->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	delete gauss;
	delete gauss_base;
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	return Ke;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorL1L2(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum: case Domain2DverticalEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*compute all load vectors for this element*/
	ElementVector* pe1=CreatePVectorL1L2DrivingStress(basalelement);
	ElementVector* pe2=CreatePVectorL1L2Front(basalelement);
	ElementVector* pe =new ElementVector(pe1,pe2);

	/*clean-up and return*/
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	delete pe1;
	delete pe2;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorL1L2DrivingStress(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble  thickness,Jdet,slope[2];
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and vectors*/
	ElementVector* pe    = element->NewElementVector(L1L2ApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input*     thickness_input=element->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input*     surface_input  =element->GetInput(SurfaceEnum);   _assert_(surface_input);
	IssmDouble rhog = element->FindParam(MaterialsRhoIceEnum)*element->FindParam(ConstantsGEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis, gauss);

		thickness_input->GetInputValue(&thickness,gauss);
		surface_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);

		IssmDouble factor = rhog*thickness*Jdet*gauss->weight;
		for(int i=0;i<numnodes;i++){
			pe->values[i*2+0]+=-factor*slope[0]*basis[i];
			pe->values[i*2+1]+=-factor*slope[1]*basis[i];
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,XYEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorL1L2Front(Element* element){/*{{{*/

	/*If no front, return NULL*/
	if(!element->IsIcefront()) return NULL;

	/*Intermediaries*/
	IssmDouble  Jdet,thickness,base,sealevel,water_pressure,ice_pressure;
	IssmDouble  surface_under_water,base_under_water,pressure;
	IssmDouble *xyz_list = NULL;
	IssmDouble *xyz_list_front = NULL;
	IssmDouble  normal[2];

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector(L1L2ApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	Input* thickness_input = element->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* base_input       = element->GetInput(BaseEnum);       _assert_(base_input);
	Input* sealevel_input       = element->GetInput(SealevelEnum);       _assert_(sealevel_input);
	IssmDouble rho_water   = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice     = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble gravity     = element->FindParam(ConstantsGEnum);
	element->GetVerticesCoordinates(&xyz_list);
	element->GetIcefrontCoordinates(&xyz_list_front,xyz_list,MaskIceLevelsetEnum);
	element->NormalSection(&normal[0],xyz_list_front);

	/*Start looping on Gaussian points*/
	Gauss* gauss=element->NewGauss(xyz_list,xyz_list_front,3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		base_input->GetInputValue(&base,gauss);
		sealevel_input->GetInputValue(&sealevel,gauss);
		element->JacobianDeterminantSurface(&Jdet,xyz_list_front,gauss);
		element->NodalFunctions(basis,gauss);

		surface_under_water = min(0.,thickness+base-sealevel); // 0 if the top of the glacier is above water level
		base_under_water    = min(0.,base-sealevel);           // 0 if the bottom of the glacier is above water level
		water_pressure = 1.0/2.0*gravity*rho_water*(surface_under_water*surface_under_water - base_under_water*base_under_water);
		ice_pressure   = 1.0/2.0*gravity*rho_ice*thickness*thickness;
		pressure = ice_pressure + water_pressure;

		IssmDouble factor = pressure*Jdet*gauss->weight;
		for (int i=0;i<numnodes;i++){
			pe->values[2*i+0]+= factor*normal[0]*basis[i];
			pe->values[2*i+1]+= factor*normal[1]*basis[i];
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,XYEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_front);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolutionL1L2(IssmDouble* solution,Element* element){/*{{{*/

	int         i,dim,domaintype;
	IssmDouble  rho_ice,g;
	int*        doflist=NULL;
	IssmDouble* xyz_list=NULL;
	Element*    basalelement=NULL;

	/*Deal with pressure first*/
	int numvertices = element->GetNumberOfVertices();
	IssmDouble* pressure  = xNew<IssmDouble>(numvertices);
	IssmDouble* thickness = xNew<IssmDouble>(numvertices);
	IssmDouble* surface   = xNew<IssmDouble>(numvertices);

	element->FindParam(&dim,DomainDimensionEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	rho_ice =element->FindParam(MaterialsRhoIceEnum);
	g       =element->FindParam(ConstantsGEnum);
	if(dim==2){
		element->GetInputListOnVertices(thickness,ThicknessEnum);
		for(i=0;i<numvertices;i++) pressure[i]=rho_ice*g*thickness[i];
	}
	else{
		element->GetVerticesCoordinates(&xyz_list);
		element->GetInputListOnVertices(surface,SurfaceEnum);
		for(i=0;i<numvertices;i++) pressure[i]=rho_ice*g*(surface[i]-xyz_list[i*3+2]);
	}
	element->AddInput(PressureEnum,pressure,P1Enum);
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(thickness);
	xDelete<IssmDouble>(surface);

	/*Get basal element*/
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()){xDelete<IssmDouble>(xyz_list); return;}
			basalelement=element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();
	int numdof   = numnodes*2;

	/*Fetch dof list and allocate solution vectors*/
	basalelement->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values    = xNew<IssmDouble>(numdof);
	IssmDouble* vx        = xNew<IssmDouble>(numnodes);
	IssmDouble* vy        = xNew<IssmDouble>(numnodes);
	IssmDouble* vz        = xNew<IssmDouble>(numnodes);
	IssmDouble* vel       = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numdof;i++) values[i]=solution[doflist[i]];

	/*Transform solution in Cartesian Space*/
	basalelement->TransformSolutionCoord(&values[0],XYEnum);
	basalelement->FindParam(&domaintype,DomainTypeEnum);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	for(i=0;i<numnodes;i++){
		vx[i]=values[i*2+0];
		vy[i]=values[i*2+1];

		/*Check solution*/
		if(xIsNan<IssmDouble>(vx[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vx[i])) _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(vy[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vy[i])) _error_("Inf found in solution vector");
	}

	/*Get Vz and compute vel*/
	basalelement->GetInputListOnNodes(&vz[0],VzEnum,0.);
	for(i=0;i<numnodes;i++) vel[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

	/*Add vx and vy as inputs to the tria element: */
	element->AddBasalInput(VxEnum,vx,element->GetElementType());
	element->AddBasalInput(VyEnum,vy,element->GetElementType());
	element->AddBasalInput(VelEnum,vel,element->GetElementType());

	/*Also add surface vx and vy for the misfit, and base vx and vy for friction*/
	element->AddBasalInput(VxSurfaceEnum,vx,element->GetElementType());
	element->AddBasalInput(VySurfaceEnum,vy,element->GetElementType());
	element->AddBasalInput(VxBaseEnum,vx,element->GetElementType());
	element->AddBasalInput(VyBaseEnum,vy,element->GetElementType());

	/*Free resources:*/
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vz);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/

/*MOLHO*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixMOLHO(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum: 
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement(true);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixMOLHOViscous(basalelement);
	ElementMatrix* Ke2=CreateKMatrixMOLHOFriction(basalelement);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);

	/*clean-up and return*/
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	delete Ke1;
	delete Ke2;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixMOLHOFriction(Element* element){/*{{{*/

	if(element->IsAllFloating()) return NULL;

	/*Intermediaries*/
	int         dim,domaintype;
	bool        mainlyfloating;
	int         friction_style,point1;
	IssmDouble  alpha2,Jdet,fraction1,fraction2;
	IssmDouble  gllevelset,phi=1.;
	IssmDouble *xyz_list  = NULL;
	Gauss*      gauss     = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1;break;
		case Domain2DhorizontalEnum: dim = 2;break;
		case Domain3DEnum:           dim = 2;break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*dim;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke = element->NewElementMatrix(MOLHOApproximationEnum);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&friction_style,GroundinglineFrictionInterpolationEnum);
	Input* surface_input    = element->GetInput(SurfaceEnum); _assert_(surface_input);
	Input* gllevelset_input = NULL;

	/*build friction object, used later on: */
	Friction* friction=new Friction(element, dim);

	/*Recover portion of element that is grounded*/
	if(!(friction_style==SubelementFriction2Enum)) phi=element->GetGroundedPortion(xyz_list);
	if(friction_style==SubelementFriction2Enum){
		gllevelset_input=element->GetInput(MaskOceanLevelsetEnum); _assert_(gllevelset_input);
		element->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating,MaskOceanLevelsetEnum,0);
	   gauss = element->NewGauss(point1,fraction1,fraction2,mainlyfloating,2);
	}
	else{
		gauss = element->NewGauss(2);
	}

	/* Start  looping on the number of gaussian points: */
	while(gauss->next()){

		friction->GetAlpha2(&alpha2,gauss);
		if(friction_style==SubelementFriction1Enum) alpha2=phi*alpha2;
		else if(friction_style==SubelementFriction2Enum){
			gllevelset_input->GetInputValue(&gllevelset, gauss);
			if(gllevelset<0.) alpha2=0.;
		}
		else if(friction_style==NoFrictionOnPartiallyFloatingEnum){
			if (phi<0.99999999) alpha2=0.;
		}
		else  _error_("friction interpolation "<<EnumToStringx(friction_style)<<" not implemented yet");

		element->NodalFunctions(basis,gauss);
		element->JacobianDeterminant(&Jdet, xyz_list,gauss);

		IssmDouble factor = alpha2*gauss->weight*Jdet;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[(4*i+0)*2*2*numnodes+4*j+0] += factor*basis[i]*basis[j];
				Ke->values[(4*i+2)*2*2*numnodes+4*j+2] += factor*basis[i]*basis[j];
			}
		}
	}

	/*Clean up and return*/
	delete gauss;
	delete friction;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixMOLHOViscous(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	IssmDouble  viscosity[4],Jdet,thickness,n;
	IssmDouble *xyz_list = NULL;
	int      domaintype;
	int dim=2;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix(MOLHOApproximationEnum);
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes); // like SSA
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes); // like SSA

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* thickness_input= element->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* surface_input  = element->GetInput(SurfaceEnum);   _assert_(surface_input);
	Input* vxbase_input   = element->GetInput(VxBaseEnum);    _assert_(vxbase_input); 
	Input* vybase_input   = element->GetInput(VyBaseEnum);    _assert_(vybase_input);
	Input* vxshear_input  = element->GetInput(VxShearEnum);   _assert_(vxshear_input); //shear vx at the surface 
	Input* vyshear_input  = element->GetInput(VyShearEnum);   _assert_(vyshear_input); //shear vy at the surface
	Input* n_input			 = element->GetInput(MaterialsRheologyNEnum); _assert_(n_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss      = element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		element->NodalFunctions(basis, gauss);

		thickness_input->GetInputValue(&thickness, gauss);
		n_input->GetInputValue(&n,gauss);
		element->material->ViscosityMOLHO(&viscosity[0],dim,xyz_list,gauss,vxbase_input,vybase_input,vxshear_input,vyshear_input,thickness_input,n_input);

		IssmDouble factor = gauss->weight*Jdet*thickness;
		for(int i=0;i<numnodes;i++){//shape functions on tria element
			for(int j=0;j<numnodes;j++){
				Ke->values[(4*i+0)*2*2*numnodes+4*j+0] += factor*viscosity[0]*(
							4.*dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]
							);//K11
				Ke->values[(4*i+0)*2*2*numnodes+4*j+1] += factor*viscosity[1]*(
							4.*dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]
							)*(n+1)/(n+2);//K12
				Ke->values[(4*i+0)*2*2*numnodes+4*j+2] += factor*viscosity[0]*(
							2.*dbasis[1*numnodes+j]*dbasis[0*numnodes+i] + dbasis[0*numnodes+j]*dbasis[1*numnodes+i]
							);//K13
				Ke->values[ (4*i+0)*2*2*numnodes+4*j+3] += factor*viscosity[1]*(
                     2.*dbasis[1*numnodes+j]*dbasis[0*numnodes+i] + dbasis[0*numnodes+j]*dbasis[1*numnodes+i]
                     )*(n+1)/(n+2);//K14

				Ke->values[(4*i+1)*2*2*numnodes+4*j+0] += factor*viscosity[1]*(
							4.*dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]
							)*(n+1)/(n+2);//K21
				Ke->values[(4*i+1)*2*2*numnodes+4*j+1] += factor*viscosity[2]*(
							4.*dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]
							)*2*pow(n+1,2)/((2*n+3)*(n+2)) 
							+ 
							gauss->weight*Jdet*viscosity[3]*basis[j]*basis[i]*pow(n+1,2)/(thickness*(2*n+1))
							;//K22
				Ke->values[(4*i+1)*2*2*numnodes+4*j+2] += factor*viscosity[1]*(
							2.*dbasis[1*numnodes+j]*dbasis[0*numnodes+i] + dbasis[0*numnodes+j]*dbasis[1*numnodes+i]
							)*(n+1)/(n+2);//K23
				Ke->values[(4*i+1)*2*2*numnodes+4*j+3] += factor*viscosity[2]*(
							2.*dbasis[1*numnodes+j]*dbasis[0*numnodes+i] + dbasis[0*numnodes+j]*dbasis[1*numnodes+i]
							)*2*pow(n+1,2)/((2*n+3)*(n+2));//K24

				Ke->values[(4*i+2)*2*2*numnodes+4*j+0] += factor*viscosity[0]*(
							2.*dbasis[0*numnodes+j]*dbasis[1*numnodes+i] + dbasis[1*numnodes+j]*dbasis[0*numnodes+i]
							);//K31
				Ke->values[(4*i+2)*2*2*numnodes+4*j+1] += factor*viscosity[1]*(
							2.*dbasis[0*numnodes+j]*dbasis[1*numnodes+i] + dbasis[1*numnodes+j]*dbasis[0*numnodes+i]
							)*(n+1)/(n+2);//K32
				Ke->values[(4*i+2)*2*2*numnodes+4*j+2] += factor*viscosity[0]*(
							4.*dbasis[1*numnodes+j]*dbasis[1*numnodes+i] + dbasis[0*numnodes+j]*dbasis[0*numnodes+i]
							);//K33
				Ke->values[(4*i+2)*2*2*numnodes+4*j+3] += factor*viscosity[1]*(
							4.*dbasis[1*numnodes+j]*dbasis[1*numnodes+i] + dbasis[0*numnodes+j]*dbasis[0*numnodes+i]
							)*(n+1)/(n+2);//K34

				Ke->values[(4*i+3)*2*2*numnodes+4*j+0] += factor*viscosity[1]*(
                     2.*dbasis[0*numnodes+j]*dbasis[1*numnodes+i] + dbasis[1*numnodes+j]*dbasis[0*numnodes+i]
                     )*(n+1)/(n+2);//K41
				Ke->values[(4*i+3)*2*2*numnodes+4*j+1] += factor*viscosity[2]*(
                     2.*dbasis[0*numnodes+j]*dbasis[1*numnodes+i] + dbasis[1*numnodes+j]*dbasis[0*numnodes+i]
                     )*2*pow(n+1,2)/((2*n+3)*(n+2));//K42
				Ke->values[(4*i+3)*2*2*numnodes+4*j+2] += factor*viscosity[1]*(
							4.*dbasis[1*numnodes+j]*dbasis[1*numnodes+i] + dbasis[0*numnodes+j]*dbasis[0*numnodes+i]
							)*(n+1)/(n+2);//K43
				Ke->values[(4*i+3)*2*2*numnodes+4*j+3] += factor*viscosity[2]*(
							4.*dbasis[1*numnodes+j]*dbasis[1*numnodes+i] + dbasis[0*numnodes+j]*dbasis[0*numnodes+i]
							)*2*pow(n+1,2)/((2*n+3)*(n+2))
							+ 
							gauss->weight*Jdet*viscosity[3]*basis[j]*basis[i]*pow(n+1,2)/(thickness*(2*n+1))
							;//K44
			}
		}
	}

	/*Transform Coordinate System*/
	//basalelement->TransformStiffnessMatrixCoord(Ke,XYMOLHOEnum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	xDelete<IssmDouble>(basis);
	return Ke;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorMOLHO(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int      domaintype;
	Element* basalelement;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum: case Domain2DverticalEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*compute all load vectors for this element*/
	ElementVector* pe1=CreatePVectorMOLHODrivingStress(basalelement);
	ElementVector* pe2=CreatePVectorMOLHOFront(basalelement);
	ElementVector* pe =new ElementVector(pe1,pe2);

	/*clean-up and return*/
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	delete pe1;
	delete pe2;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorMOLHODrivingStress(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble  thickness,Jdet,slope[2],n;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and vectors*/
	ElementVector* pe    = element->NewElementVector(MOLHOApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input*     thickness_input=element->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input*     surface_input  =element->GetInput(SurfaceEnum);   _assert_(surface_input);
	Input*     n_input        =element->GetInput(MaterialsRheologyNEnum); _assert_(n_input);
	IssmDouble rhog = element->FindParam(MaterialsRhoIceEnum)*element->FindParam(ConstantsGEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis, gauss);

		thickness_input->GetInputValue(&thickness,gauss);
		surface_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);
		n_input->GetInputValue(&n,gauss);

		IssmDouble factor = rhog*thickness*Jdet*gauss->weight;
		for(int i=0;i<numnodes;i++){// per node: vx (basal vx), vshx, vy (basal vy), vshy
			pe->values[i*4+0]+=-factor*slope[0]*basis[i]; //F1
			pe->values[i*4+1]+=-factor*slope[0]*basis[i]*(n+1)/(n+2); //F2
			pe->values[i*4+2]+=-factor*slope[1]*basis[i]; //F3
			pe->values[i*4+3]+=-factor*slope[1]*basis[i]*(n+1)/(n+2); //F4
		}
	}

	/*Transform coordinate system*/
	//element->TransformLoadVectorCoord(pe,XYMOLHOEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorMOLHOFront(Element* element){/*{{{*/

	/*If no front, return NULL*/
	if(!element->IsIcefront()) return NULL;

	/*Intermediaries*/
	IssmDouble  Jdet,thickness,base,sealevel,water_pressure,ice_pressure;
	IssmDouble  water_pressure_sh,ice_pressure_sh,pressure_sh,pressure,n,s,b;
	IssmDouble *xyz_list = NULL;
	IssmDouble *xyz_list_front = NULL;
	IssmDouble  normal[2];

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector(MOLHOApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	Input* thickness_input = element->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* base_input      = element->GetInput(BaseEnum);       _assert_(base_input);
	Input* sealevel_input  = element->GetInput(SealevelEnum);       _assert_(sealevel_input);
	Input* n_input         = element->GetInput(MaterialsRheologyNEnum); _assert_(n_input);
	IssmDouble rho_water   = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice     = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble gravity     = element->FindParam(ConstantsGEnum);
	element->GetVerticesCoordinates(&xyz_list);
	element->GetIcefrontCoordinates(&xyz_list_front,xyz_list,MaskIceLevelsetEnum);
	element->NormalSection(&normal[0],xyz_list_front);

	/*Start looping on Gaussian points*/
	Gauss* gauss=element->NewGauss(xyz_list,xyz_list_front,3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		base_input->GetInputValue(&base,gauss);
		n_input->GetInputValue(&n,gauss);
		sealevel_input->GetInputValue(&sealevel,gauss);
		element->JacobianDeterminantSurface(&Jdet,xyz_list_front,gauss);
		element->NodalFunctions(basis,gauss);

		b = base-sealevel; // ice base shifted by the sea level
		s = thickness+b;   // ice surface shifted by the sea level
		/*Vertically integrated pressure - SSA type*/
		water_pressure = -(1.0/2.0)*gravity*rho_water*(b*b);
		ice_pressure   =  (1.0/2.0)*gravity*rho_ice*thickness*thickness;
		pressure = ice_pressure + water_pressure;
		/*Vertically integrated pressure - HO type*/
		water_pressure_sh = gravity*rho_water*( -b*b/2 - (s*thickness/(n+2))*(1-pow(s/thickness,n+2)) + (thickness*thickness/(n+3))*(1-pow(s/thickness,n+3)) );
		ice_pressure_sh   = gravity*rho_ice*thickness*thickness*(n+1)/(2*(n+3));
		pressure_sh = ice_pressure_sh + water_pressure_sh;
		if (b>0) {
			pressure = ice_pressure;
			pressure_sh = ice_pressure_sh;
		}

		IssmDouble factor = Jdet*gauss->weight*pressure;
		IssmDouble factor_sh = Jdet*gauss->weight*pressure_sh;
		for (int i=0;i<numnodes;i++){
			pe->values[i*4+0]+= factor*normal[0]*basis[i]; // F1
			pe->values[i*4+1]+= factor_sh*normal[0]*basis[i]; // F2
			pe->values[i*4+2]+= factor*normal[1]*basis[i]; // F3
			pe->values[i*4+3]+= factor_sh*normal[1]*basis[i]; // F4
		}
	}

	/*Transform coordinate system*/
	//element->TransformLoadVectorCoord(pe,XYMOLHOEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_front);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolutionMOLHO(IssmDouble* solution,Element* element){/*{{{*/

	int         i,dim,domaintype;
	IssmDouble  rho_ice,g;
	int*        doflist=NULL;
	IssmDouble* xyz_list=NULL;
	Element*    basalelement=NULL;

	/*Deal with pressure first*/
	int numvertices = element->GetNumberOfVertices();
	IssmDouble* pressure  = xNew<IssmDouble>(numvertices);
	IssmDouble* thickness = xNew<IssmDouble>(numvertices);
	IssmDouble* surface   = xNew<IssmDouble>(numvertices);

	element->FindParam(&domaintype,DomainTypeEnum);
	rho_ice =element->FindParam(MaterialsRhoIceEnum);
	g       =element->FindParam(ConstantsGEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->GetInputListOnVertices(thickness,ThicknessEnum);
			for(i=0;i<numvertices;i++) pressure[i]=rho_ice*g*thickness[i];
			dim=2;
			break;
		case Domain3DEnum:
			element->GetVerticesCoordinates(&xyz_list);
			element->GetInputListOnVertices(surface,SurfaceEnum);
			for(i=0;i<numvertices;i++) pressure[i]=rho_ice*g*(surface[i]-xyz_list[i*3+2]);
			dim=2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	element->AddInput(PressureEnum,pressure,P1Enum);
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(thickness);
	xDelete<IssmDouble>(surface);

	/*Get basal element*/
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()){xDelete<IssmDouble>(xyz_list); return;}
			basalelement=element->SpawnBasalElement();
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();
	int numdof   = numnodes*dim*2; //2xdim DOFs per node

	/*Fetch dof list and allocate solution vectors*/
	basalelement->GetDofListLocal(&doflist,MOLHOApproximationEnum,GsetEnum); 
	IssmDouble* values    = xNew<IssmDouble>(numdof);
	IssmDouble* vbx       = xNew<IssmDouble>(numnodes);
	IssmDouble* vby       = xNew<IssmDouble>(numnodes);
   IssmDouble* vshx      = xNew<IssmDouble>(numnodes); 
   IssmDouble* vshy      = xNew<IssmDouble>(numnodes);
	IssmDouble* vsx       = xNew<IssmDouble>(numnodes);
	IssmDouble* vsy       = xNew<IssmDouble>(numnodes);
	IssmDouble* vx        = xNew<IssmDouble>(numnodes);
	IssmDouble* vy        = xNew<IssmDouble>(numnodes);
	IssmDouble* vz        = xNew<IssmDouble>(numnodes);
	IssmDouble* vel       = xNew<IssmDouble>(numnodes);
	IssmDouble* n			 = xNew<IssmDouble>(numnodes);
	IssmDouble* H			 = xNew<IssmDouble>(numnodes);
	IssmDouble* s			 = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numdof;i++) values[i]=solution[doflist[i]];

	/*Transform solution in Cartesian Space*/
	if(dim==2) basalelement->TransformSolutionCoord(&values[0],XYEnum);

   /*Ok, we have vx and vy in values, fill in vx and vy arrays: */
   for(i=0;i<numnodes;i++){ //numnodes of the 2D mesh in which the MOLHO is written
      vbx[i] =values[i*4+0]; //base vx
      vshx[i]=values[i*4+1]; //shear vx
		vsx[i] =vbx[i]+vshx[i]; //surface vx
		if(xIsNan<IssmDouble>(vbx[i]))	_error_("NaN found in solution vector");
      if(xIsInf<IssmDouble>(vbx[i]))	_error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(vshx[i]))	_error_("NaN found in solution vector");
      if(xIsInf<IssmDouble>(vshx[i]))	_error_("Inf found in solution vector");
		vby[i] =values[i*4+2]; //base vy
		vshy[i]=values[i*4+3]; //shear vy
		vsy[i] =vby[i]+vshy[i]; //surface vy
		if(xIsNan<IssmDouble>(vby[i]))	_error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vby[i]))	_error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(vshy[i]))	_error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vshy[i]))	_error_("Inf found in solution vector");
	}

	/*Add vx and vy as inputs to the tria element (shear velocities): */
	element->AddBasalInput(VxShearEnum,vshx,element->GetElementType());
	element->AddBasalInput(VyShearEnum,vshy,element->GetElementType());

	/*Add vx and vy as inputs to the tria element (base velocities): */
	element->AddBasalInput(VxBaseEnum,vbx,element->GetElementType());
	element->AddBasalInput(VyBaseEnum,vby,element->GetElementType());

	/*Add vx and vy as inputs to the tria element (surface velocities): */
	element->AddBasalInput(VxSurfaceEnum,vsx,element->GetElementType());
	element->AddBasalInput(VySurfaceEnum,vsy,element->GetElementType());

	/*Compute the vertically averaged velocities on each node*/
	basalelement->GetInputListOnNodes(&n[0],MaterialsRheologyNEnum,0.); 

	/* Reconstruct vx, vy and vz solutions for 3D problem
	 * Add vx and vy as inputs to the tria element (vertically averaged velocities): */

   switch(domaintype){
      case Domain2DhorizontalEnum:
			for(i=0;i<numnodes;i++){ //numnodes of the 2D mesh in which the MOLHO is written
				vx[i]=vbx[i]+vshx[i]*(n[i]+1)/(n[i]+2);
				vy[i]=vby[i]+vshy[i]*(n[i]+1)/(n[i]+2);
				vel[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i]); 
			}
			element->AddBasalInput(VxEnum,vx,element->GetElementType());
			element->AddBasalInput(VyEnum,vy,element->GetElementType());
			element->AddBasalInput(VelEnum,vel,element->GetElementType()); 
         break;
      case Domain3DEnum:
		   basalelement->GetInputListOnNodes(&H[0],ThicknessEnum,0.);
		   basalelement->GetInputListOnNodes(&s[0],SurfaceEnum,0.);
			element->Recover3DMOLHOInput(VxEnum, numnodes, vbx, vshx, n, H, s);
			element->Recover3DMOLHOInput(VyEnum, numnodes, vby, vshy, n, H, s);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Free resources:*/
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vz);
	xDelete<IssmDouble>(vby);
	xDelete<IssmDouble>(vbx);
	xDelete<IssmDouble>(vsy);
	xDelete<IssmDouble>(vsx);
	xDelete<IssmDouble>(vshy);
	xDelete<IssmDouble>(vshx);
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(n);
	xDelete<IssmDouble>(s);
	xDelete<IssmDouble>(H);
	xDelete<int>(doflist);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           StressbalanceAnalysis::GetSolutionFromInputsMOLHO(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	IssmDouble   vbx,vby,vshx,vshy;
	int          domaintype,dim,approximation,dofpernode;
	int*         doflist = NULL;

	/*Get some parameters*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum: dim = 2; dofpernode = 4; break;
		case Domain3DEnum: dim = 2; dofpernode = 4; break;
		case Domain2DverticalEnum: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet"); break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*dofpernode;
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation!=MOLHOApproximationEnum) _error_("mesh "<<EnumToStringx(approximation)<<" not supported here");

	/*Fetch dof list and allocate solution vector*/
	element->GetDofList(&doflist,approximation,GsetEnum); 
	IssmDouble* values = xNew<IssmDouble>(numdof);

	/*Get inputs*/
	Input* vxbase_input	=element->GetInput(VxBaseEnum);	_assert_(vxbase_input);
	Input* vxshear_input	=element->GetInput(VxShearEnum);	_assert_(vxshear_input);
	Input* vybase_input	=element->GetInput(VyBaseEnum);	_assert_(vybase_input);
	Input* vyshear_input	=element->GetInput(VyShearEnum);	_assert_(vyshear_input);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	Gauss* gauss=element->NewGauss();
	for(int i=0;i<numnodes;i++){
		gauss->GaussNode(element->FiniteElement(),i);

		/*Recover vx and vy*/
		vxbase_input->GetInputValue(&vbx,gauss);	//base vx
		vxshear_input->GetInputValue(&vshx,gauss);//shear vx
		values[i*dofpernode+0]=vbx;  //base vx
		values[i*dofpernode+1]=vshx; //shear vx
		vybase_input->GetInputValue(&vby,gauss);	//base vy
		vyshear_input->GetInputValue(&vshy,gauss);//shear vy
		values[i*dofpernode+2]=vby; //base vy
		values[i*dofpernode+3]=vshy;//shear vy  
	}

	solution->SetValues(numdof,doflist,values,INS_VAL);

	/*Free resources:*/
	delete gauss;
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);
}/*}}}*/

/*HO*/
ElementMatrix* StressbalanceAnalysis::CreateJacobianMatrixHO(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble Jdet;
	IssmDouble eps1dotdphii,eps1dotdphij;
	IssmDouble eps2dotdphii,eps2dotdphij;
	IssmDouble mu_prime;
	IssmDouble epsilon[5]; /* epsilon=[exx,eyy,exy,exz,eyz];*/
	IssmDouble eps1[3],eps2[3];
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element matrix, vectors and Gaussian points*/
	ElementMatrix* Ke=this->CreateKMatrixHO(element); //Initialize Jacobian with regular HO (first part of the Gateau derivative)
	IssmDouble*    dbasis = xNew<IssmDouble>(3*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vx_input = element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = element->GetInput(VyEnum); _assert_(vy_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		element->StrainRateHO(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
		element->material->ViscosityHODerivativeEpsSquare(&mu_prime,&epsilon[0],gauss);
		eps1[0]=2*epsilon[0]+epsilon[1];   eps2[0]=epsilon[2];
		eps1[1]=epsilon[2];                eps2[1]=epsilon[0]+2*epsilon[1];
		eps1[2]=epsilon[3];                eps2[2]=epsilon[4];

		IssmDouble factor = gauss->weight*Jdet*2.*mu_prime;
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				eps1dotdphii=eps1[0]*dbasis[0*numnodes+i]+eps1[1]*dbasis[1*numnodes+i]+eps1[2]*dbasis[2*numnodes+i];
				eps1dotdphij=eps1[0]*dbasis[0*numnodes+j]+eps1[1]*dbasis[1*numnodes+j]+eps1[2]*dbasis[2*numnodes+j];
				eps2dotdphii=eps2[0]*dbasis[0*numnodes+i]+eps2[1]*dbasis[1*numnodes+i]+eps2[2]*dbasis[2*numnodes+i];
				eps2dotdphij=eps2[0]*dbasis[0*numnodes+j]+eps2[1]*dbasis[1*numnodes+j]+eps2[2]*dbasis[2*numnodes+j];

				Ke->values[2*numnodes*(2*i+0)+2*j+0]+=factor*eps1dotdphij*eps1dotdphii;
				Ke->values[2*numnodes*(2*i+0)+2*j+1]+=factor*eps2dotdphij*eps1dotdphii;
				Ke->values[2*numnodes*(2*i+1)+2*j+0]+=factor*eps1dotdphij*eps2dotdphii;
				Ke->values[2*numnodes*(2*i+1)+2*j+1]+=factor*eps2dotdphij*eps2dotdphii;
			}
		}
	}

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixHO(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixHOViscous(element);
	ElementMatrix* Ke2=CreateKMatrixHOFriction(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixHOFriction(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	if(element->IsAllFloating() || !element->IsOnBase()) return NULL;

	/*Intermediaries*/
	int         dim;
	bool        mainlyfloating;
	int         friction_style,point1;
	IssmDouble  alpha2,Jdet,fraction1,fraction2;
	IssmDouble  gllevelset,phi=1.;
	IssmDouble *xyz_list_base = NULL;
	Gauss*      gauss         = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*(dim-1);

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix(HOApproximationEnum);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->FindParam(&friction_style,GroundinglineFrictionInterpolationEnum);
	Input* gllevelset_input = NULL;

	/*build friction object, used later on: */
	/*dim=4 is special for HO, which is actually 2.5D*/
	Friction* friction=NULL;
	if(dim==3)
	 friction = new Friction(element, 2.5); 
	else
	 friction = new Friction(element, 1); 

	/*Recover portion of element that is grounded*/
	if(!(friction_style==SubelementFriction2Enum)) phi=element->GetGroundedPortion(xyz_list_base);
	if(friction_style==SubelementFriction2Enum){
		gllevelset_input=element->GetInput(MaskOceanLevelsetEnum); _assert_(gllevelset_input);
		element->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating, MaskOceanLevelsetEnum,0);
		gauss = element->NewGauss(point1,fraction1,fraction2,mainlyfloating,2);
	}
	else{
		gauss=element->NewGaussBase(2);
	}

	/* Start  looping on the number of gaussian points: */
	while(gauss->next()){

		friction->GetAlpha2(&alpha2,gauss);
		if(friction_style==SubelementFriction1Enum) alpha2=phi*alpha2;
		else if(friction_style==SubelementFriction2Enum){
			gllevelset_input->GetInputValue(&gllevelset, gauss);
			if(gllevelset<0.) alpha2=0.;
		}
		else if(friction_style==NoFrictionOnPartiallyFloatingEnum){
			if (phi<0.99999999) alpha2=0.;
		}
		else  _error_("friction interpolation "<<EnumToStringx(friction_style)<<" not implemented yet");

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctions(basis,gauss);

		IssmDouble factor = alpha2*gauss->weight*Jdet;
		if(dim==3){
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[2*i*2*numnodes+2*j]       += factor*basis[i]*basis[j];
					Ke->values[(2*i+1)*2*numnodes+2*j+1] += factor*basis[i]*basis[j];
				}
			}
		}
		else{
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += factor*basis[i]*basis[j];
				}
			}
		}
	}

	/*Transform Coordinate System*/
	if(dim==3) element->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	delete gauss;
	delete friction;
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(basis);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixHOViscous(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int         dim,bsize;
	IssmDouble  viscosity,Jdet;
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*(dim-1);

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix(HOApproximationEnum);
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vx_input    = element->GetInput(VxEnum);       _assert_(vx_input);
	Input* vy_input    = NULL;
	if(dim==3){
		vy_input=element->GetInput(VyEnum);          _assert_(vy_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		element->material->ViscosityHO(&viscosity,dim,xyz_list,gauss,vx_input,vy_input);

		IssmDouble factor =  gauss->weight*Jdet*viscosity;
		if(dim==3){
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[2*i*2*numnodes+2*j] += factor*(
								4.*dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i] + dbasis[2*numnodes+j]*dbasis[2*numnodes+i]
								);
					Ke->values[2*i*2*numnodes+2*j+1] += factor*(
								2.*dbasis[1*numnodes+j]*dbasis[0*numnodes+i] + dbasis[0*numnodes+j]*dbasis[1*numnodes+i]
								);
					Ke->values[(2*i+1)*2*numnodes+2*j] += factor*(
								2.*dbasis[0*numnodes+j]*dbasis[1*numnodes+i] + dbasis[1*numnodes+j]*dbasis[0*numnodes+i]
								);
					Ke->values[(2*i+1)*2*numnodes+2*j+1] += factor*(
								dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + 4.*dbasis[1*numnodes+j]*dbasis[1*numnodes+i] + dbasis[2*numnodes+j]*dbasis[2*numnodes+i]
								);
				}
			}
		}
		else{
			for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += factor*(
								4.*dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]
								);
				}
			}
		}
	}

	/*Transform Coordinate System*/
	if(dim==3) element->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	return Ke;
}/*}}}*/
#ifdef FSANALYTICAL
ElementVector* StressbalanceAnalysis::CreatePVectorHO(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         dim;
	IssmDouble  x_coord,y_coord,z_coord;
	IssmDouble  Jdet,forcex,forcey,forcez;
	IssmDouble* xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and vectors*/
	ElementVector* pe=element->NewElementVector(HOApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(3);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis, gauss);

		x_coord=element->GetXcoord(xyz_list,gauss);
		y_coord=element->GetYcoord(xyz_list,gauss);
		if(dim==3) z_coord=element->GetZcoord(xyz_list,gauss);
		else z_coord=0.;

		forcex=fx(x_coord,y_coord,z_coord,FSANALYTICAL);
		forcey=fy(x_coord,y_coord,z_coord,FSANALYTICAL);

		IssmDouble Jdet*gauss->weight;
		for(int i=0;i<numnodes;i++){
			pe->values[i*(dim-1)+0]+=forcex*factor*basis[i];
			pe->values[i*(dim-1)+1]+=forcey*factor*basis[i];
		}
	}

	/*Transform coordinate system*/
	if(dim==3) element->TransformLoadVectorCoord(pe,XYEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return pe;
}/*}}}*/
#else
ElementVector* StressbalanceAnalysis::CreatePVectorHO(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*compute all load vectors for this element*/
	ElementVector* pe1=CreatePVectorHODrivingStress(element);
	ElementVector* pe2=CreatePVectorHOFront(element);
	ElementVector* pe =new ElementVector(pe1,pe2);

	/*clean-up and return*/
	delete pe1;
	delete pe2;
	return pe;
}/*}}}*/
#endif
ElementVector* StressbalanceAnalysis::CreatePVectorHODrivingStress(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         dim;
	IssmDouble  Jdet,slope[3];
	IssmDouble* xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and vectors*/
	ElementVector* pe=element->NewElementVector(HOApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input*     surface_input = element->GetInput(SurfaceEnum);   _assert_(surface_input);
	IssmDouble rhog = element->FindParam(MaterialsRhoIceEnum)*element->FindParam(ConstantsGEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(3);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis, gauss);
		surface_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);

		IssmDouble factor = -rhog*Jdet*gauss->weight;
		for(int i=0;i<numnodes;i++){
			pe->values[i*(dim-1)+0]+=factor*slope[0]*basis[i];
			if(dim==3) pe->values[i*(dim-1)+1]+=factor*slope[1]*basis[i];
		}
	}

	/*Transform coordinate system*/
	if(dim==3) element->TransformLoadVectorCoord(pe,XYEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorHOFront(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*If no front, return NULL*/
	if(!element->IsIcefront()) return NULL;

	/*Intermediaries*/
	int         dim;
	IssmDouble  Jdet,surface,sealevel,z,water_pressure,ice_pressure;
	IssmDouble  surface_under_water,base_under_water,pressure;
	IssmDouble* xyz_list       = NULL;
	IssmDouble* xyz_list_front = NULL;
	IssmDouble  normal[3];
	Gauss*      gauss = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = element->GetNumberOfNodes();
	int numvertices = element->GetNumberOfVertices();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector(HOApproximationEnum);
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	Input* surface_input = element->GetInput(SurfaceEnum); _assert_(surface_input);
	Input* sealevel_input       = element->GetInput(SealevelEnum);       _assert_(sealevel_input);
	IssmDouble rho_water = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble gravity   = element->FindParam(ConstantsGEnum);
	element->GetVerticesCoordinates(&xyz_list);
	element->GetIcefrontCoordinates(&xyz_list_front,xyz_list,MaskIceLevelsetEnum);
	element->NormalSection(&normal[0],xyz_list_front);

	/*Initialize gauss points*/
	IssmDouble zmax=xyz_list[0*3+(dim-1)]; for(int i=1;i<numvertices;i++) if(xyz_list[i*3+(dim-1)]>zmax) zmax=xyz_list[i*3+(dim-1)];
	IssmDouble zmin=xyz_list[0*3+(dim-1)]; for(int i=1;i<numvertices;i++) if(xyz_list[i*3+(dim-1)]<zmin) zmin=xyz_list[i*3+(dim-1)];
	if(zmax>0. && zmin<0.) gauss=element->NewGauss(xyz_list,xyz_list_front,3,10);//refined in vertical because of the sea level discontinuity
	else                   gauss=element->NewGauss(xyz_list,xyz_list_front,3,3);

	/* Start  looping on the number of gaussian points: */
	while(gauss->next()){
		surface_input->GetInputValue(&surface,gauss);
		sealevel_input->GetInputValue(&sealevel,gauss);
		if(dim==3) z=element->GetZcoord(xyz_list,gauss);
		else       z=element->GetYcoord(xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->JacobianDeterminantSurface(&Jdet,xyz_list_front,gauss);

		water_pressure = rho_water*gravity*min(0.,z-sealevel);//0 if the gaussian point is above water level
		ice_pressure   = rho_ice*gravity*(surface-z);
		pressure       = ice_pressure + water_pressure;

		IssmDouble factor = pressure*Jdet*gauss->weight;
		for (int i=0;i<numnodes;i++){
			pe->values[(dim-1)*i+0]+= factor*normal[0]*basis[i];
			if(dim==3) pe->values[(dim-1)*i+1]+= factor*normal[1]*basis[i];
		}
	}

	/*Transform coordinate system*/
	if(dim==3)element->TransformLoadVectorCoord(pe,XYEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_front);
	delete gauss;
	return pe;
}/*}}}*/
void           StressbalanceAnalysis::GetBHOFriction(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3] where Bi is square and of size 2.
	 * For node i, Bi can be expressed in the actual coordinate system
	 * by:
	 *                       3D           2D
	 *                 Bi=[ N   0 ]    Bi=N
	 *                    [ 0   N ]
	 * where N is the finiteelement function for node i.
	 *
	 * We assume B has been allocated already, of size: 2 x (numdof*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* basis=xNew<IssmDouble>(numnodes);
	element->NodalFunctions(basis,gauss);

	/*Build L: */
	if(dim==3){
		for(int i=0;i<numnodes;i++){
			B[2*numnodes*0+2*i+0] = basis[i];
			B[2*numnodes*0+2*i+1] = 0.;
			B[2*numnodes*1+2*i+0] = 0.;
			B[2*numnodes*1+2*i+1] = basis[i];
		}
	}
	else{
		for(int i=0;i<numnodes;i++){
			B[i] = basis[i];
		}
	}

	/*Clean-up*/
	xDelete<IssmDouble>(basis);
}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolutionHO(IssmDouble* solution,Element* element){/*{{{*/

	int         i,domaintype,dim;
	int*        doflist=NULL;
	IssmDouble* xyz_list=NULL;

	/*Deal with pressure first*/
	int numvertices = element->GetNumberOfVertices();
	IssmDouble* pressure  = xNew<IssmDouble>(numvertices);
	IssmDouble* surface   = xNew<IssmDouble>(numvertices);

	element->FindParam(&domaintype,DomainTypeEnum);
	IssmDouble rho_ice =element->FindParam(MaterialsRhoIceEnum);
	IssmDouble g       =element->FindParam(ConstantsGEnum);
	switch(domaintype){
		case Domain3DEnum:
			element->GetVerticesCoordinates(&xyz_list);
			element->GetInputListOnVertices(surface,SurfaceEnum);
			for(i=0;i<numvertices;i++) pressure[i]=rho_ice*g*(surface[i]-xyz_list[i*3+2]);
			dim=3;
			break;
		case Domain2DverticalEnum:
			element->GetVerticesCoordinates(&xyz_list);
			element->GetInputListOnVertices(surface,SurfaceEnum);
			for(i=0;i<numvertices;i++) pressure[i]=rho_ice*g*(surface[i]-xyz_list[i*3+1]);
			dim=2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	element->AddInput(PressureEnum,pressure,P1Enum);
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(surface);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*(dim-1);

	/*Fetch dof list and allocate solution vectors*/
	element->GetDofListLocal(&doflist,HOApproximationEnum,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numdof);
	IssmDouble* vx     = xNew<IssmDouble>(numnodes);
	IssmDouble* vy     = xNew<IssmDouble>(numnodes);
	IssmDouble* vz     = xNew<IssmDouble>(numnodes);
	IssmDouble* vel    = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numdof;i++) values[i]=solution[doflist[i]];

	/*Transform solution in Cartesian Space*/
	if(dim==3) element->TransformSolutionCoord(&values[0],XYEnum);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	for(i=0;i<numnodes;i++){
		vx[i]=values[i*(dim-1)+0];
		if(xIsNan<IssmDouble>(vx[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vx[i])) _error_("Inf found in solution vector");
		if(dim==3){
			vy[i]=values[i*(dim-1)+1];
			if(xIsNan<IssmDouble>(vy[i])) _error_("NaN found in solution vector");
			if(xIsInf<IssmDouble>(vy[i])) _error_("Inf found in solution vector");
		}
	}

	/*Get Vz and compute vel*/
	if(dim==3){
		element->GetInputListOnNodes(&vz[0],VzEnum,0.);
		for(i=0;i<numnodes;i++) vel[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	}
	else{
		element->GetInputListOnNodes(&vy[0],VyEnum,0.);
		for(i=0;i<numnodes;i++) vel[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
	}

	/*Add vx and vy as inputs to the element: */
	element->AddInput(VxEnum,vx,element->GetElementType());
	if(dim==3)element->AddInput(VyEnum,vy,element->GetElementType());
	element->AddInput(VelEnum,vel,element->GetElementType());

	/*Free resources:*/
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vz);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(values);
	xDelete<IssmDouble>(xyz_list);
	xDelete<int>(doflist);
}/*}}}*/

/*FS*/
ElementVector* StressbalanceAnalysis::CreateDVectorFS(Element* element){/*{{{*/

	int dim;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();

	/*Initialize output vector*/
	ElementVector* de = element->NewElementVector(FSvelocityEnum);

	for(int i=0;i<vnumnodes;i++){
		for(int j=0;j<dim;j++) de->values[i*dim+j]=VelocityEnum;
	}
	for(int i=0;i<pnumnodes;i++){
		de->values[vnumnodes*dim+i]=PressureEnum;
	}

	return de;

}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateJacobianMatrixFS(Element* element){/*{{{*/

	/*Intermediaries */
	int        i,j;
	IssmDouble Jdet;
	IssmDouble eps1dotdphii,eps1dotdphij;
	IssmDouble eps2dotdphii,eps2dotdphij;
	IssmDouble eps3dotdphii,eps3dotdphij;
	IssmDouble mu_prime;
	IssmDouble epsilon[6]; /* epsilon=[exx,eyy,ezz,exy,exz,eyz];*/
	IssmDouble eps1[3],eps2[3],eps3[3];
	IssmDouble *xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*3 + pnumnodes*1;

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	for(i=0;i<vnumnodes;i++) cs_list[i]           = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Initialize Element matrix, vectors and Gaussian points*/
	ElementMatrix* Ke=this->CreateKMatrixFS(element); //Initialize Jacobian with regular FS (first part of the Gateau derivative)
	IssmDouble*    dbasis = xNew<IssmDouble>(3*vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vx_input = element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = element->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input = element->GetInput(VzEnum); _assert_(vz_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivativesVelocity(dbasis,xyz_list,gauss);

		//element->StrainRateFS(&epsilon[0],xyz_list,gauss,vx_input,vy_input,vz_input);
		//eps1[0]=epsilon[0];   eps2[0]=epsilon[3];   eps3[0]=epsilon[4];
		//eps1[1]=epsilon[3];   eps2[1]=epsilon[1];   eps3[1]=epsilon[5];
		//eps1[2]=epsilon[4];   eps2[2]=epsilon[5];   eps3[2]=epsilon[2];
		element->StrainRateHO(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
		eps1[0]=epsilon[0];   eps2[0]=epsilon[2];   eps3[0]=epsilon[3];
		eps1[1]=epsilon[2];   eps2[1]=epsilon[1];   eps3[1]=epsilon[4];
		eps1[2]=epsilon[3];   eps2[2]=epsilon[4];   eps3[2]= -epsilon[0] -epsilon[1];
		element->material->ViscosityFSDerivativeEpsSquare(&mu_prime,&epsilon[0],gauss);

		IssmDouble factor = gauss->weight*Jdet*2*mu_prime;
		for(i=0;i<vnumnodes;i++){
			for(j=0;j<vnumnodes;j++){
				eps1dotdphii=eps1[0]*dbasis[0*vnumnodes+i]+eps1[1]*dbasis[1*vnumnodes+i]+eps1[2]*dbasis[2*vnumnodes+i];
				eps1dotdphij=eps1[0]*dbasis[0*vnumnodes+j]+eps1[1]*dbasis[1*vnumnodes+j]+eps1[2]*dbasis[2*vnumnodes+j];
				eps2dotdphii=eps2[0]*dbasis[0*vnumnodes+i]+eps2[1]*dbasis[1*vnumnodes+i]+eps2[2]*dbasis[2*vnumnodes+i];
				eps2dotdphij=eps2[0]*dbasis[0*vnumnodes+j]+eps2[1]*dbasis[1*vnumnodes+j]+eps2[2]*dbasis[2*vnumnodes+j];
				eps3dotdphii=eps3[0]*dbasis[0*vnumnodes+i]+eps3[1]*dbasis[1*vnumnodes+i]+eps3[2]*dbasis[2*vnumnodes+i];
				eps3dotdphij=eps3[0]*dbasis[0*vnumnodes+j]+eps3[1]*dbasis[1*vnumnodes+j]+eps3[2]*dbasis[2*vnumnodes+j];

				Ke->values[numdof*(3*i+0)+3*j+0]+=factor*eps1dotdphij*eps1dotdphii;
				Ke->values[numdof*(3*i+0)+3*j+1]+=factor*eps2dotdphij*eps1dotdphii;
				Ke->values[numdof*(3*i+0)+3*j+2]+=factor*eps3dotdphij*eps1dotdphii;

				Ke->values[numdof*(3*i+1)+3*j+0]+=factor*eps1dotdphij*eps2dotdphii;
				Ke->values[numdof*(3*i+1)+3*j+1]+=factor*eps2dotdphij*eps2dotdphii;
				Ke->values[numdof*(3*i+1)+3*j+2]+=factor*eps3dotdphij*eps2dotdphii;

				Ke->values[numdof*(3*i+2)+3*j+0]+=factor*eps1dotdphij*eps3dotdphii;
				Ke->values[numdof*(3*i+2)+3*j+1]+=factor*eps2dotdphij*eps3dotdphii;
				Ke->values[numdof*(3*i+2)+3*j+2]+=factor*eps3dotdphij*eps3dotdphii;
			}
		}
	}

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	xDelete<int>(cs_list);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixFS(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Get type of algorithm*/
	int fe_FS;
	bool isNitsche;

	element->FindParam(&fe_FS,FlowequationFeFSEnum);
	element->FindParam(&isNitsche,FlowequationIsNitscheEnum);

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=NULL;
	if(fe_FS==XTaylorHoodEnum)
	 Ke1=CreateKMatrixFSViscousXTH(element);
	else if(fe_FS==LATaylorHoodEnum || fe_FS==LACrouzeixRaviartEnum)
	 Ke1=CreateKMatrixFSViscousLA(element);
	else if(fe_FS==P1P1GLSEnum)
	 Ke1=CreateKMatrixFSViscousGLS(element);
	else
	 Ke1=CreateKMatrixFSViscous(element);

	ElementMatrix* Ke2;
	if (isNitsche) {
		Ke2 = CreateKMatrixFSFrictionNitsche(element);
	}
	else {
		Ke2 = CreateKMatrixFSFriction(element);
	}
	ElementMatrix* Ke3=CreateKMatrixFSShelf(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2,Ke3);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	delete Ke3;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixFSShelf(Element* element){/*{{{*/

	if(!element->IsAllFloating() || !element->IsOnBase()) return NULL;

	/*If on not water or not FS, skip stiffness: */
	int approximation,shelf_dampening;
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation!=FSApproximationEnum && approximation!=SSAFSApproximationEnum && approximation!=HOFSApproximationEnum) return NULL;
	element->FindParam(&shelf_dampening,StressbalanceShelfDampeningEnum);
	if(shelf_dampening==0) return NULL;

	/*Intermediaries*/
	bool        mainlyfloating;
	int         j,i,dim;
	IssmDouble  Jdet,slope2,scalar,dt;
	IssmDouble  slope[3];
	IssmDouble *xyz_list_base = NULL;
	IssmDouble *xyz_list      = NULL;
	Gauss*      gauss         = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke = element->NewElementMatrix(FSvelocityEnum);
	IssmDouble*    vbasis = xNew<IssmDouble>(vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	if(dt==0)   dt=1.e+5;
	IssmDouble  rho_water     = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble  gravity       = element->FindParam(ConstantsGEnum);
	Input*      base_input = element->GetInput(BaseEnum); _assert_(base_input);

	/* Start  looping on the number of gaussian points: */
	gauss=element->NewGaussBase(3);
	while(gauss->next()){

		base_input->GetInputDerivativeValue(&slope[0],xyz_list,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);
		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		if(dim==2) slope2=slope[0]*slope[0];
		else if(dim==3) slope2=slope[0]*slope[0]+slope[1]*slope[1];
		scalar  = rho_water*gravity*sqrt(1+slope2)*gauss->weight*Jdet*dt;
		for(i=0;i<vnumnodes;i++){
			for(j=0;j<vnumnodes;j++){
				Ke->values[numdof*((i+1)*dim-1)+(j+1)*dim-1] += scalar*vbasis[i]*vbasis[j];
			}
		}
	}

	/*DO NOT Transform Coordinate System: this stiffness matrix is already expressed in tangential coordinates*/

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(vbasis);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixFSViscous(Element* element){/*{{{*/

	/*Intermediaries*/
	int         i,dim;
	IssmDouble  viscosity,FSreconditioning,Jdet;
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke   = element->NewElementMatrix(FSvelocityEnum);
	IssmDouble* vdbasis = xNew<IssmDouble>(dim*vnumnodes);
	IssmDouble* pbasis  = xNew<IssmDouble>(pnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	Input* vx_input=element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){vz_input=element->GetInput(VzEnum); _assert_(vz_input);}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivativesVelocity(vdbasis,xyz_list,gauss);
		element->NodalFunctionsPressure(pbasis,gauss);
		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);

		IssmDouble factor = gauss->weight*Jdet*viscosity;
		IssmDouble factorrecond = gauss->weight*Jdet*FSreconditioning;
		if(dim==2 || dim==3){
			/*Stress balance*/
			for(int i=0;i<vnumnodes;i++){
				for(int j=0;j<vnumnodes;j++){
					for (int p=0;p<dim;p++){
						for (int q=0;q<dim;q++){
							/* diagonal only */
							Ke->values[(dim*i+p)*numdof+dim*j+p] += factor*(vdbasis[q*vnumnodes+j]*vdbasis[q*vnumnodes+i]);
							/* All the entries */
							Ke->values[(dim*i+p)*numdof+dim*j+q] += factor*(vdbasis[p*vnumnodes+j]*vdbasis[q*vnumnodes+i]);
						}
					}
				}
				for(int k=0;k<pnumnodes;k++){
					for (int p=0;p<dim;p++){
						Ke->values[(dim*i+p)*numdof+dim*vnumnodes+k] += -factorrecond*pbasis[k]*vdbasis[p*vnumnodes+i];
					}
				}
			}
			/*Incompressibility*/
			for(int k=0;k<pnumnodes;k++){
				for(int j=0;j<vnumnodes;j++){
					for (int p=0;p<dim;p++){
						Ke->values[(dim*vnumnodes+k)*numdof+dim*j+p] += -factorrecond*vdbasis[p*vnumnodes+j]*pbasis[k];
					}
				}
			}
		}
	}

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(pbasis);
	xDelete<IssmDouble>(vdbasis);
	xDelete<int>(cs_list);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreatePressureMassMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	int         i,dim;
	IssmDouble  FSreconditioning,Jdet;
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke  = element->NewElementMatrix(FSvelocityEnum);
	IssmDouble* pbasis = xNew<IssmDouble>(pnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	//element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsPressure(pbasis,gauss);

		IssmDouble factor = gauss->weight*Jdet;
		if(dim==3 || dim==2){
			/*Pressure mass matrix*/
			for(int k=0;k<pnumnodes;k++){
				for(int j=0;j<pnumnodes;j++){
					Ke->values[(dim*vnumnodes+k)*numdof+dim*vnumnodes+j] += factor*(pbasis[j]*pbasis[k]);
				}
			}
		}else{
			_error_("STOP");
		}
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(pbasis);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateSchurPrecondMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	int         i,dim;
	IssmDouble  viscosity,FSreconditioning,Jdet;
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke  = element->NewElementMatrix(FSvelocityEnum);
	IssmDouble* pbasis = xNew<IssmDouble>(pnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	//element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	Input* vx_input=element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){vz_input=element->GetInput(VzEnum); _assert_(vz_input);}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsPressure(pbasis,gauss);
		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);

		IssmDouble factor = gauss->weight*1./viscosity*Jdet;
		if(dim==3 || dim==2){
			/*Pressure mass matrix*/
			for(int k=0;k<pnumnodes;k++){
				for(int j=0;j<pnumnodes;j++){
					Ke->values[(dim*vnumnodes+k)*numdof+dim*vnumnodes+j] += factor*(pbasis[j]*pbasis[k]);
				}
			}
		}else{
			_error_("STOP");
		}
	}

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(pbasis);
	return Ke;
}/*}}}*/

ElementMatrix* StressbalanceAnalysis::CreateKMatrixFSViscousGLS(Element* element){/*{{{*/

	/*Intermediaries*/
	int         i,dim;
	IssmDouble  viscosity,FSreconditioning,Jdet;
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;

	/* only work for P1 elements */
	if (dim == 2) {_assert_(vnumnodes==3);}
	else if (dim == 3) {_assert_(vnumnodes==6);}
	else {_error_("GLS is not implemented except for 2D and 3D problems.");}

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke   = element->NewElementMatrix(FSvelocityEnum);
	IssmDouble* vdbasis = xNew<IssmDouble>(dim*vnumnodes);
	IssmDouble* pbasis  = xNew<IssmDouble>(pnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	Input* vx_input=element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){vz_input=element->GetInput(VzEnum); _assert_(vz_input);}

	/* prepare viscosity gradient for GLS */
	IssmDouble NodalViscosity[6];
	IssmDouble gradViscos[3];
	IssmDouble etapq, s, Tau, mk, hk;
	IssmDouble hx, hy, hz;
	IssmDouble SU[3*(3+1)*6];
    Gauss* vert = element->NewGauss();
	/* Compute the nodal values of the viscosity */
	for(int i=0;i<vnumnodes;i++){
    	vert->GaussNode(element->element_type, i);
		element->material->ViscosityFS(&NodalViscosity[i],dim,xyz_list,vert,vx_input,vy_input,vz_input);
	}
	delete vert;

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(5);
	while(gauss->next()){
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsDerivativesVelocity(vdbasis,xyz_list,gauss);
		element->NodalFunctionsPressure(pbasis,gauss);
		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);

		/* compute viscosity gradient at the gaussian point */
		element->ValueP1DerivativesOnGauss(&gradViscos[0],NodalViscosity,xyz_list,gauss);

		/*  weight*detJ */
		s = gauss->weight*Jdet;

		/* GLS Stabilization */
		mk = 1.0 /3.0;
		element->ElementSizes(&hx,&hy,&hz);
		// hk = max(max(hx, hy), hz);
		hk = max(hx, hy);
		Tau = - mk * hk * hk * 0.125 / viscosity;

		for (int q=0;q<dim;q++){
			/* column 1-3 for the velocities */
			for (int p=0;p<dim;p++){
				for(int i=0;i<vnumnodes;i++){
					SU[q*numdof+i*dim+p] = (-gradViscos[p])*vdbasis[q*vnumnodes+i];
				}
			}
			/* add the diagnal components */
			for(int i=0;i<vnumnodes;i++){
				for (int p=0;p<dim;p++){
					SU[q*numdof+i*dim+q] += (-gradViscos[p])*vdbasis[p*vnumnodes+i];
				}
			}
			/* column 4 for the pressure */
			for(int i=0;i<pnumnodes;i++){
					SU[q*numdof+dim*vnumnodes+i] = FSreconditioning*vdbasis[q*vnumnodes+i];
			}
		}

		if(dim==2 || dim==3){
			/*Stress balance*/
			for(int i=0;i<vnumnodes;i++){
				for(int j=0;j<vnumnodes;j++){
					for (int p=0;p<dim;p++){
						for (int q=0;q<dim;q++){
							/* diagonal only */
							Ke->values[(dim*i+p)*numdof+dim*j+p] += s*viscosity*(vdbasis[q*vnumnodes+j]*vdbasis[q*vnumnodes+i]);
							/* All the entries */
							Ke->values[(dim*i+p)*numdof+dim*j+q] += s*viscosity*(vdbasis[p*vnumnodes+j]*vdbasis[q*vnumnodes+i]);
						}
					}
				}
				for(int k=0;k<pnumnodes;k++){
					for (int p=0;p<dim;p++){
						Ke->values[(dim*i+p)*numdof+dim*vnumnodes+k] += s*FSreconditioning*(-pbasis[k]*vdbasis[p*vnumnodes+i]);
					}
				}
			}
			/*Incompressibility*/
			for(int k=0;k<pnumnodes;k++){
				for(int j=0;j<vnumnodes;j++){
					for (int p=0;p<dim;p++){
						Ke->values[(dim*vnumnodes+k)*numdof+dim*j+p] += s*(-FSreconditioning*vdbasis[p*vnumnodes+j]*pbasis[k]);
					}
				}
			}

			/* GLS */
			for(int i=0;i<numdof;i++){
				for(int j=0;j<numdof;j++){
					for (int p=0;p<dim;p++){
							Ke->values[i*numdof+j] += Tau * s * SU[p*numdof+i]*SU[p*numdof+j];
					}
				}
			}
		}
	}

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(pbasis);
	xDelete<IssmDouble>(vdbasis);
	xDelete<int>(cs_list);
	return Ke;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorFSViscousGLS(Element* element){/*{{{*/

	int         i,dim;
	IssmDouble  Jdet,forcex,forcey,forcez;
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Initialize vectors*/
	ElementVector* pe     = element->NewElementVector(FSvelocityEnum);
	IssmDouble*    vbasis = xNew<IssmDouble>(vnumnodes);
	IssmDouble*    vdbasis = xNew<IssmDouble>(dim*vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	IssmDouble  rho_ice =element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  gravity =element->FindParam(ConstantsGEnum);
	Input*      loadingforcex_input=element->GetInput(LoadingforceXEnum);  _assert_(loadingforcex_input);
	Input*      loadingforcey_input=element->GetInput(LoadingforceYEnum);  _assert_(loadingforcey_input);
	Input*      loadingforcez_input=NULL;
	if(dim==3){
		loadingforcez_input=element->GetInput(LoadingforceZEnum);  _assert_(loadingforcez_input);
	}

	/* prepare viscosity gradient for GLS */
	Input* vx_input=element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){vz_input=element->GetInput(VzEnum); _assert_(vz_input);}

	IssmDouble viscosity,FSreconditioning;
	IssmDouble NodalViscosity[6];
	IssmDouble gradViscos[3];
	IssmDouble etapq, s, Tau, mk, hk;
	IssmDouble hx, hy, hz;
	IssmDouble SU[3*(3+1)*6];
    Gauss* vert = element->NewGauss();

	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	/* Compute the nodal values of the viscosity */
	for(int i=0;i<vnumnodes;i++){
    	vert->GaussNode(element->element_type, i);
		element->material->ViscosityFS(&NodalViscosity[i],dim,xyz_list,vert,vx_input,vy_input,vz_input);
	}
	delete vert;

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);
		element->NodalFunctionsDerivativesVelocity(vdbasis,xyz_list,gauss);

		loadingforcex_input->GetInputValue(&forcex,gauss);
		loadingforcey_input->GetInputValue(&forcey,gauss);
		if(dim==3) loadingforcez_input->GetInputValue(&forcez,gauss);

		/* compute viscosity gradient at the gaussian point */
		element->ValueP1DerivativesOnGauss(&gradViscos[0],NodalViscosity,xyz_list,gauss);
		/*  weight*detJ */
		s = gauss->weight*Jdet;
		/* GLS Stabilization */
		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);
		mk = 1.0 /3.0;
		element->ElementSizes(&hx,&hy,&hz);
		// hk = max(max(hx, hy), hz);
		hk = max(hx, hy);
		Tau = - mk * hk * hk * 0.125 / viscosity;

		for (int q=0;q<dim;q++){
			/* column 1-3 for the velocities */
			for (int p=0;p<dim;p++){
				for(int i=0;i<vnumnodes;i++){
					SU[q*numdof+i*dim+p] = (-gradViscos[p])*vdbasis[q*vnumnodes+i];
				}
			}
			/* add the diagnal components */
			for(int i=0;i<vnumnodes;i++){
				for (int p=0;p<dim;p++){
					SU[q*numdof+i*dim+q] += (-gradViscos[p])*vdbasis[p*vnumnodes+i];
				}
			}
			/* column 4 for the pressure */
			for(int i=0;i<pnumnodes;i++){
					SU[q*numdof+dim*vnumnodes+i] = FSreconditioning*vdbasis[q*vnumnodes+i];
			}
		}
		IssmDouble factor = rho_ice*Jdet*gauss->weight;
		for(i=0;i<vnumnodes;i++){
			pe->values[i*dim+0] += factor*forcex *vbasis[i];
			pe->values[i*dim+1] += factor*forcey *vbasis[i];
			if(dim==3) pe->values[i*dim+2] += factor*forcez*vbasis[i];

			pe->values[i*dim+dim-1] += -factor*gravity*vbasis[i];
		}

		for(int i=0;i<numdof;i++){
			pe->values[i] += Tau*factor*(-gravity)*SU[(dim-1)*numdof+i];
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<int>(cs_list);
	xDelete<IssmDouble>(vbasis);
	xDelete<IssmDouble>(vdbasis);
	xDelete<IssmDouble>(xyz_list);
	return pe;
}/*}}}*/

ElementMatrix* StressbalanceAnalysis::CreateKMatrixFSViscousLA(Element* element){/*{{{*/

	/*Intermediaries*/
	int         i,dim,epssize;
	IssmDouble  r,rl,Jdet,viscosity,DU,DUl;
	IssmDouble	normal[3];
	IssmDouble *xyz_list = NULL;
	IssmDouble *xyz_list_base = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);
	element->FindParam(&r,AugmentedLagrangianREnum);
	if(dim==2) epssize = 3;
	else       epssize = 6;

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->GetNumberOfNodes(P1Enum);
	int lnumnodes = element->GetNumberOfNodes(P2Enum);
	int numdof    = vnumnodes*dim;
	int pnumdof   = pnumnodes;
	int lnumdof   = lnumnodes;

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke       = element->NewElementMatrix(FSvelocityEnum);
	IssmDouble*    B        = xNew<IssmDouble>(epssize*numdof);
	IssmDouble*    Bprime   = xNew<IssmDouble>(epssize*numdof);
	IssmDouble*    BtBUzawa = xNewZeroInit<IssmDouble>(numdof*pnumdof);
	IssmDouble*    BU       = xNew<IssmDouble>(pnumdof);
	IssmDouble*    BprimeU  = xNew<IssmDouble>(numdof);
	IssmDouble*    D        = xNewZeroInit<IssmDouble>(epssize*epssize);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vx_input = element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input = element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){vz_input = element->GetInput(VzEnum); _assert_(vz_input);}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		this->GetBFSvel(B,element,dim,xyz_list,gauss);
		this->GetBFSprimevel(Bprime,element,dim,xyz_list,gauss);

		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);
		for(i=0;i<epssize;i++)   D[i*epssize+i] = 2*viscosity*gauss->weight*Jdet;

		TripleMultiply(B,epssize,numdof,1,
					D,epssize,epssize,0,
					Bprime,epssize,numdof,0,
					&Ke->values[0],1);

		this->GetBFSUzawa(BU,element,dim,xyz_list,gauss);
		this->GetBFSprimeUzawa(BprimeU,element,dim,xyz_list,gauss);

		DU = gauss->weight*Jdet*sqrt(r);

		TripleMultiply(BU,1,pnumdof,1,
					&DU,1,1,0,
					BprimeU,1,numdof,0,
					BtBUzawa,1);
	}

	/*The pressure augmentation should not be transformed*/
	MatrixMultiply(BtBUzawa,pnumdof,numdof,1,
				BtBUzawa,pnumdof,numdof,0,
				&Ke->values[0],1);

	if(element->IsOnBase() && 0){
		element->FindParam(&rl,AugmentedLagrangianRlambdaEnum);
		element->GetVerticesCoordinatesBase(&xyz_list_base);
		element->NormalBase(&normal[0],xyz_list_base);

		IssmDouble* Dlambda  = xNewZeroInit<IssmDouble>(dim*dim);
		IssmDouble* C        = xNewZeroInit<IssmDouble>(dim*lnumdof);
		IssmDouble* Cprime   = xNewZeroInit<IssmDouble>(dim*numdof);
		IssmDouble* CtCUzawa = xNewZeroInit<IssmDouble>(numdof*lnumdof);

		delete gauss;
		gauss = element->NewGaussBase(5);
		while(gauss->next()){

			element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
			this->GetCFS(C,element,dim,xyz_list,gauss);
			this->GetCFSprime(Cprime,element,dim,xyz_list,gauss);
			IssmDouble factor = gauss->weight*Jdet;
			for(i=0;i<dim;i++) Dlambda[i*dim+i] = factor*sqrt(normal[i]*normal[i])*sqrt(rl);
			TripleMultiply(C,dim,lnumdof,1,
						Dlambda,dim,dim,0,
						Cprime,dim,numdof,0,
						CtCUzawa,1);
		}

		/*The sigma naugmentation should not be transformed*/
		MatrixMultiply(CtCUzawa,lnumdof,numdof,1,
					CtCUzawa,lnumdof,numdof,0,
					&Ke->values[0],1);

		/*Delete base part*/
		xDelete<IssmDouble>(Dlambda);
		xDelete<IssmDouble>(C);
		xDelete<IssmDouble>(Cprime);
		xDelete<IssmDouble>(CtCUzawa);
		xDelete<IssmDouble>(xyz_list_base);
	}

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(D);
	xDelete<IssmDouble>(Bprime);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(BprimeU);
	xDelete<IssmDouble>(BU);
	xDelete<IssmDouble>(BtBUzawa);
	xDelete<int>(cs_list);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixFSViscousXTH(Element* element){/*{{{*/

	/*Intermediaries*/
	int         i,dim,epssize;
	IssmDouble  r,FSreconditioning,Jdet;
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);
	element->FindParam(&r,AugmentedLagrangianREnum);
	if(dim==2) epssize = 3;
	else       epssize = 6;

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;
	int bsize     = epssize + 2;

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix(FSvelocityEnum);
	IssmDouble*    B      = xNew<IssmDouble>(bsize*numdof);
	IssmDouble*    Bprime = xNew<IssmDouble>(bsize*numdof);
	IssmDouble*    D      = xNewZeroInit<IssmDouble>(bsize*bsize);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	Input* vx_input=element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input;
	if(dim==3){vz_input=element->GetInput(VzEnum); _assert_(vz_input);}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss = element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		this->GetBFS(B,element,dim,xyz_list,gauss);
		this->GetBFSprime(Bprime,element,dim,xyz_list,gauss);

		IssmDouble factor = r*gauss->weight*Jdet;
		IssmDouble factor2 = - FSreconditioning*gauss->weight*Jdet;
		for(i=0;i<epssize;i++)     D[i*bsize+i] = factor;
		for(i=epssize;i<bsize;i++) D[i*bsize+i] = factor2;

		TripleMultiply(B,bsize,numdof,1,
					D,bsize,bsize,0,
					Bprime,bsize,numdof,0,
					&Ke->values[0],1);
	}

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(D);
	xDelete<IssmDouble>(Bprime);
	xDelete<IssmDouble>(B);
	xDelete<int>(cs_list);
	return Ke;
}/*}}}*/
#ifdef FSANALYTICAL
ElementMatrix* StressbalanceAnalysis::CreateKMatrixFSFriction(Element* element){/*{{{*/

	if(element->IsAllFloating() || !element->IsOnBase()) return NULL;

	/*If on water or not FS, skip stiffness: */
	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation!=FSApproximationEnum && approximation!=SSAFSApproximationEnum && approximation!=HOFSApproximationEnum) return NULL;

	/*Intermediaries*/
	int         i,dim;
	IssmDouble  alpha2,Jdet;
	IssmDouble  x_coord,y_coord,z_coord;
	IssmDouble *xyz_list_base = NULL;
	IssmDouble *xyz_list      = NULL;
	Gauss*      gauss         = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke = element->NewElementMatrix(FSvelocityEnum);
	IssmDouble*    B  = xNew<IssmDouble>(dim*numdof);
	IssmDouble*    D  = xNewZeroInit<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->GetVerticesCoordinates(&xyz_list);
	Input* vx_input         = element->GetInput(VxEnum);      _assert_(vx_input);
	Input* vy_input         = element->GetInput(VyEnum);      _assert_(vy_input);
	Input* vz_input         = NULL;
	if(dim==3){    vz_input = element->GetInput(VzEnum);      _assert_(vz_input);}

	/* Start  looping on the number of gaussian points: */
	gauss=element->NewGaussBase(10);
	while(gauss->next()){

		x_coord=element->GetXcoord(xyz_list,gauss);
		y_coord=element->GetYcoord(xyz_list,gauss);
		if(dim==3) z_coord=element->GetZcoord(xyz_list,gauss);
		else z_coord=0.;

		alpha2=alpha(x_coord,y_coord,z_coord,FSANALYTICAL);

		this->GetBFSFriction(B,element,dim,xyz_list_base,gauss);
		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		IssmDouble factor = alpha2*gauss->weight*Jdet;
		for(int i=0;i<dim;i++) D[i*dim+i] = factor; //taub_x = -alpha2 v_x (same for y)

		TripleMultiply(B,dim,numdof,1,
					D,dim,dim,0,
					B,dim,numdof,0,
					&Ke->values[0],1);
	}

	/*DO NOT Transform Coordinate System: this stiffness matrix is already expressed in tangential coordinates*/

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(D);
	return Ke;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorFS(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	ElementVector* pe = NULL;

	ElementVector* pe1=CreatePVectorFSViscous(element);
	ElementVector* pe2=CreatePVectorFSFriction(element);
	ElementVector* pe3=CreatePVectorFSStress(element);
	pe =new ElementVector(pe1,pe2,pe3);
	delete pe1;
	delete pe2;
	delete pe3;

	/*clean-up and return*/
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorFSFriction(Element* element){/*{{{*/

	if(!element->IsOnBase()) return NULL;

	/*Intermediaries*/
	int         dim;
	IssmDouble  alpha2,Jdet;
	IssmDouble  bed_normal[3];
	IssmDouble *xyz_list_base = NULL;
	Gauss*      gauss         = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();

	/*Initialize Element matrix and vectors*/
	ElementVector* pe = element->NewElementVector(FSvelocityEnum);
	IssmDouble*    vbasis = xNew<IssmDouble>(vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	Input*  alpha2_input=element->GetInput(FrictionCoefficientEnum); _assert_(alpha2_input);

	/* Start  looping on the number of gaussian points: */
	gauss=element->NewGaussBase(3);
	while(gauss->next()){

		alpha2_input->GetInputValue(&alpha2, gauss);
		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);
		element->NormalBase(&bed_normal[0],xyz_list_base);

		IssmDouble factor = alpha2*gauss->weight*Jdet;
		for(int i=0;i<vnumnodes;i++){
			pe->values[i*dim+0] += - factor*vbasis[i]*bed_normal[1];
			pe->values[i*dim+1] += factor*vbasis[i]*bed_normal[0];
			if(dim==3){
				pe->values[i*dim+2]+= factor*vbasis[i];
			}
		}

	}

	/*DO NOT Transform Coordinate System: this stiffness matrix is already expressed in tangential coordinates*/

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(vbasis);
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorFSStress(Element* element){/*{{{*/

	/*Skipping for now*/
	return NULL;
	if(!element->IsOnBase()) return NULL;

	/*Intermediaries*/
	int         dim;
	IssmDouble  sigmann,sigmant,Jdet,bedslope,beta;
	IssmDouble *xyz_list_base = NULL;
	Gauss*      gauss         = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();

	/*Initialize Element matrix and vectors*/
	ElementVector* pe = element->NewElementVector(FSvelocityEnum);
	IssmDouble*    vbasis = xNew<IssmDouble>(vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	Input*  sigmann_input=element->GetInput(VzEnum); _assert_(sigmann_input);
	Input*  sigmant_input=element->GetInput(TemperatureEnum); _assert_(sigmant_input);
	Input*  bedslope_input=element->GetInput(BedSlopeXEnum);     _assert_(bedslope_input);

	/* Start  looping on the number of gaussian points: */
	gauss=element->NewGaussBase(3);
	while(gauss->next()){

		sigmann_input->GetInputValue(&sigmann, gauss);
		sigmant_input->GetInputValue(&sigmant, gauss);
		bedslope_input->GetInputValue(&bedslope, gauss);
		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);

		beta=sqrt(1+bedslope*bedslope);
		IssmDouble factor = - (1./beta)*gauss->weight*Jdet;
		for(int i=0;i<vnumnodes;i++){
			pe->values[i*dim+0] += factor*(-bedslope*sigmann + sigmant)*vbasis[i];
			pe->values[i*dim+1] += factor*(sigmann + bedslope*sigmant)*vbasis[i];
			if(dim==3){
				//pe->values[i*dim+2]+= alpha2*gauss->weight*Jdet*vbasis[i];
				_error_("3d not supported yet");
			}
		}

	}

	/*DO NOT Transform Coordinate System: this stiffness matrix is already expressed in tangential coordinates*/

	/*Clean up and return*/
	delete gauss;
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(vbasis);
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorFSViscous(Element* element){/*{{{*/

	int         i,dim,fe_FS;
	IssmDouble  x_coord,y_coord,z_coord;
	IssmDouble  Jdet,forcex,forcey,forcez;
	IssmDouble *xyz_list = NULL;

	element->FindParam(&fe_FS,FlowequationFeFSEnum);
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Initialize vectors*/
	ElementVector* pe     = element->NewElementVector(FSvelocityEnum);
	IssmDouble*    vbasis = xNew<IssmDouble>(vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);

		x_coord=element->GetXcoord(xyz_list,gauss);
		y_coord=element->GetYcoord(xyz_list,gauss);
		if(dim==3) z_coord=element->GetZcoord(xyz_list,gauss);
		else z_coord=0.;

		forcex=fx(x_coord,y_coord,z_coord,FSANALYTICAL);
		forcey=fy(x_coord,y_coord,z_coord,FSANALYTICAL);
		forcez=fz(x_coord,y_coord,z_coord,FSANALYTICAL);

		IssmDouble factor = Jdet*gauss->weight;
		for(i=0;i<vnumnodes;i++){
			pe->values[i*dim+0] += forcex *factor*vbasis[i];
			pe->values[i*dim+1] += forcey *factor*vbasis[i];
			if(dim==3) pe->values[i*dim+2] += forcez *factor*vbasis[i];
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<int>(cs_list);
	xDelete<IssmDouble>(vbasis);
	xDelete<IssmDouble>(xyz_list);
	if(fe_FS==XTaylorHoodEnum){
		ElementVector* pe2=CreatePVectorFSViscousXTH(element);
		ElementVector* pe3 = new ElementVector(pe,pe2);
		delete pe;
		delete pe2;
		return pe3;
	}
	else if(fe_FS==LATaylorHoodEnum || fe_FS==LACrouzeixRaviartEnum){
		ElementVector* pe2=CreatePVectorFSViscousLA(element);
		ElementVector* pe3 = new ElementVector(pe,pe2);
		delete pe;
		delete pe2;
		return pe3;
	}
	return pe;
}/*}}}*/
#else
ElementMatrix* StressbalanceAnalysis::CreateKMatrixFSFriction(Element* element){/*{{{*/

	if(element->IsAllFloating() || !element->IsOnBase()) return NULL;

	/*If on water or not FS, skip stiffness: */
	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation!=FSApproximationEnum && approximation!=SSAFSApproximationEnum && approximation!=HOFSApproximationEnum) return NULL;

	/*Intermediaries*/
	bool        mainlyfloating;
	int         dim,domaintype;
	int         friction_style,point1;
	IssmDouble  alpha2,Jdet,fraction1,fraction2;
	IssmDouble  gllevelset,phi=1.;
	IssmDouble *xyz_list_base = NULL;
	Gauss*      gauss         = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);
	element->FindParam(&domaintype,DomainTypeEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke = element->NewElementMatrix(FSvelocityEnum);
	IssmDouble* vbasis=xNew<IssmDouble>(vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->FindParam(&friction_style,GroundinglineFrictionInterpolationEnum);
	Input* gllevelset_input = NULL;

	/*build friction object, used later on: */
	Friction* friction=new Friction(element,dim==3?3:1);

	/*Recover portion of element that is grounded*/
	if(!(friction_style==SubelementFriction2Enum)) phi=element->GetGroundedPortion(xyz_list_base);
	if(friction_style==SubelementFriction2Enum){
		if(domaintype==Domain2DverticalEnum) _error_("Subelement Friction 2 not implemented yet for Flowline");
		gllevelset_input=element->GetInput(MaskOceanLevelsetEnum); _assert_(gllevelset_input);
		element->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating, MaskOceanLevelsetEnum,0);
		//gauss = element->NewGauss(point1,fraction1,fraction2,mainlyfloating,2);
		gauss=element->NewGaussBase(3);
	}
	else{
		gauss=element->NewGaussBase(3);
	}

	/* Start  looping on the number of gaussian points: */
	while(gauss->next()){

		friction->GetAlpha2(&alpha2,gauss);
		if(friction_style==SubelementFriction1Enum) alpha2=phi*alpha2;
		else if(friction_style==SubelementFriction2Enum){
			gllevelset_input->GetInputValue(&gllevelset, gauss);
			if(gllevelset<0.) alpha2=0.;
		}
		else if(friction_style==NoFrictionOnPartiallyFloatingEnum){
			if (phi<0.99999999) alpha2=0.;
		}
		else  _error_("friction interpolation "<<EnumToStringx(friction_style)<<" not implemented yet");

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);

		IssmDouble factor = alpha2*gauss->weight*Jdet;
		if(dim==3){
			/*Stress balance*/
			for(int i=0;i<vnumnodes;i++){
				for(int j=0;j<vnumnodes;j++){
					Ke->values[(3*i+0)*numdof+3*j+0] += factor*vbasis[i]*vbasis[j];
					Ke->values[(3*i+1)*numdof+3*j+1] += factor*vbasis[i]*vbasis[j];
					Ke->values[(3*i+2)*numdof+3*j+2] += factor*vbasis[i]*vbasis[j];
				}
			}
		}
		else{
			/*Stress balance*/
			for(int i=0;i<vnumnodes;i++){
				for(int j=0;j<vnumnodes;j++){
					Ke->values[(2*i+0)*numdof+2*j+0] += factor*vbasis[i]*vbasis[j];
					Ke->values[(2*i+1)*numdof+2*j+1] += factor*vbasis[i]*vbasis[j];
				}
			}
		}
	}

	/*DO NOT Transform Coordinate System: this stiffness matrix is already expressed in tangential coordinates*/

	/*Clean up and return*/
	delete gauss;
	delete friction;
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(vbasis);
	return Ke;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorFS(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	ElementVector* pe = NULL;
	int fe_FS;
	element->FindParam(&fe_FS,FlowequationFeFSEnum);

	if(fe_FS==XTaylorHoodEnum){
		ElementVector* pe1=CreatePVectorFSViscous(element);
		ElementVector* pe2=CreatePVectorFSShelf(element);
		ElementVector* pe3=CreatePVectorFSFront(element);
		ElementVector* petemp =new ElementVector(pe1,pe2,pe3);
		ElementVector* pe4=CreatePVectorFSViscousXTH(element);
		pe = new ElementVector(petemp,pe4);
		delete pe1;
		delete pe2;
		delete pe3;
		delete petemp;
		delete pe4;
	}
	else if(fe_FS==LATaylorHoodEnum || fe_FS==LACrouzeixRaviartEnum){
		ElementVector* pe1=CreatePVectorFSViscous(element);
		ElementVector* pe2=CreatePVectorFSShelf(element);
		ElementVector* pe3=CreatePVectorFSFront(element);
		ElementVector* petemp =new ElementVector(pe1,pe2,pe3);
		ElementVector* pe4=CreatePVectorFSViscousLA(element);
		pe = new ElementVector(petemp,pe4);
		delete pe1;
		delete pe2;
		delete pe3;
		delete petemp;
		delete pe4;
	}
	else if(fe_FS==P1P1GLSEnum) {
		ElementVector* pe1=CreatePVectorFSViscousGLS(element);
		ElementVector* pe2=CreatePVectorFSShelf(element);
		ElementVector* pe3=CreatePVectorFSFront(element);
		pe =new ElementVector(pe1,pe2,pe3);
		delete pe1;
		delete pe2;
		delete pe3;
	}
	else{
		ElementVector* pe1=CreatePVectorFSViscous(element);
		ElementVector* pe2=CreatePVectorFSShelf(element);
		ElementVector* pe3=CreatePVectorFSFront(element);
		pe =new ElementVector(pe1,pe2,pe3);
		delete pe1;
		delete pe2;
		delete pe3;
	}

	/*clean-up and return*/
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorFSViscous(Element* element){/*{{{*/

	int         i,dim;
	IssmDouble  Jdet,forcex,forcey,forcez;
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Initialize vectors*/
	ElementVector* pe     = element->NewElementVector(FSvelocityEnum);
	IssmDouble*    vbasis = xNew<IssmDouble>(vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	IssmDouble  rho_ice =element->FindParam(MaterialsRhoIceEnum);
	IssmDouble  gravity =element->FindParam(ConstantsGEnum);
	Input*      loadingforcex_input=element->GetInput(LoadingforceXEnum);  _assert_(loadingforcex_input);
	Input*      loadingforcey_input=element->GetInput(LoadingforceYEnum);  _assert_(loadingforcey_input);
	Input*      loadingforcez_input=NULL;
	if(dim==3){
		loadingforcez_input=element->GetInput(LoadingforceZEnum);  _assert_(loadingforcez_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);

		loadingforcex_input->GetInputValue(&forcex,gauss);
		loadingforcey_input->GetInputValue(&forcey,gauss);
		if(dim==3) loadingforcez_input->GetInputValue(&forcez,gauss);

		IssmDouble factor = rho_ice*Jdet*gauss->weight;
		for(i=0;i<vnumnodes;i++){
			pe->values[i*dim+0] += factor*forcex *vbasis[i];
			pe->values[i*dim+1] += factor*forcey *vbasis[i];
			if(dim==3) pe->values[i*dim+2] += factor*forcez*vbasis[i];

			pe->values[i*dim+dim-1] += -factor*gravity*vbasis[i];
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<int>(cs_list);
	xDelete<IssmDouble>(vbasis);
	xDelete<IssmDouble>(xyz_list);
	return pe;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixFSFrictionNitsche(Element* element){/*{{{*/

	if(element->IsAllFloating() || !element->IsOnBase()) return NULL;

	/*If on water or not FS, skip stiffness: */
	int approximation;
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation!=FSApproximationEnum && approximation!=SSAFSApproximationEnum && approximation!=HOFSApproximationEnum) return NULL;

	/*Intermediaries*/
	bool        mainlyfloating;
	int         dim,domaintype;
	int         friction_style,point1;
	IssmDouble  alpha2,Jdet,fraction1,fraction2;
	IssmDouble  viscosity, FSreconditioning;
	IssmDouble  gllevelset,phi=1.;
	IssmDouble *xyz_list_base = NULL;
	IssmDouble *xyz_list = NULL;
	Gauss*      gauss         = NULL;
	/*coefficient of Nitsche's method*/
	IssmDouble gamma, hx, hy, hz;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numdof    = vnumnodes*dim + pnumnodes;

	/*Initialize Element matrix and vectors*/
	ElementMatrix* Ke = element->NewElementMatrix(FSvelocityEnum);
	IssmDouble* vbasis = xNew<IssmDouble>(vnumnodes);
	IssmDouble* vdbasis = xNew<IssmDouble>(dim*vnumnodes);
	IssmDouble* pbasis  = xNew<IssmDouble>(pnumnodes);

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(int i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(int i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(int i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&friction_style,GroundinglineFrictionInterpolationEnum);
	Input* gllevelset_input = NULL;
	Input* vx_input=element->GetInput(VxEnum);     _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum);     _assert_(vy_input);
	Input* vz_input = NULL;
	if(dim==3){vz_input=element->GetInput(VzEnum); _assert_(vz_input);}

	/*build friction object, used later on: */
	Friction* friction=new Friction(element,dim==3?3:1);

	/*Recover portion of element that is grounded*/
	if(!(friction_style==SubelementFriction2Enum)) phi=element->GetGroundedPortion(xyz_list_base);
	if(friction_style==SubelementFriction2Enum){
		if(domaintype==Domain2DverticalEnum) _error_("Subelement Friction 2 not implemented yet for Flowline");
		gllevelset_input=element->GetInput(MaskOceanLevelsetEnum); _assert_(gllevelset_input);
		element->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating, MaskOceanLevelsetEnum,0);
		gauss=element->NewGaussBase(3);
	}
	else{
		gauss=element->NewGaussBase(3);
	}

	element->ElementSizes(&hx,&hy,&hz);
	element->FindParam(&gamma,FeFSNitscheGammaEnum);

	gamma = gamma / sqrt(hx*hx+hy*hy+hz*hz);

	/* Start  looping on the number of gaussian points: */
	while(gauss->next()){
		friction->GetAlpha2(&alpha2,gauss);
		if(friction_style==SubelementFriction1Enum) alpha2=phi*alpha2;
		else if(friction_style==SubelementFriction2Enum){
			gllevelset_input->GetInputValue(&gllevelset, gauss);
			if(gllevelset<0.) alpha2=0.;
		}
		else if(friction_style==NoFrictionOnPartiallyFloatingEnum){
			if (phi<0.99999999) alpha2=0.;
		}
		else  _error_("friction interpolation "<<EnumToStringx(friction_style)<<" not implemented yet");

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);
		/* The full element is needed for calculating derivative of the basis functions */
		element->NodalFunctionsDerivativesVelocity(vdbasis,xyz_list,gauss);
		element->NodalFunctionsPressure(pbasis,gauss);
		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);

		IssmDouble factor =  alpha2*gauss->weight*Jdet;
		for(int i=0;i<vnumnodes;i++){
			for(int j=0;j<vnumnodes;j++){
				for (int d=0; d<dim; d++) {
						Ke->values[(dim*i+d)*numdof+dim*j+d] += factor*vbasis[i]*vbasis[j];
				}
			}
		}
		/* -------- Nitsche terms -------- */
		/* Boundary terms for integration by parts.
			The coefficient matrix of n*sigma*n-gamma*n*u is stored in the following order:
			rows--dimensions, columns--number of nodes.
			If we consider nsigma as a 1d vector, it has exactly the same order as the unknown vector.
		*/
		factor = gauss->weight*Jdet*viscosity;
		IssmDouble factor2 = gauss->weight*Jdet*FSreconditioning;
		for(int i=0;i<vnumnodes;i++){
			for(int j=0;j<vnumnodes;j++){
				/* gamma*(n*u)*(n*v) */
				Ke->values[(dim*i+(dim-1))*numdof+dim*j+(dim-1)] += factor * gamma * vbasis[j] * vbasis[i];
				/* -sigma(v)*(n*u) */
				Ke->values[(dim*i+(dim-1))*numdof+dim*j+(dim-1)] += - factor * 2.0 * (-vbasis[j]) * vdbasis[(dim-1)*vnumnodes+i];
				/* -sigma(u)*(n*v) */
				Ke->values[(dim*i+(dim-1))*numdof+dim*j+(dim-1)] += - factor * 2.0 * vdbasis[(dim-1)*vnumnodes+j] * (-vbasis[i]);
			}
		}
		/* pressure x velocity  component A12, +p*(n*v) */
		for(int k=0;k<pnumnodes;k++){
			for(int i=0;i<vnumnodes;i++){
				Ke->values[(dim*i+dim-1)*numdof+dim*vnumnodes+k] += factor*pbasis[k] * (-vbasis[i]);
			}
		}
		/* velocity x pressure component A21, +(n*u)*q */
		for(int k=0;k<pnumnodes;k++){
			for(int j=0;j<vnumnodes;j++){
				Ke->values[(dim*vnumnodes+k)*numdof+dim*j+dim-1] += factor * (-vbasis[j]) * pbasis[k];
			}
		}
	}

	/*Clean up and return*/
	delete gauss;
	delete friction;
	xDelete<IssmDouble>(xyz_list_base);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(vbasis);
	xDelete<IssmDouble>(vdbasis);
	xDelete<IssmDouble>(pbasis);
	return Ke;
}/*}}}*/
#endif
ElementVector* StressbalanceAnalysis::CreatePVectorFSFront(Element* element){/*{{{*/

	/*If no front, return NULL*/
	if(!element->IsIcefront()) return NULL;

	/*Intermediaries*/
	int         i,dim;
	IssmDouble  Jdet,pressure,surface,sealevel,z;
	IssmDouble	normal[3];
	IssmDouble *xyz_list       = NULL;
	IssmDouble *xyz_list_front = NULL;
	Gauss      *gauss          = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes   = element->NumberofNodesVelocity();
	int pnumnodes   = element->NumberofNodesPressure();
	int numvertices = element->GetNumberOfVertices();

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Initialize vectors*/
	ElementVector* pe     = element->NewElementVector(FSvelocityEnum);
	IssmDouble*    vbasis = xNew<IssmDouble>(vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetIcefrontCoordinates(&xyz_list_front,xyz_list,MaskIceLevelsetEnum);
	element->NormalSection(&normal[0],xyz_list_front);
	Input* surface_input  = element->GetInput(SurfaceEnum); _assert_(surface_input);
	Input* sealevel_input       = element->GetInput(SealevelEnum);       _assert_(sealevel_input);
	IssmDouble  rho_water = element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble  gravity   = element->FindParam(ConstantsGEnum);

	/*Initialize gauss points*/
	IssmDouble zmax=xyz_list[0*3+(dim-1)]; for(int i=1;i<numvertices;i++) if(xyz_list[i*3+(dim-1)]>zmax) zmax=xyz_list[i*3+(dim-1)];
	IssmDouble zmin=xyz_list[0*3+(dim-1)]; for(int i=1;i<numvertices;i++) if(xyz_list[i*3+(dim-1)]<zmin) zmin=xyz_list[i*3+(dim-1)];
	if(zmax>0. && zmin<0.) gauss=element->NewGauss(xyz_list,xyz_list_front,3,30);//refined in vertical because of the sea level discontinuity
	else                   gauss=element->NewGauss(xyz_list,xyz_list_front,3,3);

	/* Start  looping on the number of gaussian points: */
	while(gauss->next()){

		element->JacobianDeterminantSurface(&Jdet,xyz_list_front,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);
		surface_input->GetInputValue(&surface,gauss);
		sealevel_input->GetInputValue(&sealevel,gauss);
		if(dim==3) z=element->GetZcoord(xyz_list,gauss);
		else       z=element->GetYcoord(xyz_list,gauss);
		pressure = rho_water*gravity*min(0.,z-sealevel);//0 if the gaussian point is above water level

		IssmDouble factor = pressure*Jdet*gauss->weight;
		for (int i=0;i<vnumnodes;i++){
			pe->values[dim*i+0]+= factor*normal[0]*vbasis[i];
			pe->values[dim*i+1]+= factor*normal[1]*vbasis[i];
			if(dim==3) pe->values[dim*i+2]+= factor*normal[2]*vbasis[i];
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<int>(cs_list);
	xDelete<IssmDouble>(vbasis);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_front);
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorFSShelf(Element* element){/*{{{*/

	int         i,dim;
	IssmDouble  Jdet,water_pressure,base;
	IssmDouble *xyz_list_base = NULL;

	/*Get basal element*/
	if(!element->IsOnBase() || !element->IsAllFloating()) return NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Initialize vectors*/
	ElementVector* pe     = element->NewElementVector(FSvelocityEnum);
	IssmDouble*    vbasis = xNew<IssmDouble>(vnumnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinatesBase(&xyz_list_base);
	Input*      base_input=element->GetInput(BaseEnum); _assert_(base_input);
	IssmDouble  rho_water=element->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble  gravity  =element->FindParam(ConstantsGEnum);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(5);
	while(gauss->next()){

		element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
		element->NodalFunctionsVelocity(vbasis,gauss);

		base_input->GetInputValue(&base, gauss);
		water_pressure=gravity*rho_water*base;

		IssmDouble factor = -water_pressure*gauss->weight*Jdet;
		for(i=0;i<vnumnodes;i++){
				pe->values[i*dim+(dim-1)]+=factor*vbasis[i];
		}
	}

	/* shelf dampening*/
	int shelf_dampening;
	element->FindParam(&shelf_dampening,StressbalanceShelfDampeningEnum);
	if(shelf_dampening) {
		Input*      mb_input=element->GetInput(BasalforcingsFloatingiceMeltingRateEnum); _assert_(mb_input);
		IssmDouble  dt,mb;
		element->FindParam(&dt,TimesteppingTimeStepEnum);
		gauss->Reset();
		while(gauss->next()){
			element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
			element->NodalFunctionsVelocity(vbasis,gauss);
			mb_input->GetInputValue(&mb, gauss);
			IssmDouble factor = -dt*rho_water*gravity*mb*gauss->weight*Jdet;
			for(i=0;i<vnumnodes;i++){
				pe->values[i*dim+(dim-1)] += factor*vbasis[i];
			}
		}
	}

	/*DO NOT Transform Coordinate System: this stiffness matrix is already expressed in tangential coordinates*/

	/*Clean up and return*/
	delete gauss;
	xDelete<int>(cs_list);
	xDelete<IssmDouble>(vbasis);
	xDelete<IssmDouble>(xyz_list_base);
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorFSViscousLA(Element* element){/*{{{*/

	int         i,dim;
	IssmDouble  Jdet,pressure;
	IssmDouble  bed_normal[3];
	IssmDouble *xyz_list      = NULL;
	IssmDouble *xyz_list_base = NULL;
	Gauss*      gauss         = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(numnodes);
	if(dim==2) for(i=0;i<numnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<numnodes;i++) cs_list[i] = XYZEnum;

	/*Initialize vectors*/
	ElementVector* pe      = element->NewElementVector(FSvelocityEnum);
	IssmDouble*    dbasis  = xNew<IssmDouble>(3*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/*Get pressure and sigmann*/
	Input* pressure_input=element->GetInput(PressureEnum); _assert_(pressure_input);
	Input* sigmann_input =element->GetInput(SigmaNNEnum);  _assert_(sigmann_input);

	gauss=element->NewGauss(5);
	while(gauss->next()){
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);

		pressure_input->GetInputValue(&pressure, gauss);
		element->NodalFunctionsDerivativesVelocity(dbasis,xyz_list,gauss);

		IssmDouble factor = pressure*gauss->weight*Jdet;
		for(i=0;i<numnodes;i++){
			pe->values[i*dim+0] += factor*dbasis[0*numnodes+i];
			pe->values[i*dim+1] += factor*dbasis[1*numnodes+i];
			if(dim==3) pe->values[i*dim+2]+= factor*dbasis[2*numnodes+i];
		}
	}

	if(element->IsOnBase() && 0){
		IssmDouble   sigmann;
		IssmDouble*  vbasis = xNew<IssmDouble>(numnodes);

		element->GetVerticesCoordinatesBase(&xyz_list_base);
		element->NormalBase(&bed_normal[0],xyz_list_base);

		delete gauss;
		gauss=element->NewGaussBase(5);
		while(gauss->next()){

			element->JacobianDeterminantBase(&Jdet,xyz_list_base,gauss);
			element->NodalFunctionsVelocity(vbasis,gauss);
			sigmann_input->GetInputValue(&sigmann, gauss);

			IssmDouble factor = sigmann*gauss->weight*Jdet;
			for(i=0;i<numnodes;i++){
				pe->values[i*dim+0] += factor*bed_normal[0]*vbasis[i];
				pe->values[i*dim+1] += factor*bed_normal[1]*vbasis[i];
				if(dim==3) pe->values[i*dim+2] += factor*bed_normal[2]*vbasis[i];
			}
		}
		xDelete<IssmDouble>(xyz_list_base);
		xDelete<IssmDouble>(vbasis);
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<int>(cs_list);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorFSViscousXTH(Element* element){/*{{{*/

	int         i,tausize,dim;
	IssmDouble  Jdet,r;
	IssmDouble  epsxx,epsyy,epszz,epsxy,epsxz,epsyz;
	IssmDouble  sigmapxx,sigmapyy,sigmapzz,sigmapxy,sigmapxz,sigmapyz;
	IssmDouble *xyz_list = NULL;
	Gauss*      gauss    = NULL;

	/*Get problem dimension*/
	element->FindParam(&dim,DomainDimensionEnum);
	if(dim==2) tausize = 3;
	else       tausize = 6;

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int tnumnodes = element->GetNumberOfVertices();      //Tensors, P1 DG

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i]  = PressureEnum;

	/*Initialize vectors*/
	ElementVector* pe      = element->NewElementVector(FSvelocityEnum);
	IssmDouble*    Dstar   = xNewZeroInit<IssmDouble>((dim*vnumnodes)*(tausize*tnumnodes));
	IssmDouble*    tau     = xNew<IssmDouble>(tausize*tnumnodes);
	IssmDouble*    d       = xNew<IssmDouble>(tausize*tnumnodes);
	IssmDouble*    vdbasis = xNew<IssmDouble>(dim*vnumnodes);
	IssmDouble*    tbasis  = xNew<IssmDouble>(tnumnodes);
	IssmDouble*    D       = xNewZeroInit<IssmDouble>(tausize*tnumnodes*tausize*tnumnodes);

	/*Retrieve all inputs and parameters*/
	element->FindParam(&r,AugmentedLagrangianREnum);
	element->GetVerticesCoordinates(&xyz_list);

	/*Get d and tau*/
	Input* epsxx_input=element->GetInput(StrainRatexxEnum); _assert_(epsxx_input);
	Input* epsyy_input=element->GetInput(StrainRateyyEnum); _assert_(epsyy_input);
	Input* epsxy_input=element->GetInput(StrainRatexyEnum); _assert_(epsxy_input);
	Input* epszz_input=NULL; Input* epsxz_input=NULL; Input* epsyz_input=NULL;
	Input* sigmapxx_input=element->GetInput(DeviatoricStressxxEnum); _assert_(sigmapxx_input);
	Input* sigmapyy_input=element->GetInput(DeviatoricStressyyEnum); _assert_(sigmapyy_input);
	Input* sigmapxy_input=element->GetInput(DeviatoricStressxyEnum); _assert_(sigmapxy_input);
	Input* sigmapzz_input=NULL; Input* sigmapxz_input=NULL; Input* sigmapyz_input=NULL;
	if(dim==3){
		epszz_input=element->GetInput(StrainRatezzEnum); _assert_(epszz_input);
		epsxz_input=element->GetInput(StrainRatexzEnum); _assert_(epsxz_input);
		epsyz_input=element->GetInput(StrainRateyzEnum); _assert_(epsyz_input);
		sigmapzz_input=element->GetInput(DeviatoricStresszzEnum); _assert_(sigmapzz_input);
		sigmapxz_input=element->GetInput(DeviatoricStressxzEnum); _assert_(sigmapxz_input);
		sigmapyz_input=element->GetInput(DeviatoricStressyzEnum); _assert_(sigmapyz_input);
	}

	gauss = element->NewGauss();
	for(int i=0;i<tnumnodes;i++){
		gauss->GaussNode(P1DGEnum,i);

		epsxx_input->GetInputValue(&epsxx,gauss); sigmapxx_input->GetInputValue(&sigmapxx,gauss);
		epsyy_input->GetInputValue(&epsyy,gauss); sigmapyy_input->GetInputValue(&sigmapyy,gauss);
		epsxy_input->GetInputValue(&epsxy,gauss); sigmapxy_input->GetInputValue(&sigmapxy,gauss);
		if(dim==2){
			d[i*tausize+0]=epsxx;  tau[i*tausize+0]=sigmapxx;
			d[i*tausize+1]=epsyy;  tau[i*tausize+1]=sigmapyy;
			d[i*tausize+2]=epsxy;  tau[i*tausize+2]=sigmapxy;
		}
		else{
			epszz_input->GetInputValue(&epszz,gauss); sigmapzz_input->GetInputValue(&sigmapzz,gauss);
			epsxz_input->GetInputValue(&epsxz,gauss); sigmapxz_input->GetInputValue(&sigmapxz,gauss);
			epsyz_input->GetInputValue(&epsyz,gauss); sigmapyz_input->GetInputValue(&sigmapyz,gauss);
			d[i*tausize+0]=epsxx;  tau[i*tausize+0]=sigmapxx;
			d[i*tausize+1]=epsyy;  tau[i*tausize+1]=sigmapyy;
			d[i*tausize+2]=epszz;  tau[i*tausize+2]=sigmapzz;
			d[i*tausize+3]=epsxy;  tau[i*tausize+3]=sigmapxy;
			d[i*tausize+4]=epsxz;  tau[i*tausize+4]=sigmapxz;
			d[i*tausize+5]=epsyz;  tau[i*tausize+5]=sigmapyz;
		}
	}

	/* Start  looping on the number of gaussian points: */
	delete gauss;
	gauss=element->NewGauss(5);
	while(gauss->next()){
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Create Dstar*/
		/*In dim = 2
		 *
		 *       <----------------- tausize ---------------> x tnumnodes
		 *       |  gamma_ij^x         0       gamma_ij^y  | ^
		 * Dij = |                                         | dim
		 *       |     0          gamma_ij^y   gamma_ij^x  | v
		 *                                                   x
		 *                                                   vnumnodes
		 *
		 *In dim = 3
		 *
		 *       |  gamma_ij^x         0          0         gamma_ij^y  gamma_ij^z      0      |
		 *       |                                                                             |
		 * Dij = |     0          gamma_ij^y      0         gamma_ij^x     0        gamma_ij^z |
		 *       |                                                                             |
		 *       |     0               0      gamma_ij^z        0       gamma_ij^x  gamma_ij^y |
		 *
		 * gamma_ij^x = zeta_j dphi_i/dx
		 *
		 * where:
		 *   - zeta_j is the nodal function for the j^th node of the tensor (P1DG)
		 *   - phi_i  is the nodal function for the i^th node of the velocity (P2)*/
		element->NodalFunctionsDerivativesVelocity(vdbasis,xyz_list,gauss);
		element->NodalFunctionsTensor(tbasis,gauss);
		IssmDouble factor = gauss->weight*Jdet;
		if(dim==2){
			for(int i=0;i<vnumnodes;i++){
				for(int j=0;j<tnumnodes;j++){
					Dstar[(i*dim+0)*tausize*tnumnodes + j*tausize+0] += factor*tbasis[j]*vdbasis[0*vnumnodes+i];
					Dstar[(i*dim+0)*tausize*tnumnodes + j*tausize+2] += factor*tbasis[j]*vdbasis[1*vnumnodes+i];

					Dstar[(i*dim+1)*tausize*tnumnodes + j*tausize+1] += factor*tbasis[j]*vdbasis[1*vnumnodes+i];
					Dstar[(i*dim+1)*tausize*tnumnodes + j*tausize+2] += factor*tbasis[j]*vdbasis[0*vnumnodes+i];
				}
			}
		}
		else{
			for(int i=0;i<vnumnodes;i++){
				for(int j=0;j<tnumnodes;j++){
					Dstar[(i*dim+0)*tausize*tnumnodes + j*tausize+0] += factor*tbasis[j]*vdbasis[0*vnumnodes+i];
					Dstar[(i*dim+0)*tausize*tnumnodes + j*tausize+3] += factor*tbasis[j]*vdbasis[1*vnumnodes+i];
					Dstar[(i*dim+0)*tausize*tnumnodes + j*tausize+4] += factor*tbasis[j]*vdbasis[2*vnumnodes+i];

					Dstar[(i*dim+1)*tausize*tnumnodes + j*tausize+1] += factor*tbasis[j]*vdbasis[1*vnumnodes+i];
					Dstar[(i*dim+1)*tausize*tnumnodes + j*tausize+3] += factor*tbasis[j]*vdbasis[0*vnumnodes+i];
					Dstar[(i*dim+1)*tausize*tnumnodes + j*tausize+5] += factor*tbasis[j]*vdbasis[2*vnumnodes+i];

					Dstar[(i*dim+2)*tausize*tnumnodes + j*tausize+2] += factor*tbasis[j]*vdbasis[2*vnumnodes+i];
					Dstar[(i*dim+2)*tausize*tnumnodes + j*tausize+4] += factor*tbasis[j]*vdbasis[0*vnumnodes+i];
					Dstar[(i*dim+2)*tausize*tnumnodes + j*tausize+5] += factor*tbasis[j]*vdbasis[1*vnumnodes+i];
				}
			}
		}
	}

	/*contribution -Dstar tau*/
	for(i=0;i<tausize*tnumnodes;i++) D[i*(tausize*tnumnodes)+i] = -1.;
	TripleMultiply(Dstar,dim*vnumnodes,tausize*tnumnodes,0,
				D,tausize*tnumnodes,tausize*tnumnodes,0,
				tau,tausize*tnumnodes,1,0,
				&pe->values[0],1);

	/*contribution + r Dstar d*/
	for(i=0;i<tausize*tnumnodes;i++) D[i*(tausize*tnumnodes)+i] = +r;
	TripleMultiply(Dstar,dim*vnumnodes,tausize*tnumnodes,0,
				D,tausize*tnumnodes,tausize*tnumnodes,0,
				d,tausize*tnumnodes,1,0,
				&pe->values[0],1);

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,cs_list);

	/*Clean up and return*/
	delete gauss;
	xDelete<int>(cs_list);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(Dstar);
	xDelete<IssmDouble>(d);
	xDelete<IssmDouble>(D);
	xDelete<IssmDouble>(tau);
	xDelete<IssmDouble>(vdbasis);
	xDelete<IssmDouble>(tbasis);
	return pe;
}/*}}}*/
void           StressbalanceAnalysis::GetBFS(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[Bv1 Bv2 ... Bp1 Bp2 ...] where Bvi is of size 3*3.
	 * For node i, Bvi can be expressed in the actual coordinate system
	 * by: 	   Bvi=[ dphi/dx          0        ]
	 *					 [   0           dphi/dy     ]
	 *					 [ 1/2*dphi/dy    1/2*dphi/dx]
	 *					 [   0             0         ]
	 *					 [ dphi/dx         dphi/dy   ]
	 *
	 *         Bpi=[  0    ]
	 *					[  0    ]
	 *					[  0    ]
	 *					[ phi_p ]
	 *					[  0    ]
	 *
	 *	In 3d:
	 *     	   Bvi=[ dh/dx          0             0      ]
	 *					[   0           dh/dy           0      ]
	 *					[   0             0           dh/dz    ]
	 *					[ 1/2*dh/dy    1/2*dh/dx        0      ]
	 *					[ 1/2*dh/dz       0         1/2*dh/dx  ]
	 *					[   0          1/2*dh/dz    1/2*dh/dy  ]
	 *					[   0             0             0      ]
	 *					[ dh/dx         dh/dy         dh/dz    ]
	 *
	 *         Bpi=[ 0 ]
	 *					[ 0 ]
	 *					[ 0 ]
	 *					[ 0 ]
	 *					[ 0 ]
	 *					[ 0 ]
	 *					[ h ]
	 *					[ 0 ]
	 *	where phi is the finiteelement function for node i.
	 *	Same thing for Bb except the last column that does not exist.
	 */

	/*Fetch number of nodes for this finite element*/
	int pnumnodes = element->NumberofNodesPressure();
	int vnumnodes = element->NumberofNodesVelocity();

	/*Get nodal functions derivatives*/
	IssmDouble* vdbasis=xNew<IssmDouble>(dim*vnumnodes);
	IssmDouble* pbasis =xNew<IssmDouble>(pnumnodes);
	element->NodalFunctionsDerivativesVelocity(vdbasis,xyz_list,gauss);
	element->NodalFunctionsPressure(pbasis,gauss);

	/*Build B: */
	if(dim==2){
		for(int i=0;i<vnumnodes;i++){
			B[(dim*vnumnodes+pnumnodes)*0+dim*i+0] = vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*0+dim*i+1] = 0.;
			B[(dim*vnumnodes+pnumnodes)*1+dim*i+0] = 0.;
			B[(dim*vnumnodes+pnumnodes)*1+dim*i+1] = vdbasis[1*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*2+dim*i+0] = .5*vdbasis[1*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*2+dim*i+1] = .5*vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*3+dim*i+0] = 0.;
			B[(dim*vnumnodes+pnumnodes)*3+dim*i+1] = 0.;
			B[(dim*vnumnodes+pnumnodes)*4+dim*i+0] = vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*4+dim*i+1] = vdbasis[1*vnumnodes+i];
		}
		for(int i=0;i<pnumnodes;i++){
			B[(dim*vnumnodes+pnumnodes)*0+(dim*vnumnodes)+i] = 0.;
			B[(dim*vnumnodes+pnumnodes)*1+(dim*vnumnodes)+i] = 0.;
			B[(dim*vnumnodes+pnumnodes)*2+(dim*vnumnodes)+i] = 0.;
			B[(dim*vnumnodes+pnumnodes)*3+(dim*vnumnodes)+i] = pbasis[i];
			B[(dim*vnumnodes+pnumnodes)*4+(dim*vnumnodes)+i] = 0.;
		}
	}
	else{
		for(int i=0;i<vnumnodes;i++){
			B[(dim*vnumnodes+pnumnodes)*0+dim*i+0] = vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*0+dim*i+1] = 0.;
			B[(dim*vnumnodes+pnumnodes)*0+dim*i+2] = 0.;
			B[(dim*vnumnodes+pnumnodes)*1+dim*i+0] = 0.;
			B[(dim*vnumnodes+pnumnodes)*1+dim*i+1] = vdbasis[1*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*1+dim*i+2] = 0.;
			B[(dim*vnumnodes+pnumnodes)*2+dim*i+0] = 0.;
			B[(dim*vnumnodes+pnumnodes)*2+dim*i+1] = 0.;
			B[(dim*vnumnodes+pnumnodes)*2+dim*i+2] = vdbasis[2*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*3+dim*i+0] = .5*vdbasis[1*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*3+dim*i+1] = .5*vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*3+dim*i+2] = 0.;
			B[(dim*vnumnodes+pnumnodes)*4+dim*i+0] = .5*vdbasis[2*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*4+dim*i+1] = 0.;
			B[(dim*vnumnodes+pnumnodes)*4+dim*i+2] = .5*vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*5+dim*i+0] = 0.;
			B[(dim*vnumnodes+pnumnodes)*5+dim*i+1] = .5*vdbasis[2*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*5+dim*i+2] = .5*vdbasis[1*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*6+dim*i+0] = 0.;
			B[(dim*vnumnodes+pnumnodes)*6+dim*i+1] = 0.;
			B[(dim*vnumnodes+pnumnodes)*6+dim*i+2] = 0.;
			B[(dim*vnumnodes+pnumnodes)*7+dim*i+0] = vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*7+dim*i+1] = vdbasis[1*vnumnodes+i];
			B[(dim*vnumnodes+pnumnodes)*7+dim*i+2] = vdbasis[2*vnumnodes+i];
		}
		for(int i=0;i<pnumnodes;i++){
			B[(dim*vnumnodes+pnumnodes)*0+(dim*vnumnodes)+i] = 0.;
			B[(dim*vnumnodes+pnumnodes)*1+(dim*vnumnodes)+i] = 0.;
			B[(dim*vnumnodes+pnumnodes)*2+(dim*vnumnodes)+i] = 0.;
			B[(dim*vnumnodes+pnumnodes)*3+(dim*vnumnodes)+i] = 0.;
			B[(dim*vnumnodes+pnumnodes)*4+(dim*vnumnodes)+i] = 0.;
			B[(dim*vnumnodes+pnumnodes)*5+(dim*vnumnodes)+i] = 0.;
			B[(dim*vnumnodes+pnumnodes)*6+(dim*vnumnodes)+i] = pbasis[i];
			B[(dim*vnumnodes+pnumnodes)*7+(dim*vnumnodes)+i] = 0.;
		}
	}

	/*Clean up*/
	xDelete<IssmDouble>(vdbasis);
	xDelete<IssmDouble>(pbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetBFSFriction(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/* Compute L  matrix. L=[L1 L2 L3] where Li is square and of size numdof.
	 * For node i, Li can be expressed in the actual coordinate system
	 * by in 3d
	 *       Li=[ h 0 0 0 ]
	 *	 	      [ 0 h 0 0 ]
	 *	in 2d:
	 *       Li=[ h 0 0 ]
	 * where h is the interpolation function for node i.
	 */

	/*Fetch number of nodes for this finite element*/
	int pnumnodes = element->NumberofNodesPressure();
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumdof   = pnumnodes;
	int vnumdof   = vnumnodes*dim;

	/*Get nodal functions derivatives*/
	IssmDouble* vbasis=xNew<IssmDouble>(vnumnodes);
	element->NodalFunctionsVelocity(vbasis,gauss);

	/*Build B: */
	if(dim==3){
		for(int i=0;i<vnumnodes;i++){
			B[(vnumdof+pnumdof)*0+3*i+0] = vbasis[i];
			B[(vnumdof+pnumdof)*0+3*i+1] = 0.;
			B[(vnumdof+pnumdof)*0+3*i+2] = 0.;

			B[(vnumdof+pnumdof)*1+3*i+0] = 0.;
			B[(vnumdof+pnumdof)*1+3*i+1] = vbasis[i];
			B[(vnumdof+pnumdof)*1+3*i+2] = 0.;

			B[(vnumdof+pnumdof)*2+3*i+0] = 0.;
			B[(vnumdof+pnumdof)*2+3*i+1] = 0.;
			B[(vnumdof+pnumdof)*2+3*i+2] = vbasis[i];
		}
		for(int i=0;i<pnumnodes;i++){
			B[(vnumdof+pnumdof)*0+i+vnumdof+0] = 0.;
			B[(vnumdof+pnumdof)*1+i+vnumdof+0] = 0.;
			B[(vnumdof+pnumdof)*2+i+vnumdof+0] = 0.;
		}
	}
	else{
		for(int i=0;i<vnumnodes;i++){
			B[(vnumdof+pnumdof)*0+2*i+0] = vbasis[i];
			B[(vnumdof+pnumdof)*0+2*i+1] = 0.;

			B[(vnumdof+pnumdof)*1+2*i+0] = 0.;
			B[(vnumdof+pnumdof)*1+2*i+1] = vbasis[i];
		}

		for(int i=0;i<pnumnodes;i++){
			B[(vnumdof+pnumdof)*0+i+vnumdof+0] = 0.;
			B[(vnumdof+pnumdof)*1+i+vnumdof+0] = 0.;
		}
	}

	/*Clean-up*/
	xDelete<IssmDouble>(vbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetBFSprime(IssmDouble* Bprime,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*	Compute B'  matrix. B'=[B1' B2' B3' B4' B5' B6' Bb'] where Bi' is of size 3*2.
	 *	For node i, Bi' can be expressed in the actual coordinate system
	 *	by:
	 *			Bvi' = [  dphi/dx     0     ]
	 *					 [     0      dphi/dy ]
	 *					 [  dphi/dy   dphi/dx ]
	 *					 [  dphi/dx   dphi/dy ]
	 *					 [     0      0       ]
	 *
	 * by: 	  Bpi=[  0  ]
	 *					[  0  ]
	 *					[  0  ]
	 *					[  0  ]
	 *					[ phi ]
	 *
	 *	In 3d
	 *     	   Bvi=[ dh/dx     0        0    ]
	 *					[   0      dh/dy      0    ]
	 *					[   0        0      dh/dz  ]
	 *					[ dh/dy    dh/dx      0    ]
	 *					[ dh/dz      0      dh/dx  ]
	 *					[   0      dh/dz    dh/dy  ]
	 *					[ dh/dx    dh/dy    dh/dz  ]
	 *					[   0        0        0    ]
	 *
	 *         Bpi=[ 0 ]
	 *					[ 0 ]
	 *					[ 0 ]
	 *					[ 0 ]
	 *					[ 0 ]
	 *					[ 0 ]
	 *					[ 0 ]
	 *					[ h ]
	 *	where phi is the finiteelement function for node i.
	 *	In 3d:
	 */

	/*Fetch number of nodes for this finite element*/
	int pnumnodes = element->NumberofNodesPressure();
	int vnumnodes = element->NumberofNodesVelocity();

	/*Get nodal functions derivatives*/
	IssmDouble* vdbasis=xNew<IssmDouble>(dim*vnumnodes);
	IssmDouble* pbasis =xNew<IssmDouble>(pnumnodes);
	element->NodalFunctionsDerivativesVelocity(vdbasis,xyz_list,gauss);
	element->NodalFunctionsPressure(pbasis,gauss);

	/*Build B_prime: */
	if(dim==2){
		for(int i=0;i<vnumnodes;i++){
			Bprime[(dim*vnumnodes+pnumnodes)*0+dim*i+0] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*0+dim*i+1] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*1+dim*i+0] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*1+dim*i+1] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*2+dim*i+0] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*2+dim*i+1] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*3+dim*i+0] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*3+dim*i+1] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*4+dim*i+0] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*4+dim*i+1] = 0.;
		}
		for(int i=0;i<pnumnodes;i++){
			Bprime[(dim*vnumnodes+pnumnodes)*0+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*1+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*2+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*3+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*4+(dim*vnumnodes)+i] = pbasis[i];
		}
	}
	else{
		for(int i=0;i<vnumnodes;i++){
			Bprime[(dim*vnumnodes+pnumnodes)*0+dim*i+0] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*0+dim*i+1] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*0+dim*i+2] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*1+dim*i+0] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*1+dim*i+1] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*1+dim*i+2] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*2+dim*i+0] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*2+dim*i+1] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*2+dim*i+2] = vdbasis[2*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*3+dim*i+0] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*3+dim*i+1] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*3+dim*i+2] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*4+dim*i+0] = vdbasis[2*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*4+dim*i+1] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*4+dim*i+2] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*5+dim*i+0] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*5+dim*i+1] = vdbasis[2*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*5+dim*i+2] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*6+dim*i+0] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*6+dim*i+1] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*6+dim*i+2] = vdbasis[2*vnumnodes+i];
			Bprime[(dim*vnumnodes+pnumnodes)*7+dim*i+0] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*7+dim*i+1] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*7+dim*i+2] = 0.;
		}
		for(int i=0;i<pnumnodes;i++){
			Bprime[(dim*vnumnodes+pnumnodes)*0+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*1+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*2+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*3+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*4+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*5+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*6+(dim*vnumnodes)+i] = 0.;
			Bprime[(dim*vnumnodes+pnumnodes)*7+(dim*vnumnodes)+i] = pbasis[i];
		}
	}

	/*Clean up*/
	xDelete<IssmDouble>(vdbasis);
	xDelete<IssmDouble>(pbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetBFSprimeUzawa(IssmDouble* Bprime,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*	Compute B'  matrix. B'=[B1' B2' B3' B4' B5' B6'] where Bi' is of size 3*2.
	 *	For node i, Bi' can be expressed in the actual coordinate system
	 *	by:
	 *			Bvi' = [  dphi/dx   dphi/dy ]
	 *
	 *	In 3d
	 *     	   Bvi=[ dh/dx   dh/dy    dh/dz  ]
	 *	where phi is the finiteelement function for node i.
	 */

	/*Fetch number of nodes for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();

	/*Get nodal functions derivatives*/
	IssmDouble* vdbasis=xNew<IssmDouble>(dim*vnumnodes);
	element->NodalFunctionsDerivativesVelocity(vdbasis,xyz_list,gauss);

	/*Build B_prime: */
	if(dim==2){
		for(int i=0;i<vnumnodes;i++){
			Bprime[dim*i+0] = vdbasis[0*vnumnodes+i];
			Bprime[dim*i+1] = vdbasis[1*vnumnodes+i];
		}
	}
	else{
		for(int i=0;i<vnumnodes;i++){
			Bprime[dim*i+0] = vdbasis[0*vnumnodes+i];
			Bprime[dim*i+1] = vdbasis[1*vnumnodes+i];
			Bprime[dim*i+2] = vdbasis[2*vnumnodes+i];
		}
	}

	/*Clean up*/
	xDelete<IssmDouble>(vdbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetBFSprimevel(IssmDouble* Bprime,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*	Compute B'  matrix. B'=[B1' B2' B3' B4' B5' B6' Bb'] where Bi' is of size 3*2.
	 *	For node i, Bi' can be expressed in the actual coordinate system
	 *	by:
	 *			Bvi' = [  dphi/dx     0     ]
	 *					 [     0      dphi/dy ]
	 *					 [  dphi/dy   dphi/dx ]
	 *
	 *	In 3d
	 *     	   Bvi=[ dh/dx     0        0    ]
	 *					[   0      dh/dy      0    ]
	 *					[   0        0      dh/dz  ]
	 *					[ dh/dy    dh/dx      0    ]
	 *					[ dh/dz      0      dh/dx  ]
	 *					[   0      dh/dz    dh/dy  ]
	 *	where phi is the finiteelement function for node i.
	 *	In 3d:
	 */

	/*Fetch number of nodes for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();

	/*Get nodal functions derivatives*/
	IssmDouble* vdbasis=xNew<IssmDouble>(dim*vnumnodes);
	element->NodalFunctionsDerivativesVelocity(vdbasis,xyz_list,gauss);

	/*Build B_prime: */
	if(dim==2){
		for(int i=0;i<vnumnodes;i++){
			Bprime[(dim*vnumnodes)*0+dim*i+0] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes)*0+dim*i+1] = 0.;
			Bprime[(dim*vnumnodes)*1+dim*i+0] = 0.;
			Bprime[(dim*vnumnodes)*1+dim*i+1] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes)*2+dim*i+0] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes)*2+dim*i+1] = vdbasis[0*vnumnodes+i];
		}
	}
	else{
		for(int i=0;i<vnumnodes;i++){
			Bprime[(dim*vnumnodes)*0+dim*i+0] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes)*0+dim*i+1] = 0.;
			Bprime[(dim*vnumnodes)*0+dim*i+2] = 0.;
			Bprime[(dim*vnumnodes)*1+dim*i+0] = 0.;
			Bprime[(dim*vnumnodes)*1+dim*i+1] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes)*1+dim*i+2] = 0.;
			Bprime[(dim*vnumnodes)*2+dim*i+0] = 0.;
			Bprime[(dim*vnumnodes)*2+dim*i+1] = 0.;
			Bprime[(dim*vnumnodes)*2+dim*i+2] = vdbasis[2*vnumnodes+i];
			Bprime[(dim*vnumnodes)*3+dim*i+0] = vdbasis[1*vnumnodes+i];
			Bprime[(dim*vnumnodes)*3+dim*i+1] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes)*3+dim*i+2] = 0.;
			Bprime[(dim*vnumnodes)*4+dim*i+0] = vdbasis[2*vnumnodes+i];
			Bprime[(dim*vnumnodes)*4+dim*i+1] = 0.;
			Bprime[(dim*vnumnodes)*4+dim*i+2] = vdbasis[0*vnumnodes+i];
			Bprime[(dim*vnumnodes)*5+dim*i+0] = 0.;
			Bprime[(dim*vnumnodes)*5+dim*i+1] = vdbasis[2*vnumnodes+i];
			Bprime[(dim*vnumnodes)*5+dim*i+2] = vdbasis[1*vnumnodes+i];
		}
	}

	/*Clean up*/
	xDelete<IssmDouble>(vdbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetBFSUzawa(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[Bp1 Bp2 ...] where Bpi=phi_pi.
	 */

	/*Fetch number of nodes for this finite element*/
	int pnumnodes;
	if(dim==2) pnumnodes=3;
	else pnumnodes=6;
	//int pnumnodes = element->NumberofNodes(P1Enum);

	/*Get nodal functions derivatives*/
	IssmDouble* basis =xNew<IssmDouble>(pnumnodes);
	element->NodalFunctionsP1(basis,gauss);

	/*Build B: */
	for(int i=0;i<pnumnodes;i++){
		B[i] = basis[i];
	}

	/*Clean up*/
	xDelete<IssmDouble>(basis);
}/*}}}*/
void           StressbalanceAnalysis::GetBFSvel(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[Bv1 Bv2 ... Bp1 Bp2 ...] where Bvi is of size 3*3.
	 * For node i, Bvi can be expressed in the actual coordinate system
	 * by: 	   Bvi=[ dphi/dx          0        ]
	 *					 [   0           dphi/dy     ]
	 *					 [ 1/2*dphi/dy    1/2*dphi/dx]
	 *
	 *
	 *	In 3d:
	 *     	   Bvi=[ dh/dx          0             0      ]
	 *					[   0           dh/dy           0      ]
	 *					[   0             0           dh/dz    ]
	 *					[ 1/2*dh/dy    1/2*dh/dx        0      ]
	 *					[ 1/2*dh/dz       0         1/2*dh/dx  ]
	 *					[   0          1/2*dh/dz    1/2*dh/dy  ]
	 *
	 *	where phi is the finiteelement function for node i.
	 *	Same thing for Bb except the last column that does not exist.
	 */

	/*Fetch number of nodes for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();

	/*Get nodal functions derivatives*/
	IssmDouble* vdbasis=xNew<IssmDouble>(dim*vnumnodes);
	element->NodalFunctionsDerivativesVelocity(vdbasis,xyz_list,gauss);

	/*Build B: */
	if(dim==2){
		for(int i=0;i<vnumnodes;i++){
			B[(dim*vnumnodes)*0+dim*i+0] = vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes)*0+dim*i+1] = 0.;
			B[(dim*vnumnodes)*1+dim*i+0] = 0.;
			B[(dim*vnumnodes)*1+dim*i+1] = vdbasis[1*vnumnodes+i];
			B[(dim*vnumnodes)*2+dim*i+0] = .5*vdbasis[1*vnumnodes+i];
			B[(dim*vnumnodes)*2+dim*i+1] = .5*vdbasis[0*vnumnodes+i];
		}
	}
	else{
		for(int i=0;i<vnumnodes;i++){
			B[(dim*vnumnodes)*0+dim*i+0] = vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes)*0+dim*i+1] = 0.;
			B[(dim*vnumnodes)*0+dim*i+2] = 0.;
			B[(dim*vnumnodes)*1+dim*i+0] = 0.;
			B[(dim*vnumnodes)*1+dim*i+1] = vdbasis[1*vnumnodes+i];
			B[(dim*vnumnodes)*1+dim*i+2] = 0.;
			B[(dim*vnumnodes)*2+dim*i+0] = 0.;
			B[(dim*vnumnodes)*2+dim*i+1] = 0.;
			B[(dim*vnumnodes)*2+dim*i+2] = vdbasis[2*vnumnodes+i];
			B[(dim*vnumnodes)*3+dim*i+0] = .5*vdbasis[1*vnumnodes+i];
			B[(dim*vnumnodes)*3+dim*i+1] = .5*vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes)*3+dim*i+2] = 0.;
			B[(dim*vnumnodes)*4+dim*i+0] = .5*vdbasis[2*vnumnodes+i];
			B[(dim*vnumnodes)*4+dim*i+1] = 0.;
			B[(dim*vnumnodes)*4+dim*i+2] = .5*vdbasis[0*vnumnodes+i];
			B[(dim*vnumnodes)*5+dim*i+0] = 0.;
			B[(dim*vnumnodes)*5+dim*i+1] = .5*vdbasis[2*vnumnodes+i];
			B[(dim*vnumnodes)*5+dim*i+2] = .5*vdbasis[1*vnumnodes+i];
		}
	}

	/*Clean up*/
	xDelete<IssmDouble>(vdbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetCFS(IssmDouble* C,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute C  matrix. C=[Cp1 Cp2 ...] where:
	 *     Cpi=[phi phi].
	 */

	/*Fetch number of nodes for this finite element*/
	int lnumnodes = element->GetNumberOfNodes(P2Enum);

	/*Get nodal functions derivatives*/
	IssmDouble* basis =xNew<IssmDouble>(lnumnodes);
	element->NodalFunctionsP2(basis,gauss);

	/*Build B: */
	for(int i=0;i<lnumnodes;i++){
		C[lnumnodes*0+i] = basis[i];
		C[lnumnodes*1+i] = basis[i];
	}

	/*Clean up*/
	xDelete<IssmDouble>(basis);
}/*}}}*/
void           StressbalanceAnalysis::GetCFSprime(IssmDouble* Cprime,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*	Compute C'  matrix. C'=[C1' C2' ...]
	 *			Ci' = [  phi  0  ]
	 *			      [   0  phi ]
	 *
	 *	In 3d
	 *			Ci' = [  phi  0   0  ]
	 *			      [   0  phi  0  ]
	 *			      [   0   0  phi ]
	 *	where phi is the finiteelement function for node i.
	 */

	/*Fetch number of nodes for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int vnumdof   = vnumnodes*dim;

	IssmDouble* vbasis=xNew<IssmDouble>(vnumnodes);
	element->NodalFunctionsVelocity(vbasis,gauss);

	/*Build B: */
	if(dim==3){
		for(int i=0;i<vnumnodes;i++){
			Cprime[vnumdof*0+3*i+0] = vbasis[i];
			Cprime[vnumdof*0+3*i+1] = 0.;
			Cprime[vnumdof*0+3*i+2] = 0.;

			Cprime[vnumdof*1+3*i+0] = 0.;
			Cprime[vnumdof*1+3*i+1] = vbasis[i];
			Cprime[vnumdof*1+3*i+2] = 0.;

			Cprime[vnumdof*2+3*i+0] = 0.;
			Cprime[vnumdof*2+3*i+1] = 0.;
			Cprime[vnumdof*2+3*i+2] = vbasis[i];
		}
	}
	else{
		for(int i=0;i<vnumnodes;i++){
			Cprime[vnumdof*0+2*i+0] = vbasis[i];
			Cprime[vnumdof*0+2*i+1] = 0.;

			Cprime[vnumdof*1+2*i+0] = 0.;
			Cprime[vnumdof*1+2*i+1] = vbasis[i];
		}
	}

	/*Clean-up*/
	xDelete<IssmDouble>(vbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetSolutionFromInputsFS(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	int*         vdoflist=NULL;
	int*         pdoflist=NULL;
	Input*       vz_input=NULL;
	int          dim;
	IssmDouble   vx,vy,vz,p;
	IssmDouble   FSreconditioning;

	/*Get some parameters*/
	element->FindParam(&dim,DomainDimensionEnum);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int vnumdof   = vnumnodes*dim;
	int pnumdof   = pnumnodes*1;

	/*Initialize values*/
	IssmDouble* vvalues = xNew<IssmDouble>(vnumdof);
	IssmDouble* pvalues = xNew<IssmDouble>(pnumdof);

	/*Get dof list: */
	element->GetDofListVelocity(&vdoflist,GsetEnum);
	element->GetDofListPressure(&pdoflist,GsetEnum);
	Input*     vx_input=element->GetInput(VxEnum);       _assert_(vx_input);
	Input*     vy_input=element->GetInput(VyEnum);       _assert_(vy_input);
	if(dim==3){vz_input=element->GetInput(VzEnum);       _assert_(vz_input);}
	Input*     p_input =element->GetInput(PressureEnum); _assert_(p_input);

	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);

	/*Ok, we have the velocities in inputs, fill in solution */
	Gauss* gauss = element->NewGauss();
	for(int i=0;i<vnumnodes;i++){
		gauss->GaussNode(element->VelocityInterpolation(),i);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vvalues[i*dim+0]=vx;
		vvalues[i*dim+1]=vy;
		if(dim==3){
			vz_input->GetInputValue(&vz,gauss);
			vvalues[i*dim+2]=vz;
		}
	}
	for(int i=0;i<pnumnodes;i++){
		gauss->GaussNode(element->PressureInterpolation(),i);
		p_input->GetInputValue(&p ,gauss);
		pvalues[i]=p/FSreconditioning;
	}

	/*Do NOT account for bubble using GetInputValue! This is wrong*/
	if(element->VelocityInterpolation()==P1bubblecondensedEnum ||
				element->VelocityInterpolation()==P1bubbleEnum){
		vvalues[(vnumnodes-1)*dim+0]=0.;
		vvalues[(vnumnodes-1)*dim+1]=0.;
		if(dim==3) vvalues[(vnumnodes-1)*dim+2]=0.;
	}

	/*Add value to global vector*/
	solution->SetValues(vnumdof,vdoflist,vvalues,INS_VAL);
	if(pnumdof>0) solution->SetValues(pnumdof,pdoflist,pvalues,INS_VAL);

	/*Free resources:*/
	delete gauss;
	xDelete<int>(pdoflist);
	xDelete<int>(vdoflist);
	xDelete<IssmDouble>(pvalues);
	xDelete<IssmDouble>(vvalues);
}/*}}}*/
void           StressbalanceAnalysis::InitializeXTH(Elements* elements,Parameters* parameters){/*{{{*/

	/*Intermediaries*/
	int        dim;
	IssmDouble dvx[3],dvy[3],dvz[3];
	IssmDouble viscosity;
	IssmDouble *xyz_list = NULL;

	/*Get problem dimension*/
	parameters->FindParam(&dim,DomainDimensionEnum);

	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);

		/*Get inputs and parameters*/
		element->GetVerticesCoordinates(&xyz_list);
		Input* vx_input=element->GetInput(VxEnum); _assert_(vx_input);
		Input* vy_input=element->GetInput(VyEnum); _assert_(vy_input);
		Input* vz_input;
		if(dim==3){vz_input=element->GetInput(VzEnum); _assert_(vz_input);}

		/*Allocate new inputs*/
		int tnumnodes = element->GetNumberOfVertices();      //Tensors, P1 DG
		IssmDouble* epsxx = xNew<IssmDouble>(tnumnodes); IssmDouble* sigmapxx = xNew<IssmDouble>(tnumnodes);
		IssmDouble* epsyy = xNew<IssmDouble>(tnumnodes); IssmDouble* sigmapyy = xNew<IssmDouble>(tnumnodes);
		IssmDouble* epsxy = xNew<IssmDouble>(tnumnodes); IssmDouble* sigmapxy = xNew<IssmDouble>(tnumnodes);
		IssmDouble* epszz = NULL;                        IssmDouble* sigmapzz = NULL;
		IssmDouble* epsxz = NULL;                        IssmDouble* sigmapxz = NULL;
		IssmDouble* epsyz = NULL;                        IssmDouble* sigmapyz = NULL;
		if(dim==3){
			epszz = xNew<IssmDouble>(tnumnodes); sigmapzz = xNew<IssmDouble>(tnumnodes);
			epsxz = xNew<IssmDouble>(tnumnodes); sigmapxz = xNew<IssmDouble>(tnumnodes);
			epsyz = xNew<IssmDouble>(tnumnodes); sigmapyz = xNew<IssmDouble>(tnumnodes);
		}

		/*Get d and tau*/
		Gauss* gauss = element->NewGauss();
		for(int i=0;i<tnumnodes;i++){
			gauss->GaussNode(P1DGEnum,i);

			vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
			vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
			if(dim==3){
				vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
			}

			element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);
			epsxx[i] = dvx[0];                sigmapxx[i] = 2.*viscosity*epsxx[i];
			epsyy[i] = dvy[1];                sigmapyy[i] = 2.*viscosity*epsyy[i];
			epsxy[i] = 0.5*(dvx[1] + dvy[0]); sigmapxy[i] = 2.*viscosity*epsxy[i];
			if(dim==3){
				epszz[i] = dvz[2];                sigmapzz[i] = 2.*viscosity*epszz[i];
				epsxz[i] = 0.5*(dvx[2] + dvz[0]); sigmapxz[i] = 2.*viscosity*epsxz[i];
				epsyz[i] = 0.5*(dvy[2] + dvz[1]); sigmapyz[i] = 2.*viscosity*epsyz[i];
			}
		}

		/*Add inputs*/
		element->AddInput(StrainRatexxEnum,epsxx,P1DGEnum); element->AddInput(DeviatoricStressxxEnum,sigmapxx,P1DGEnum);
		element->AddInput(StrainRateyyEnum,epsyy,P1DGEnum); element->AddInput(DeviatoricStressyyEnum,sigmapyy,P1DGEnum);
		element->AddInput(StrainRatexyEnum,epsxy,P1DGEnum); element->AddInput(DeviatoricStressxyEnum,sigmapxy,P1DGEnum);
		if(dim==3){
			element->AddInput(StrainRatezzEnum,epszz,P1DGEnum); element->AddInput(DeviatoricStresszzEnum,sigmapzz,P1DGEnum);
			element->AddInput(StrainRatexzEnum,epsxz,P1DGEnum); element->AddInput(DeviatoricStressxzEnum,sigmapxz,P1DGEnum);
			element->AddInput(StrainRateyzEnum,epsyz,P1DGEnum); element->AddInput(DeviatoricStressyzEnum,sigmapyz,P1DGEnum);
		}

		/*Clean up*/
		delete gauss;
		xDelete<IssmDouble>(xyz_list);
		xDelete<IssmDouble>(epsxx); xDelete<IssmDouble>(sigmapxx);
		xDelete<IssmDouble>(epsyy); xDelete<IssmDouble>(sigmapyy);
		xDelete<IssmDouble>(epszz); xDelete<IssmDouble>(sigmapzz);
		xDelete<IssmDouble>(epsxy); xDelete<IssmDouble>(sigmapxy);
		xDelete<IssmDouble>(epsxz); xDelete<IssmDouble>(sigmapxz);
		xDelete<IssmDouble>(epsyz); xDelete<IssmDouble>(sigmapyz);
	}

}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolutionFS(IssmDouble* solution,Element* element){/*{{{*/

	int          i,dim;
	int*         vdoflist=NULL;
	int*         pdoflist=NULL;
	IssmDouble   FSreconditioning;

	element->FindParam(&dim,DomainDimensionEnum);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);

	/*Fetch number of nodes and dof for this finite element*/
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int vnumdof   = vnumnodes*dim;
	int pnumdof   = pnumnodes*1;

	/*Initialize values*/
	IssmDouble* values   = xNew<IssmDouble>(vnumdof+pnumdof);
	IssmDouble* vx       = xNew<IssmDouble>(vnumnodes);
	IssmDouble* vy       = xNew<IssmDouble>(vnumnodes);
	IssmDouble* vz       = xNew<IssmDouble>(vnumnodes);
	IssmDouble* vel      = xNew<IssmDouble>(vnumnodes);
	IssmDouble* pressure = xNew<IssmDouble>(pnumnodes);

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	if(dim==2) for(i=0;i<vnumnodes;i++) cs_list[i] = XYEnum;
	else       for(i=0;i<vnumnodes;i++) cs_list[i] = XYZEnum;
	for(i=0;i<pnumnodes;i++) cs_list[vnumnodes+i] = PressureEnum;

	/*Get dof list: */
	element->GetDofListLocalVelocity(&vdoflist,GsetEnum);
	element->GetDofListLocalPressure(&pdoflist,GsetEnum);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<vnumdof;i++) values[i]        =solution[vdoflist[i]];
	for(i=0;i<pnumdof;i++) values[vnumdof+i]=solution[pdoflist[i]];

	/*Transform solution in Cartesian Space*/
	element->TransformSolutionCoord(values,cs_list);

	/*Ok, we have vx and vy in values, fill in all arrays: */
	for(i=0;i<vnumnodes;i++){
		vx[i] = values[i*dim+0];
		vy[i] = values[i*dim+1];
		if(xIsNan<IssmDouble>(vx[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vx[i])) _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(vy[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vy[i])) _error_("Inf found in solution vector");

		if(dim==3){
			vz[i] = values[i*dim+2];
			if(xIsNan<IssmDouble>(vz[i])) _error_("NaN found in solution vector");
			if(xIsInf<IssmDouble>(vz[i])) _error_("Inf found in solution vector");
		}
	}
	for(i=0;i<pnumnodes;i++){
		pressure[i] = values[vnumdof+i];
		if(xIsNan<IssmDouble>(pressure[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(pressure[i])) _error_("Inf found in solution vector");
	}

	/*Recondition pressure and compute vel: */
	for(i=0;i<pnumnodes;i++) pressure[i] = pressure[i]*FSreconditioning;
	if(dim==3) for(i=0;i<vnumnodes;i++) vel[i] = sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	else       for(i=0;i<vnumnodes;i++) vel[i] = sqrt(vx[i]*vx[i] + vy[i]*vy[i]);

	/*Add vx and vy as inputs to the tria element: */
	int v_interp =  element->VelocityInterpolation();
	if(v_interp==P1bubbleEnum) v_interp=P1Enum;
	if(v_interp == P1bubblecondensedEnum) v_interp = P1Enum;
	element->AddInput(VxEnum, vx, v_interp);
	element->AddInput(VyEnum, vy, v_interp);
	element->AddInput(VelEnum,vel,v_interp);
	if(pnumdof>0) element->AddInput(PressureEnum,pressure,element->PressureInterpolation());
	if(dim==3) element->AddInput(VzEnum,vz,v_interp);

	/*Free resources:*/
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vz);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(values);
	xDelete<int>(vdoflist);
	xDelete<int>(pdoflist);
	xDelete<int>(cs_list);
}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolutionFSXTH_d(Elements* elements,Parameters* parameters){/*{{{*/

	/*Intermediaries*/
	int         dim,tausize;
	IssmDouble  epsxx,epsyy,epszz,epsxy,epsxz,epsyz,D_scalar;
	IssmDouble  epsxx_old,epsyy_old,epszz_old,epsxy_old,epsxz_old,epsyz_old;
	IssmDouble  sigmapxx,sigmapyy,sigmapzz,sigmapxy,sigmapxz,sigmapyz;
	IssmDouble  dvx[3],dvy[3],dvz[3],B,n;
	IssmDouble *xyz_list = NULL;
	IssmDouble  Jdet,r;

	parameters->FindParam(&r,AugmentedLagrangianREnum);
	parameters->FindParam(&dim,DomainDimensionEnum);
	if(dim==2) tausize = 3;
	else       tausize = 6;

	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);

		/*Get inputs and parameters*/
		element->GetVerticesCoordinates(&xyz_list);
		Input*  B_input=element->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
		Input*  n_input=element->GetInput(MaterialsRheologyNEnum); _assert_(n_input);
		Input* vx_input=element->GetInput(VxEnum);                 _assert_(vx_input);
		Input* vy_input=element->GetInput(VyEnum);                 _assert_(vy_input);
		Input* vz_input;
		if(dim==3){vz_input=element->GetInput(VzEnum); _assert_(vz_input);}

		/*Fetch number of nodes and dof for this finite element*/
		int tnumnodes = element->GetNumberOfVertices();      //Tensors, P1 DG

		/*Initialize vectors*/
		IssmDouble* tbasis = xNew<IssmDouble>(tnumnodes);
		IssmDouble* Ke     = xNewZeroInit<IssmDouble>(tnumnodes*tnumnodes);
		IssmDouble* pe_xx  = xNewZeroInit<IssmDouble>(tnumnodes);
		IssmDouble* pe_yy  = xNewZeroInit<IssmDouble>(tnumnodes);
		IssmDouble* pe_xy  = xNewZeroInit<IssmDouble>(tnumnodes);
		IssmDouble* pe_zz  = NULL; IssmDouble* pe_xz  = NULL; IssmDouble* pe_yz  = NULL;
		if(dim==3){
			pe_zz = xNewZeroInit<IssmDouble>(tnumnodes);
			pe_xz = xNewZeroInit<IssmDouble>(tnumnodes);
			pe_yz = xNewZeroInit<IssmDouble>(tnumnodes);
		}

		/*Get previous d*/
		Input* epsxx_input=element->GetInput(StrainRatexxEnum); _assert_(epsxx_input);
		Input* epsyy_input=element->GetInput(StrainRateyyEnum); _assert_(epsyy_input);
		Input* epsxy_input=element->GetInput(StrainRatexyEnum); _assert_(epsxy_input);
		Input* epszz_input=NULL; Input* epsxz_input=NULL; Input* epsyz_input=NULL;
		if(dim==3){
			epszz_input=element->GetInput(StrainRatezzEnum); _assert_(epszz_input);
			epsxz_input=element->GetInput(StrainRatexzEnum); _assert_(epsxz_input);
			epsyz_input=element->GetInput(StrainRateyzEnum); _assert_(epsyz_input);
		}

		/*Get tau*/
		Input* sigmapxx_input=element->GetInput(DeviatoricStressxxEnum); _assert_(sigmapxx_input);
		Input* sigmapyy_input=element->GetInput(DeviatoricStressyyEnum); _assert_(sigmapyy_input);
		Input* sigmapxy_input=element->GetInput(DeviatoricStressxyEnum); _assert_(sigmapxy_input);
		Input* sigmapzz_input=NULL; Input* sigmapxz_input=NULL; Input* sigmapyz_input=NULL;
		if(dim==3){
			sigmapzz_input=element->GetInput(DeviatoricStresszzEnum); _assert_(sigmapzz_input);
			sigmapxz_input=element->GetInput(DeviatoricStressxzEnum); _assert_(sigmapxz_input);
			sigmapyz_input=element->GetInput(DeviatoricStressyzEnum); _assert_(sigmapyz_input);
		}

		Gauss* gauss=element->NewGauss(5);
		while(gauss->next()){
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);
			element->NodalFunctionsTensor(tbasis,gauss);

			/*Get tau from inputs*/
			sigmapxx_input->GetInputValue(&sigmapxx,gauss);
			sigmapyy_input->GetInputValue(&sigmapyy,gauss);
			sigmapxy_input->GetInputValue(&sigmapxy,gauss);
			if(dim==3){
				sigmapzz_input->GetInputValue(&sigmapzz,gauss);
				sigmapxz_input->GetInputValue(&sigmapxz,gauss);
				sigmapyz_input->GetInputValue(&sigmapyz,gauss);
			}

			/*Get previous d*/
			epsxx_input->GetInputValue(&epsxx_old,gauss);
			epsyy_input->GetInputValue(&epsyy_old,gauss);
			epsxy_input->GetInputValue(&epsxy_old,gauss);
			if(dim==3){
				epszz_input->GetInputValue(&epszz_old,gauss);
				epsxz_input->GetInputValue(&epsxz_old,gauss);
				epsyz_input->GetInputValue(&epsyz_old,gauss);
			}

			/*Calculate d from previous results*/
			vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
			vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
			if(dim==3){
				vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
			}
			epsxx = dvx[0];
			epsyy = dvy[1];
			epsxy = 0.5*(dvx[1] + dvy[0]);
			if(dim==3){
				epszz = dvz[2];
				epsxz = 0.5*(dvx[2] + dvz[0]);
				epsyz = 0.5*(dvy[2] + dvz[1]);
			}

			/*Solve 2 eta_0 |d|^s-1 + r |d| = |rD(u) + tau|*/
			IssmDouble coef1,coef2,coef3;
			B_input->GetInputValue(&B,gauss);
			n_input->GetInputValue(&n,gauss);
			coef1 = B*pow(1./sqrt(2.),(1.-n)/n); //2 eta_0 = 2 * B/(2* (1/sqrt(2)  )^(n-1)/n )
			coef2 = r;
			if(dim==2){
				coef3 = sqrt(
							  (r*epsxx + sigmapxx)*(r*epsxx + sigmapxx)
							+ (r*epsyy + sigmapyy)*(r*epsyy + sigmapyy)
							+ 2*(r*epsxy + sigmapxy)*(r*epsxy + sigmapxy)
							);
			}
			else{
				coef3 = sqrt(
					  		  (r*epsxx + sigmapxx)*(r*epsxx + sigmapxx)
							+ (r*epsyy + sigmapyy)*(r*epsyy + sigmapyy)
							+ (r*epszz + sigmapzz)*(r*epszz + sigmapzz)
							+ 2*(r*epsxy + sigmapxy)*(r*epsxy + sigmapxy)
							+ 2*(r*epsxz + sigmapxz)*(r*epsxz + sigmapxz)
							+ 2*(r*epsyz + sigmapyz)*(r*epsyz + sigmapyz)
							);
			}
			IssmDouble dnorm;
			if(dim==2){
				dnorm = sqrt( epsxx_old*epsxx_old + epsyy_old*epsyy_old + 2.*epsxy_old*epsxy_old );
			}
			else{
				dnorm = sqrt( epsxx_old*epsxx_old + epsyy_old*epsyy_old + epszz_old*epszz_old
							+2.*(epsxy_old*epsxy_old + epsxz_old*epsxz_old + epsyz_old*epsyz_old));
			}
			/*Initial guess cannot be 0 otherwise log(0)  - inf*/
			if(dnorm==0.) dnorm=1.;
			NewtonSolveDnorm(&dnorm,coef1,coef2,coef3,n,dnorm);
			_assert_(dnorm>=0.);
			_assert_(!xIsNan<IssmDouble>(dnorm));

			/*Create Ke*/
			D_scalar=(coef1*pow(dnorm,(1.-n)/n)+r)*gauss->weight*Jdet;
			TripleMultiply(tbasis,tnumnodes,1,0,
						&D_scalar,1,1,0,
						tbasis,1,tnumnodes,0,
						Ke,1);

			IssmDouble factor = gauss->weight*Jdet;
			/*Create Right hand sides*/
			for(int ii=0;ii<tnumnodes;ii++) pe_xx[ii] += (r*epsxx+sigmapxx)*tbasis[ii]*factor;
			for(int ii=0;ii<tnumnodes;ii++) pe_yy[ii] += (r*epsyy+sigmapyy)*tbasis[ii]*factor;
			for(int ii=0;ii<tnumnodes;ii++) pe_xy[ii] += (r*epsxy+sigmapxy)*tbasis[ii]*factor;
			if(dim==3){
				for(int ii=0;ii<tnumnodes;ii++) pe_zz[ii] += (r*epszz+sigmapzz)*tbasis[ii]*factor;
				for(int ii=0;ii<tnumnodes;ii++) pe_xz[ii] += (r*epsxz+sigmapxz)*tbasis[ii]*factor;
				for(int ii=0;ii<tnumnodes;ii++) pe_yz[ii] += (r*epsyz+sigmapyz)*tbasis[ii]*factor;
			}
		}

		/*Solve the systems*/
		IssmDouble* d_xx = xNew<IssmDouble>(tnumnodes);
		IssmDouble* d_yy = xNew<IssmDouble>(tnumnodes);
		IssmDouble* d_xy = xNew<IssmDouble>(tnumnodes);
		IssmDouble* d_zz = NULL;
		IssmDouble* d_xz = NULL;
		IssmDouble* d_yz = NULL;
		if(dim==2){
			_assert_(tnumnodes==3);
			Matrix3x3Solve(&d_xx[0],Ke,pe_xx);
			Matrix3x3Solve(&d_yy[0],Ke,pe_yy);
			Matrix3x3Solve(&d_xy[0],Ke,pe_xy);
			for(int i=0;i<3;i++) _assert_(!xIsNan<IssmDouble>(d_xx[i]));
			for(int i=0;i<3;i++) _assert_(!xIsNan<IssmDouble>(d_yy[i]));
			for(int i=0;i<3;i++) _assert_(!xIsNan<IssmDouble>(d_xx[i]));
			element->AddInput(StrainRatexxEnum,d_xx,P1DGEnum);
			element->AddInput(StrainRateyyEnum,d_yy,P1DGEnum);
			element->AddInput(StrainRatexyEnum,d_xy,P1DGEnum);
		}
		else{
			_assert_(tnumnodes==4);
			d_zz = xNew<IssmDouble>(tnumnodes);
			d_xz = xNew<IssmDouble>(tnumnodes);
			d_yz = xNew<IssmDouble>(tnumnodes);
			Matrix4x4Solve(&d_xx[0],Ke,pe_xx);
			Matrix4x4Solve(&d_yy[0],Ke,pe_yy);
			Matrix4x4Solve(&d_xy[0],Ke,pe_xy);
			Matrix4x4Solve(&d_zz[0],Ke,pe_zz);
			Matrix4x4Solve(&d_xz[0],Ke,pe_xz);
			Matrix4x4Solve(&d_yz[0],Ke,pe_yz);
			element->AddInput(StrainRatexxEnum,d_xx,P1DGEnum);
			element->AddInput(StrainRateyyEnum,d_yy,P1DGEnum);
			element->AddInput(StrainRatexyEnum,d_xy,P1DGEnum);
			element->AddInput(StrainRatezzEnum,d_zz,P1DGEnum);
			element->AddInput(StrainRatexzEnum,d_xz,P1DGEnum);
			element->AddInput(StrainRateyzEnum,d_yz,P1DGEnum);
		}

		/*Clean up*/
		delete gauss;
		xDelete<IssmDouble>(xyz_list);
		xDelete<IssmDouble>(tbasis);
		xDelete<IssmDouble>(Ke);
		xDelete<IssmDouble>(pe_xx); xDelete<IssmDouble>(d_xx);
		xDelete<IssmDouble>(pe_yy); xDelete<IssmDouble>(d_yy);
		xDelete<IssmDouble>(pe_zz); xDelete<IssmDouble>(d_zz);
		xDelete<IssmDouble>(pe_xy); xDelete<IssmDouble>(d_xy);
		xDelete<IssmDouble>(pe_xz); xDelete<IssmDouble>(d_xz);
		xDelete<IssmDouble>(pe_yz); xDelete<IssmDouble>(d_yz);
	}
}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolutionFSXTH_tau(Elements* elements,Parameters* parameters){/*{{{*/

	/*Intermediaries*/
	int         dim,tausize;
	IssmDouble  epsxx,epsyy,epszz,epsxy,epsxz,epsyz,D_scalar;
	IssmDouble  d_xx,d_yy,d_zz,d_xy,d_xz,d_yz;
	IssmDouble  sigmapxx,sigmapyy,sigmapzz,sigmapxy,sigmapxz,sigmapyz;
	IssmDouble  dvx[3],dvy[3],dvz[3];
	IssmDouble *xyz_list = NULL;
	IssmDouble  Jdet,r;

	parameters->FindParam(&r,AugmentedLagrangianREnum);
	parameters->FindParam(&dim,DomainDimensionEnum);
	if(dim==2) tausize = 3;
	else       tausize = 6;

	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);

		/*Get inputs and parameters*/
		element->GetVerticesCoordinates(&xyz_list);
		Input* vx_input=element->GetInput(VxEnum);                 _assert_(vx_input);
		Input* vy_input=element->GetInput(VyEnum);                 _assert_(vy_input);
		Input* vz_input=NULL;
		if(dim==3){vz_input=element->GetInput(VzEnum); _assert_(vz_input);}

		/*Get previous tau*/
		Input* sigmapxx_input=element->GetInput(DeviatoricStressxxEnum); _assert_(sigmapxx_input);
		Input* sigmapyy_input=element->GetInput(DeviatoricStressyyEnum); _assert_(sigmapyy_input);
		Input* sigmapxy_input=element->GetInput(DeviatoricStressxyEnum); _assert_(sigmapxy_input);
		Input* sigmapzz_input=NULL; Input* sigmapxz_input=NULL; Input* sigmapyz_input=NULL;
		if(dim==3){
			sigmapzz_input=element->GetInput(DeviatoricStresszzEnum); _assert_(sigmapzz_input);
			sigmapxz_input=element->GetInput(DeviatoricStressxzEnum); _assert_(sigmapxz_input);
			sigmapyz_input=element->GetInput(DeviatoricStressyzEnum); _assert_(sigmapyz_input);
		}

		/*Get NEW d*/
		Input* epsxx_input=element->GetInput(StrainRatexxEnum); _assert_(epsxx_input);
		Input* epsyy_input=element->GetInput(StrainRateyyEnum); _assert_(epsyy_input);
		Input* epsxy_input=element->GetInput(StrainRatexyEnum); _assert_(epsxy_input);
		Input* epszz_input=NULL; Input* epsxz_input=NULL; Input* epsyz_input=NULL;
		if(dim==3){
			epszz_input=element->GetInput(StrainRatezzEnum); _assert_(epszz_input);
			epsxz_input=element->GetInput(StrainRatexzEnum); _assert_(epsxz_input);
			epsyz_input=element->GetInput(StrainRateyzEnum); _assert_(epsyz_input);
		}

		/*Fetch number of nodes and dof for this finite element*/
		int tnumnodes = element->GetNumberOfVertices();      //Tensors, P1 DG

		/*Update tau accordingly*/
		IssmDouble* tau_xx = xNew<IssmDouble>(tnumnodes);
		IssmDouble* tau_yy = xNew<IssmDouble>(tnumnodes);
		IssmDouble* tau_xy = xNew<IssmDouble>(tnumnodes);
		IssmDouble* tau_zz = NULL;
		IssmDouble* tau_xz = NULL;
		IssmDouble* tau_yz = NULL;
		if(dim==3){
			tau_zz = xNew<IssmDouble>(tnumnodes);
			tau_xz = xNew<IssmDouble>(tnumnodes);
			tau_yz = xNew<IssmDouble>(tnumnodes);
		}
		Gauss* gauss = element->NewGauss();
		for(int ig=0;ig<tnumnodes;ig++){
			gauss->GaussNode(P1DGEnum,ig);

			/*Get D(u)*/
			vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
			vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
			if(dim==3){
				vz_input->GetInputDerivativeValue(&dvz[0],xyz_list,gauss);
			}
			epsxx = dvx[0];
			epsyy = dvy[1];
			epsxy = 0.5*(dvx[1] + dvy[0]);
			if(dim==3){
				epszz = dvz[2];
				epsxz = 0.5*(dvx[2] + dvz[0]);
				epsyz = 0.5*(dvy[2] + dvz[1]);
			}

			/*Get tau^(n-1) from inputs*/
			sigmapxx_input->GetInputValue(&sigmapxx,gauss);
			sigmapyy_input->GetInputValue(&sigmapyy,gauss);
			sigmapxy_input->GetInputValue(&sigmapxy,gauss);
			if(dim==3){
				sigmapzz_input->GetInputValue(&sigmapzz,gauss);
				sigmapxz_input->GetInputValue(&sigmapxz,gauss);
				sigmapyz_input->GetInputValue(&sigmapyz,gauss);
			}

			/*Get new d*/
			epsxx_input->GetInputValue(&d_xx,gauss);
			epsyy_input->GetInputValue(&d_yy,gauss);
			epsxy_input->GetInputValue(&d_xy,gauss);
			if(dim==3){
				epszz_input->GetInputValue(&d_zz,gauss);
				epsxz_input->GetInputValue(&d_xz,gauss);
				epsyz_input->GetInputValue(&d_yz,gauss);
			}

			/*Get d and update tau accordingly*/
			tau_xx[ig] = sigmapxx + r*(epsxx - d_xx);
			tau_yy[ig] = sigmapyy + r*(epsyy - d_yy);
			tau_xy[ig] = sigmapxy + r*(epsxy - d_xy);
			if(dim==3){
				tau_zz[ig] = sigmapzz + r*(epszz - d_zz);
				tau_xz[ig] = sigmapxz + r*(epsxz - d_xz);
				tau_yz[ig] = sigmapyz + r*(epsyz - d_yz);
			}
		}

		/*Add inputs*/
		element->AddInput(DeviatoricStressxxEnum,tau_xx,P1DGEnum);
		element->AddInput(DeviatoricStressyyEnum,tau_yy,P1DGEnum);
		element->AddInput(DeviatoricStressxyEnum,tau_xy,P1DGEnum);
		if(dim==3){
			element->AddInput(DeviatoricStresszzEnum,tau_zz,P1DGEnum);
			element->AddInput(DeviatoricStressxzEnum,tau_xz,P1DGEnum);
			element->AddInput(DeviatoricStressyzEnum,tau_yz,P1DGEnum);
		}

		/*Clean up and */
		delete gauss;
		xDelete<IssmDouble>(xyz_list);
		xDelete<IssmDouble>(tau_xx);
		xDelete<IssmDouble>(tau_yy);
		xDelete<IssmDouble>(tau_zz);
		xDelete<IssmDouble>(tau_xy);
		xDelete<IssmDouble>(tau_xz);
		xDelete<IssmDouble>(tau_yz);
	}
}/*}}}*/

/*Coupling (Tiling)*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixCouplingHOFS(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Constants*/
	int numnodes       = 3*6+1;
	int numdofp        = 2*6;
	int numdofs        = 4*6 + 3;
	int numdoftotal    = (2+4)*6+ 3;

	/*Intermediaries*/
	int   i,j,init;
	int*   cs_list     = xNew<int>(6*3+1);
	int*   cs_list2    = xNew<int>(6*2+1);
	Node  **node_list  = xNew<Node*>(6*3+1);

	/*Some parameters needed*/
	init = element->FiniteElement();

	/*prepare node list*/
	for(i=0;i<6+1;i++){
		node_list[i+6] = element->GetNode(i);
		cs_list[i+6]   = XYZEnum;
		cs_list2[i]    = XYZEnum;
	}
	for(i=0;i<6;i++){
		node_list[i]       = element->GetNode(i);
		node_list[i+2*6+1] = element->GetNode(i+6*1);
		cs_list[i]         = XYEnum;
		cs_list[i+2*6+1]   = PressureEnum;
		cs_list2[i+6+1]    = PressureEnum;
	}

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=element->NewElementMatrixCoupling(6,HOApproximationEnum);
	ElementMatrix* Ke2=element->NewElementMatrix(FSvelocityEnum);
	ElementMatrix* Ke=new ElementMatrix(Ke1,Ke2);
	delete Ke1; delete Ke2;

	/*Compute HO Matrix with P1 element type\n");*/
	Ke1=CreateKMatrixFS(element); element->TransformInvStiffnessMatrixCoord(Ke1,node_list,2*6+1,cs_list2);
	int indices[3]={18,19,20};
	Ke1->StaticCondensation(3,&indices[0]);
	element->SetTemporaryElementType(P1Enum); // P1 needed for HO
	Ke2=CreateKMatrixHO(element); element->TransformInvStiffnessMatrixCoord(Ke2,XYEnum);
	element->SetTemporaryElementType(init); // P1 needed for HO
	/*Compute FS Matrix and condense it \n");*/

	for(i=0;i<numdofs;i++) for(j=0;j<6;j++){
		Ke->values[(i+numdofp)*numdoftotal+2*j+0]+=Ke1->values[i*numdofs+3*j+0];
		Ke->values[(i+numdofp)*numdoftotal+2*j+1]+=Ke1->values[i*numdofs+3*j+1];
	}
	for(i=0;i<numdofp;i++) for(j=0;j<6;j++){
		Ke->values[i*numdoftotal+numdofp+3*j+0]+=Ke2->values[i*numdofp+2*j+0];
		Ke->values[i*numdoftotal+numdofp+3*j+1]+=Ke2->values[i*numdofp+2*j+1];
	}

	/*Transform Coordinate System*/ //Do not transform, already done in the matrices
	element->TransformStiffnessMatrixCoord(Ke,node_list,numnodes,cs_list);

	/*clean-up and return*/
	xDelete<int>(cs_list);
	xDelete<int>(cs_list2);
	xDelete<Node*>(node_list);
	delete Ke1;
	delete Ke2;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixCouplingSSAFS(Element* element){/*{{{*/

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixCouplingSSAFSViscous(element);
	ElementMatrix* Ke2=CreateKMatrixCouplingSSAFSFriction(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixCouplingSSAFSFriction(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Constants*/
	const int numdofs   = (6+1)*3 + 6*1;
	const int numdofm   = 6 *2;
	const int numdof2d  = 3 *3;
	const int numdof2dm = 3 *2;
	const int numdoftot = 6*2 + (6+1)*3 +6; // HO + FS vel + FS Pressure

	/*Intermediaries */
	int        i,j,approximation;
	int        dim=3;
	IssmDouble FSreconditioning,viscosity,alpha2_gauss,Jdet2d;
	IssmDouble bed_normal[3];
	IssmDouble LSSAFS[8][numdof2dm];
	IssmDouble LprimeSSAFS[8][numdofs];
	IssmDouble DLSSAFS[8][8]={0.0};
	IssmDouble LFSSSA[4][numdof2d];
	IssmDouble LprimeFSSSA[4][numdof2dm];
	IssmDouble DLFSSSA[4][4]={0.0};
	IssmDouble Ke_drag[numdof2dm][numdofs]={0.0};
	IssmDouble Ke_drag2[numdof2d][numdof2dm]={0.0};
	IssmDouble *xyz_list      = NULL;
	IssmDouble *xyz_list_tria = NULL;

	/*If on water or not FS, skip stiffness: */
	element->GetInputValue(&approximation,ApproximationEnum);
	if(element->IsAllFloating() || !element->IsOnBase()) return NULL;

	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numnodes  = 2*vnumnodes-1+pnumnodes;

	/*Prepare node list*/
	int* cs_list = xNew<int>(2*vnumnodes-1+pnumnodes);
	Node **node_list = xNew<Node*>(2*vnumnodes-1+pnumnodes);
	for(i=0;i<vnumnodes-1;i++){
		node_list[i] = element->GetNode(i);
		cs_list[i]   = XYEnum;
	}
	for(i=0;i<vnumnodes;i++){
		node_list[i+vnumnodes-1] = element->GetNode(i);
		cs_list[i+vnumnodes-1]   = XYZEnum;
	}
	for(i=0;i<pnumnodes;i++){
		node_list[2*vnumnodes-1+i] = element->GetNode(vnumnodes+i);
		cs_list[2*vnumnodes-1+i]   = PressureEnum;
	}

	ElementMatrix* Ke1=element->NewElementMatrixCoupling(6,SSAApproximationEnum);
	ElementMatrix* Ke2=element->NewElementMatrix(FSvelocityEnum);
	ElementMatrix* Ke=new ElementMatrix(Ke1,Ke2);
	delete Ke1; delete Ke2;

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetVerticesCoordinatesBase(&xyz_list_tria);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	Input* vx_input=element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input=element->GetInput(VzEnum); _assert_(vz_input);

	/*build friction object, used later on: */
	Friction* friction=new Friction(element,3);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(2);
	while(gauss->next()){

		element->JacobianDeterminantBase(&Jdet2d,xyz_list_tria,gauss);
		this->GetLSSAFS(&LSSAFS[0][0], element,gauss);
		this->GetLprimeSSAFS(&LprimeSSAFS[0][0], element,xyz_list, gauss);
		this->GetLFSSSA(&LFSSSA[0][0],element, gauss);
		this->GetLprimeFSSSA(&LprimeFSSSA[0][0], element,xyz_list, gauss);

		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);

		element->NormalBase(&bed_normal[0],xyz_list_tria);
		friction->GetAlpha2(&alpha2_gauss,gauss);

		IssmDouble factor = alpha2_gauss*gauss->weight*Jdet2d;
		IssmDouble factor2 = 2*viscosity*gauss->weight*Jdet2d;
		IssmDouble factor3 = FSreconditioning*gauss->weight*Jdet2d;
		DLSSAFS[0][0]=factor;
		DLSSAFS[1][1]=factor;
		DLSSAFS[2][2]=-factor*bed_normal[0]*bed_normal[2];
		DLSSAFS[3][3]=-factor*bed_normal[1]*bed_normal[2];
		DLSSAFS[4][4]=-factor2*bed_normal[0];
		DLSSAFS[5][5]=-factor2*bed_normal[1];
		DLSSAFS[6][6]=factor3*bed_normal[0];
		DLSSAFS[7][7]=factor3*bed_normal[1];

		DLFSSSA[0][0]=factor;
		DLFSSSA[1][1]=factor;
		DLFSSSA[2][2]=-factor*bed_normal[0]*bed_normal[2];
		DLFSSSA[3][3]=-factor*bed_normal[1]*bed_normal[2];

		TripleMultiply( &LSSAFS[0][0],8,numdof2dm,1,
					&DLSSAFS[0][0],8,8,0,
					&LprimeSSAFS[0][0],8,numdofs,0,
					&Ke_drag[0][0],1);

		TripleMultiply( &LFSSSA[0][0],4,numdof2d,1,
					&DLFSSSA[0][0],4,4,0,
					&LprimeFSSSA[0][0],4,numdof2dm,0,
					&Ke_drag2[0][0],1);
	}

	for(i=0;i<numdof2dm;i++) for(j=0;j<numdofs;j++) Ke->values[i*numdoftot+j+numdofm]+=Ke_drag[i][j];
	for(i=0;i<numdof2d;i++) for(j=0;j<numdof2dm;j++) Ke->values[(i+numdofm)*numdoftot+j]+=Ke_drag2[i][j];

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,node_list,numnodes,cs_list);

	/*Clean up and return*/
	xDelete<int>(cs_list);
	xDelete<Node*>(node_list);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_tria);
	delete gauss;
	delete friction;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixCouplingSSAFSViscous(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Constants*/
	const int numdofm     = 2 *3;
	const int numdofs     = 4 *6+ 3;
	const int numdoftotal = 2 *numdofm+numdofs;

	/*Intermediaries */
	int        i,j;
	int        dim=3;
	IssmDouble Jdet,viscosity,FSreconditioning,D_scalar;
	IssmDouble B[4][numdofs];
	IssmDouble Bprime[4][numdofm];
	IssmDouble B2[3][numdofm];
	IssmDouble Bprime2[3][numdofs];
	IssmDouble D[4][4]={0.0};            // material matrix, simple scalar matrix.
	IssmDouble D2[3][3]={0.0};            // material matrix, simple scalar matrix.
	IssmDouble Ke_gg[numdofs][numdofm]={0.0}; //local element stiffness matrix
	IssmDouble Ke_gg2[numdofm][numdofs]={0.0}; //local element stiffness matrix
	IssmDouble *xyz_list    = NULL;

	/*Find penta on bed as FS must be coupled to the dofs on the bed: */
	Element* pentabase=element->GetBasalElement();
	Element* basaltria=pentabase->SpawnBasalElement();

	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numnodes  = 2*vnumnodes-1+pnumnodes;

	/*Prepare node list*/
	int* cs_list     = xNew<int>(2*vnumnodes-1+pnumnodes);
	Node **node_list = xNew<Node*>(2*vnumnodes-1+pnumnodes);
	for(i=0;i<vnumnodes-1;i++){
		node_list[i] = pentabase->GetNode(i);
		cs_list[i]   = XYEnum;
	}
	for(i=0;i<vnumnodes;i++){
		node_list[i+vnumnodes-1] = element->GetNode(i);
		cs_list[i+vnumnodes-1]   = XYZEnum;
	}
	for(i=0;i<pnumnodes;i++){
		node_list[2*vnumnodes-1+i] = element->GetNode(vnumnodes+i);
		cs_list[2*vnumnodes-1+i]   = PressureEnum;
	}

	/*Initialize Element matrix and return if necessary*/
	ElementMatrix* Ke1=pentabase->NewElementMatrixCoupling(6,SSAApproximationEnum);
	ElementMatrix* Ke2=element->NewElementMatrix(FSvelocityEnum);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);
	delete Ke1; delete Ke2;

	/* Get node coordinates and dof list: */
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	Input* vx_input=element->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=element->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input=element->GetInput(VzEnum); _assert_(vz_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(5);
	Gauss* gauss_tria=new GaussTria();
	while(gauss->next()){
		gauss->SynchronizeGaussBase(gauss_tria);

		element->JacobianDeterminant(&Jdet, xyz_list,gauss);
		this->GetBSSAFS(&B[0][0],element,xyz_list, gauss);
		this->GetBprimeSSAFSTria(&Bprime[0][0], basaltria,xyz_list, gauss_tria);
		this->GetBSSAFSTria(&B2[0][0], basaltria,xyz_list, gauss_tria);
		this->GetBprimeSSAFS(&Bprime2[0][0], element,xyz_list, gauss);

		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);

		D_scalar=2*viscosity*gauss->weight*Jdet;
		for (i=0;i<3;i++) D[i][i]=D_scalar;
		D[3][3]=-gauss->weight*Jdet*FSreconditioning;
		for (i=0;i<3;i++) D2[i][i]=D_scalar;

		TripleMultiply( &B[0][0],4,numdofs,1,
					&D[0][0],4,4,0,
					&Bprime[0][0],4,numdofm,0,
					&Ke_gg[0][0],1);

		TripleMultiply( &B2[0][0],3,numdofm,1,
					&D2[0][0],3,3,0,
					&Bprime2[0][0],3,numdofs,0,
					&Ke_gg2[0][0],1);

	}
	for(i=0;i<numdofs;i++) for(j=0;j<numdofm;j++) Ke->values[(i+2*numdofm)*numdoftotal+j]+=Ke_gg[i][j];
	for(i=0;i<numdofm;i++) for(j=0;j<numdofs;j++) Ke->values[i*numdoftotal+(j+2*numdofm)]+=Ke_gg2[i][j];

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,node_list,numnodes,cs_list);

	/*Clean-up and return*/
	xDelete<int>(cs_list);
	xDelete<Node*>(node_list);
	xDelete<IssmDouble>(xyz_list);
	delete basaltria->material; delete basaltria;
	delete gauss;
	delete gauss_tria;
	return Ke;

}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixCouplingSSAHO(Element* element){/*{{{*/

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixCouplingSSAHOViscous(element);
	ElementMatrix* Ke2=CreateKMatrixCouplingSSAHOFriction(element);
	ElementMatrix* Ke=new ElementMatrix(Ke1,Ke2);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixCouplingSSAHOFriction(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	if(element->IsAllFloating() || !element->IsOnBase()) return NULL;

	/*Constants*/
	int numnodes    = element->GetNumberOfNodes();
	int numdof      = 2*numnodes;
	int numdoftotal = 4*numnodes;

	/*Intermediaries */
	int         i,j;
	IssmDouble  Jdet2d,alpha2;
	IssmDouble *xyz_list_tria = NULL;
	IssmDouble* L             = xNewZeroInit<IssmDouble>(2*numdof);
	IssmDouble  DL[2][2]      = {{ 0,0 },{0,0}}; //for basal drag
	IssmDouble  DL_scalar;
	IssmDouble* Ke_gg         = xNewZeroInit<IssmDouble>(numdof*numdof);
	Node      **node_list     = xNew<Node*>(2*numnodes);
	int*        cs_list       = xNew<int>(2*numnodes);

	/*Initialize Element matrix and return if necessary*/
	ElementMatrix* Ke1=element->NewElementMatrix(SSAApproximationEnum);
	ElementMatrix* Ke2=element->NewElementMatrix(HOApproximationEnum);
	ElementMatrix* Ke=new ElementMatrix(Ke1,Ke2);
	delete Ke1; delete Ke2;

	/*Prepare node list*/
	for(i=0;i<numnodes;i++){
		node_list[i+0*numnodes] = element->GetNode(i);
		node_list[i+1*numnodes] = element->GetNode(i);
		cs_list[i+0*numnodes] = XYEnum;
		cs_list[i+1*numnodes] = XYEnum;
	}

	/*retrieve inputs :*/
	element->GetVerticesCoordinatesBase(&xyz_list_tria);

	/*build friction object, used later on: */
	Friction* friction=new Friction(element,2);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGaussBase(2);
	while(gauss->next()){

		/*Friction: */
		friction->GetAlpha2(&alpha2,gauss);
		element->JacobianDeterminantBase(&Jdet2d, xyz_list_tria,gauss);
		this->GetBHOFriction(L,element,3,xyz_list_tria,gauss);

		DL_scalar=alpha2*gauss->weight*Jdet2d;
		for (i=0;i<2;i++) DL[i][i]=DL_scalar;

		/*  Do the triple producte tL*D*L: */
		TripleMultiply( L,2,numdof,1,
					&DL[0][0],2,2,0,
					L,2,numdof,0,
					Ke_gg,1);
	}

	for(i=0;i<numdof;i++) for(j=0;j<numdof;j++) Ke->values[i*numdoftotal+(numdof+j)]+=Ke_gg[i*numdof+j];
	for(i=0;i<numdof;i++) for(j=0;j<numdof;j++) Ke->values[(i+numdof)*numdoftotal+j]+=Ke_gg[i*numdof+j];

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,node_list,2*numnodes,cs_list);

	/*Clean up and return*/
	delete gauss;
	delete friction;
	xDelete<int>(cs_list);
	xDelete<Node*>(node_list);
	xDelete<IssmDouble>(xyz_list_tria);
	xDelete<IssmDouble>(Ke_gg);
	xDelete<IssmDouble>(L);
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixCouplingSSAHOViscous(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Constants*/
	int numnodes    = element->GetNumberOfNodes();
	int numdofm     = 1 *numnodes; //*2/2
	int numdofp     = 2 *numnodes;
	int numdoftotal = 2 *2 *numnodes;//2 dof per nodes and 2 sets of nodes for HO and SSA

	/*Intermediaries */
	int         i,j;
	IssmDouble  Jdet,viscosity;
	IssmDouble  *xyz_list      = NULL;
	IssmDouble* B              = xNew<IssmDouble>(3*numdofp);
	IssmDouble* Bprime         = xNew<IssmDouble>(3*numdofm);
	IssmDouble  D[3][3]={0.0}; // material matrix, simple scalar matrix.
	IssmDouble  D_scalar;
	IssmDouble* Ke_gg          = xNewZeroInit<IssmDouble>(numdofp*numdofm);
	Node       **node_list     = xNew<Node*>(2*numnodes);
	int*         cs_list= xNew<int>(2*numnodes);

	/*Find penta on bed as HO must be coupled to the dofs on the bed: */
	Element* pentabase=element->GetBasalElement();
	Element* basaltria=pentabase->SpawnBasalElement();

	/*prepare node list*/
	for(i=0;i<numnodes;i++){
		node_list[i+0*numnodes] = pentabase->GetNode(i);
		node_list[i+1*numnodes] = element  ->GetNode(i);
		cs_list[i+0*numnodes] = XYEnum;
		cs_list[i+1*numnodes] = XYEnum;
	}

	/*Initialize Element matrix*/
	ElementMatrix* Ke1= pentabase->NewElementMatrix(SSAApproximationEnum);
	ElementMatrix* Ke2= element  ->NewElementMatrix(HOApproximationEnum);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);
	delete Ke1; delete Ke2;

	/* Get node coordinates and dof list: */
	element->GetVerticesCoordinates(&xyz_list);
	Input* vx_input   =element->GetInput(VxEnum);       _assert_(vx_input);
	Input* vy_input   =element->GetInput(VyEnum);       _assert_(vy_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(5);
	Gauss* gauss_tria=new GaussTria();
	while(gauss->next()){
		gauss->SynchronizeGaussBase(gauss_tria);

		element->JacobianDeterminant(&Jdet, xyz_list,gauss);
		this->GetBSSAHO(B, element,xyz_list, gauss);
		this->GetBSSAprime(Bprime,basaltria,2,xyz_list, gauss_tria);
		element->material->ViscosityHO(&viscosity,3,xyz_list,gauss,vx_input,vy_input);

		D_scalar=2*viscosity*gauss->weight*Jdet;
		for (i=0;i<3;i++) D[i][i]=D_scalar;

		TripleMultiply( B,3,numdofp,1,
					&D[0][0],3,3,0,
					Bprime,3,numdofm,0,
					Ke_gg,1);
	}
	for(i=0;i<numdofp;i++) for(j=0;j<numdofm;j++) Ke->values[(i+2*numdofm)*numdoftotal+j]+=Ke_gg[i*numdofm+j];
	for(i=0;i<numdofm;i++) for(j=0;j<numdofp;j++) Ke->values[i*numdoftotal+(j+2*numdofm)]+=Ke_gg[j*numdofm+i];

	/*Transform Coordinate System*/
	element->TransformStiffnessMatrixCoord(Ke,node_list,2*numnodes,cs_list);

	/*Clean-up and return*/
	basaltria->DeleteMaterials(); delete basaltria;

	delete gauss;
	delete gauss_tria;
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(Bprime);
	xDelete<IssmDouble>(Ke_gg);
	xDelete<IssmDouble>(xyz_list);
	xDelete<Node*>(node_list);
	xDelete<int>(cs_list);
	return Ke;

}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixHOFS(Element* element){/*{{{*/

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixFS(element);
	int indices[3]={18,19,20};
	Ke1->StaticCondensation(3,&indices[0]);
	int init = element->FiniteElement();
	element->SetTemporaryElementType(P1Enum); // P1 needed for HO
	ElementMatrix* Ke2=CreateKMatrixHO(element);
	element->SetTemporaryElementType(init); // P1 needed for HO
	ElementMatrix* Ke3=CreateKMatrixCouplingHOFS(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2,Ke3);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	delete Ke3;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixSSAFS(Element* element){/*{{{*/

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixFS(element);
	int indices[3]={18,19,20};
	Ke1->StaticCondensation(3,&indices[0]);
	int init = element->FiniteElement();
	element->SetTemporaryElementType(P1Enum);
	ElementMatrix* Ke2=CreateKMatrixSSA3d(element);
	element->SetTemporaryElementType(init);
	ElementMatrix* Ke3=CreateKMatrixCouplingSSAFS(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2,Ke3);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	delete Ke3;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixSSAHO(Element* element){/*{{{*/

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixSSA3d(element);
	ElementMatrix* Ke2=CreateKMatrixHO(element);
	ElementMatrix* Ke3=CreateKMatrixCouplingSSAHO(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2,Ke3);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	delete Ke3;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixSSA3d(Element* element){/*{{{*/

	/*compute all stiffness matrices for this element*/
	ElementMatrix* Ke1=CreateKMatrixSSA3dViscous(element);
	ElementMatrix* Ke2=CreateKMatrixSSA3dFriction(element);
	ElementMatrix* Ke =new ElementMatrix(Ke1,Ke2);

	/*clean-up and return*/
	delete Ke1;
	delete Ke2;
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixSSA3dFriction(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Initialize Element matrix and return if necessary*/
	if(element->IsAllFloating() || !element->IsOnBase()) return NULL;

	/*Build a tria element using the 3 nodes of the base of the penta. Then use
	 * the tria functionality to build a friction stiffness matrix on these 3
	 * nodes: */
	Element* basalelement = element->SpawnBasalElement();
	ElementMatrix* Ke=CreateKMatrixSSAFriction(basalelement);
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};

	/*clean-up and return*/
	return Ke;
}/*}}}*/
ElementMatrix* StressbalanceAnalysis::CreateKMatrixSSA3dViscous(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Constants*/
	const int    numdof2d=2*3;

	/*Intermediaries */
	int         i,j,approximation;
	int         dim=3;
	IssmDouble  Jdet,viscosity;
	IssmDouble  epsilon[5],oldepsilon[5];       /* epsilon=[exx,eyy,exy,exz,eyz];*/
	IssmDouble  epsilons[6];                    //6 for FS
	IssmDouble  B[3][numdof2d];
	IssmDouble  Bprime[3][numdof2d];
	IssmDouble  D[3][3]= {0.0};                 // material matrix, simple scalar matrix.
	IssmDouble  D_scalar;
	IssmDouble  Ke_gg[numdof2d][numdof2d]={0.0};
	IssmDouble  *xyz_list  = NULL;

	/*Find penta on bed as this is a SSA elements: */
	Element* pentabase=element->GetBasalElement();
	Element* basaltria=pentabase->SpawnBasalElement(true);

	/*Initialize Element matrix*/
	ElementMatrix* Ke=basaltria->NewElementMatrix(SSAApproximationEnum);
	element->GetInputValue(&approximation,ApproximationEnum);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vx_input   =element->GetInput(VxEnum);       _assert_(vx_input);
	Input* vy_input   =element->GetInput(VyEnum);       _assert_(vy_input);
	Input* vz_input   =element->GetInput(VzEnum);       _assert_(vz_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(5);
	Gauss* gauss_tria=new GaussTria();
	while(gauss->next()){
		gauss->SynchronizeGaussBase(gauss_tria);

		element->JacobianDeterminant(&Jdet, xyz_list,gauss);
		this->GetBSSA(&B[0][0],basaltria,2,xyz_list, gauss_tria);
		this->GetBSSAprime(&Bprime[0][0],basaltria,2,xyz_list, gauss_tria);

		if(approximation==SSAHOApproximationEnum){
			element->material->ViscosityHO(&viscosity,dim,xyz_list,gauss,vx_input,vy_input);
		}
		else if (approximation==SSAFSApproximationEnum){
			element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);
		}
		else _error_("approximation " << approximation << " (" << EnumToStringx(approximation) << ") not supported yet");

		D_scalar=2*viscosity*gauss->weight*Jdet;
		for (i=0;i<3;i++) D[i][i]=D_scalar;

		TripleMultiply( &B[0][0],3,numdof2d,1,
					&D[0][0],3,3,0,
					&Bprime[0][0],3,numdof2d,0,
					&Ke_gg[0][0],1);

	}
	for(i=0;i<numdof2d;i++) for(j=0;j<numdof2d;j++) Ke->values[i*numdof2d+j]+=Ke_gg[i][j];

	/*Transform Coordinate System*/
	basaltria->TransformStiffnessMatrixCoord(Ke,XYEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	delete basaltria->material;
	delete basaltria;
	delete gauss_tria;
	delete gauss;
	return Ke;

}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorCouplingHOFS(Element* element){/*{{{*/

	/*compute all load vectors for this element*/
	ElementVector* pe1=CreatePVectorCouplingHOFSViscous(element);
	ElementVector* pe2=CreatePVectorCouplingHOFSFriction(element);
	ElementVector* pe =new ElementVector(pe1,pe2);

	/*clean-up and return*/
	delete pe1;
	delete pe2;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorCouplingHOFSFriction(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int         i,approximation;
	int         dim=3;
	IssmDouble  Jdet,Jdet2d,FSreconditioning;
	IssmDouble	bed_normal[3];
	IssmDouble  viscosity, w, alpha2_gauss;
	IssmDouble  dw[3];
	IssmDouble	*xyz_list_tria = NULL;
	IssmDouble  *xyz_list      = NULL;
	IssmDouble  basis[6]; //for the six nodes of the penta

	/*Initialize Element vector and return if necessary*/
	if(!element->IsOnBase() || element->IsAllFloating()) return NULL;
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation!=HOFSApproximationEnum) return NULL;

	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();
	int numnodes  = vnumnodes+pnumnodes;

	/*Prepare coordinate system list*/
	int*   cs_list   = xNew<int>(vnumnodes+pnumnodes);
	Node **node_list = xNew<Node*>(vnumnodes+pnumnodes);
	for(i=0;i<vnumnodes;i++){
		cs_list[i]           = XYZEnum;
		node_list[i]           = element->GetNode(i);
	}
	for(i=0;i<pnumnodes;i++){
		cs_list[vnumnodes+i] = PressureEnum;
		node_list[vnumnodes+i] = element->GetNode(vnumnodes+i);
	}

	ElementVector* pe=element->NewElementVector(FSvelocityEnum);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetVerticesCoordinatesBase(&xyz_list_tria);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	Input* vx_input=  element->GetInput(VxEnum);   _assert_(vx_input);
	Input* vy_input=  element->GetInput(VyEnum);   _assert_(vy_input);
	Input* vz_input=  element->GetInput(VzEnum);   _assert_(vz_input);
	Input* vzHO_input=element->GetInput(VzHOEnum); _assert_(vzHO_input);

	/*build friction object, used later on: */
	Friction* friction=new Friction(element,3);

	/* Start looping on the number of gauss 2d (nodes on the bedrock) */
	Gauss* gauss=element->NewGaussBase(2);
	while(gauss->next()){

		element->JacobianDeterminantBase(&Jdet2d,xyz_list_tria,gauss);
		element->NodalFunctionsP1(basis, gauss);

		vzHO_input->GetInputValue(&w, gauss);
		vzHO_input->GetInputDerivativeValue(&dw[0],xyz_list,gauss);

		element->NormalBase(&bed_normal[0],xyz_list_tria);
		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);
		friction->GetAlpha2(&alpha2_gauss,gauss);

		IssmDouble factor = Jdet2d*gauss->weight;
		for(i=0;i<3;i++){
			pe->values[i*3+0]+=factor*(alpha2_gauss*w*bed_normal[0]*bed_normal[2]+2*viscosity*dw[2]*bed_normal[0])*basis[i];
			pe->values[i*3+1]+=factor*(alpha2_gauss*w*bed_normal[1]*bed_normal[2]+2*viscosity*dw[2]*bed_normal[1])*basis[i];
			pe->values[i*3+2]+=factor*2*viscosity*(dw[0]*bed_normal[0]+dw[1]*bed_normal[1]+dw[2]*bed_normal[2])*basis[i];
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,node_list,vnumnodes+pnumnodes,cs_list);

	/*Clean up and return*/
	xDelete<int>(cs_list);
	xDelete<Node*>(node_list);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_tria);
	delete gauss;
	delete friction;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorCouplingHOFSViscous(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         i,approximation;
	int         dim=3;
	IssmDouble  viscosity,Jdet,FSreconditioning;
	IssmDouble  dw[3];
	IssmDouble  *xyz_list = NULL;
	IssmDouble  basis[6]; //for the six nodes of the penta
	IssmDouble  dbasis[3][6]; //for the six nodes of the penta

	/*Initialize Element vector and return if necessary*/
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation!=HOFSApproximationEnum) return NULL;
	int   vnumnodes = element->NumberofNodesVelocity();
	int   pnumnodes = element->NumberofNodesPressure();

	/*Prepare coordinate system list*/
	int*   cs_list   = xNew<int>(vnumnodes+pnumnodes);
	Node **node_list = xNew<Node*>(vnumnodes+pnumnodes);
	for(i=0;i<vnumnodes;i++){
		cs_list[i]             = XYZEnum;
		node_list[i]           = element->GetNode(i);
	}
	for(i=0;i<pnumnodes;i++){
		cs_list[vnumnodes+i]   = PressureEnum;
		node_list[vnumnodes+i] = element->GetNode(vnumnodes+i);
	}
	ElementVector* pe = element->NewElementVector(FSvelocityEnum);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	Input* vx_input   =element->GetInput(VxEnum);   _assert_(vx_input);
	Input* vy_input   =element->GetInput(VyEnum);   _assert_(vy_input);
	Input* vz_input   =element->GetInput(VzEnum);   _assert_(vz_input);
	Input* vzHO_input=element->GetInput(VzHOEnum);  _assert_(vzHO_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(5);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet, xyz_list,gauss);
		element->NodalFunctionsP1(&basis[0],gauss);
		element->NodalFunctionsP1Derivatives(&dbasis[0][0],xyz_list,gauss);

		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);
		vzHO_input->GetInputDerivativeValue(&dw[0],xyz_list,gauss);

		IssmDouble factor = -Jdet*gauss->weight*viscosity;
		IssmDouble factor2 = Jdet*gauss->weight*FSreconditioning;
		for(i=0;i<6;i++){
			pe->values[i*3+0]+=factor*dw[0]*dbasis[2][i];
			pe->values[i*3+1]+=factor*dw[1]*dbasis[2][i];
			pe->values[i*3+2]+=factor*(dw[0]*dbasis[0][i]+dw[1]*dbasis[1][i]+2*dw[2]*dbasis[2][i]);
			pe->values[3*vnumnodes+i]+=factor2*dw[2]*basis[i];
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,node_list,vnumnodes+pnumnodes,cs_list);

	/*Clean up and return*/
	xDelete<int>(cs_list);
	xDelete<Node*>(node_list);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorCouplingSSAFS(Element* element){/*{{{*/

	/*compute all load vectors for this element*/
	ElementVector* pe1=CreatePVectorCouplingSSAFSViscous(element);
	ElementVector* pe2=CreatePVectorCouplingSSAFSFriction(element);
	ElementVector* pe =new ElementVector(pe1,pe2);

	/*clean-up and return*/
	delete pe1;
	delete pe2;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorCouplingSSAFSFriction(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int         i,j,approximation;
	int         dim=3;
	IssmDouble  Jdet,Jdet2d,FSreconditioning;
	IssmDouble	bed_normal[3];
	IssmDouble  viscosity, w, alpha2_gauss;
	IssmDouble  dw[3];
	IssmDouble  basis[6]; //for the six nodes of the penta
	IssmDouble	*xyz_list_tria = NULL;
	IssmDouble  *xyz_list      = NULL;

	/*Initialize Element vector and return if necessary*/
	if(!element->IsOnBase() || element->IsAllFloating()) return NULL;
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation!=SSAFSApproximationEnum) return NULL;
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();

	/*Prepare coordinate system list*/
	int* cs_list     = xNew<int>(vnumnodes+pnumnodes);
	Node **node_list = xNew<Node*>(vnumnodes+pnumnodes);
	for(i=0;i<vnumnodes;i++){
		cs_list[i]             = XYZEnum;
		node_list[i]           = element->GetNode(i);
	}
	for(i=0;i<pnumnodes;i++){
		cs_list[vnumnodes+i]   = PressureEnum;
		node_list[vnumnodes+i] = element->GetNode(vnumnodes+i);
	}
	ElementVector* pe=element->NewElementVector(FSvelocityEnum);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetVerticesCoordinatesBase(&xyz_list_tria);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	Input* vx_input=   element->GetInput(VxEnum);    _assert_(vx_input);
	Input* vy_input=   element->GetInput(VyEnum);    _assert_(vy_input);
	Input* vz_input=   element->GetInput(VzEnum);    _assert_(vz_input);
	Input* vzSSA_input=element->GetInput(VzSSAEnum); _assert_(vzSSA_input);

	/*build friction object, used later on: */
	Friction* friction=new Friction(element,3);

	/* Start looping on the number of gauss 2d (nodes on the bedrock) */
	Gauss* gauss=element->NewGaussBase(2);
	while(gauss->next()){

		element->JacobianDeterminantBase(&Jdet2d,xyz_list_tria,gauss);
		element->NodalFunctionsP1(basis, gauss);

		vzSSA_input->GetInputValue(&w, gauss);
		vzSSA_input->GetInputDerivativeValue(&dw[0],xyz_list,gauss);

		element->NormalBase(&bed_normal[0],xyz_list_tria);
		element->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,vz_input);
		friction->GetAlpha2(&alpha2_gauss,gauss);

		IssmDouble factor = Jdet2d*gauss->weight;
		for(i=0;i<3;i++){
			pe->values[i*3+0]+=factor*(alpha2_gauss*w*bed_normal[0]*bed_normal[2]+2*viscosity*dw[2]*bed_normal[0])*basis[i];
			pe->values[i*3+1]+=factor*(alpha2_gauss*w*bed_normal[1]*bed_normal[2]+2*viscosity*dw[2]*bed_normal[1])*basis[i];
			pe->values[i*3+2]+=factor*2*viscosity*(dw[0]*bed_normal[0]+dw[1]*bed_normal[1]+dw[2]*bed_normal[2])*basis[i];
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,node_list,vnumnodes+pnumnodes,cs_list);

	/*Clean up and return*/
	xDelete<int>(cs_list);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_tria);
	xDelete<Node*>(node_list);
	delete gauss;
	delete friction;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorCouplingSSAFSViscous(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int         i,approximation;
	IssmDouble  viscosity,Jdet,FSreconditioning;
	IssmDouble  dw[3];
	IssmDouble  *xyz_list = NULL;
	IssmDouble  basis[6]; //for the six nodes of the penta
	IssmDouble  dbasis[3][6]; //for the six nodes of the penta

	/*Initialize Element vector and return if necessary*/
	element->GetInputValue(&approximation,ApproximationEnum);
	if(approximation!=SSAFSApproximationEnum) return NULL;
	int vnumnodes = element->NumberofNodesVelocity();
	int pnumnodes = element->NumberofNodesPressure();

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(vnumnodes+pnumnodes);
	Node **node_list = xNew<Node*>(vnumnodes+pnumnodes);
	for(i=0;i<vnumnodes;i++){
		cs_list[i]             = XYZEnum;
		node_list[i]           = element->GetNode(i);
	}
	for(i=0;i<pnumnodes;i++){
		cs_list[vnumnodes+i]   = PressureEnum;
		node_list[vnumnodes+i] = element->GetNode(vnumnodes+i);
	}
	ElementVector* pe=element->NewElementVector(FSvelocityEnum);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	Input* vx_input   =element->GetInput(VxEnum);      _assert_(vx_input);
	Input* vy_input   =element->GetInput(VyEnum);      _assert_(vy_input);
	Input* vz_input   =element->GetInput(VzEnum);      _assert_(vz_input);
	Input* vzSSA_input=element->GetInput(VzSSAEnum);   _assert_(vzSSA_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(5);
	while(gauss->next()){
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctionsP1(&basis[0], gauss);
		element->NodalFunctionsP1Derivatives(&dbasis[0][0],xyz_list, gauss);

		vzSSA_input->GetInputDerivativeValue(&dw[0],xyz_list,gauss);
		element->material->ViscosityFS(&viscosity,3,xyz_list,gauss,vx_input,vy_input,vz_input);

		IssmDouble factor = -Jdet*gauss->weight*viscosity;
		IssmDouble factor2 = Jdet*gauss->weight*FSreconditioning;
		for(i=0;i<6;i++){
			pe->values[i*3+0]+=factor*dw[0]*dbasis[2][i];
			pe->values[i*3+1]+=factor*dw[1]*dbasis[2][i];
			pe->values[i*3+2]+=factor*(dw[0]*dbasis[0][i]+dw[1]*dbasis[1][i]+2*dw[2]*dbasis[2][i]);
			pe->values[3*vnumnodes+i]+=factor2*dw[2]*basis[i];
		}
	}

	/*Transform coordinate system*/
	element->TransformLoadVectorCoord(pe,node_list,vnumnodes+pnumnodes,cs_list);

	/*Clean up and return*/
	xDelete<int>(cs_list);
	xDelete<Node*>(node_list);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorHOFS(Element* element){/*{{{*/

	/*compute all load vectors for this element*/
	int init = element->FiniteElement();
	element->SetTemporaryElementType(P1Enum);
	ElementVector* pe1=CreatePVectorHO(element);
	element->SetTemporaryElementType(init);
	ElementVector* pe2=CreatePVectorFS(element);
	int indices[3]={18,19,20};
	element->SetTemporaryElementType(MINIcondensedEnum);
	ElementMatrix* Ke = CreateKMatrixFS(element);
	element->SetTemporaryElementType(init);
	pe2->StaticCondensation(Ke,3,&indices[0]);
	delete Ke;
	ElementVector* pe3=CreatePVectorCouplingHOFS(element);
	ElementVector* pe =new ElementVector(pe1,pe2,pe3);

	/*clean-up and return*/
	delete pe1;
	delete pe2;
	delete pe3;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorSSAFS(Element* element){/*{{{*/

	/*compute all load vectors for this element*/
	int init = element->FiniteElement();
	element->SetTemporaryElementType(P1Enum); // P1 needed for HO
	ElementVector* pe1=CreatePVectorSSA(element);
	element->SetTemporaryElementType(init); // P1 needed for HO
	ElementVector* pe2=CreatePVectorFS(element);
	int indices[3]={18,19,20};
	element->SetTemporaryElementType(MINIcondensedEnum); // P1 needed for HO
	ElementMatrix* Ke = CreateKMatrixFS(element);
	element->SetTemporaryElementType(init); // P1 needed for HO
	pe2->StaticCondensation(Ke,3,&indices[0]);
	delete Ke;
	ElementVector* pe3=CreatePVectorCouplingSSAFS(element);
	ElementVector* pe =new ElementVector(pe1,pe2,pe3);

	/*clean-up and return*/
	delete pe1;
	delete pe2;
	delete pe3;
	return pe;
}/*}}}*/
ElementVector* StressbalanceAnalysis::CreatePVectorSSAHO(Element* element){/*{{{*/

	/*compute all load vectors for this element*/
	ElementVector* pe1=CreatePVectorSSA(element);
	ElementVector* pe2=CreatePVectorHO(element);
	ElementVector* pe =new ElementVector(pe1,pe2);

	/*clean-up and return*/
	delete pe1;
	delete pe2;
	return pe;
}/*}}}*/
void           StressbalanceAnalysis::GetBprimeSSAFS(IssmDouble* Bprime,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute Bprime  matrix. Bprime=[Bprime1 Bprime2 Bprime3 Bprime4 Bprime5 Bprime6] where Bprimei is of size 5*2.
	 * For node i, Bprimei can be expressed in the actual coordinate system
	 * by:
	 *       Bprimei=[ 2*dh/dx    dh/dy   0   0 ]
	 *               [  dh/dx    2*dh/dy  0   0 ]
	 *               [  dh/dy     dh/dx   0   0 ]
	 * where h is the interpolation function for node i.
	 *
	 * We assume Bprime has been allocated already, of size: 5x(2*NUMNODESP1)
	 */

	int    i;
	IssmDouble dbasismini[3][7];

	/*Get dbasis in actual coordinate system: */
	element->NodalFunctionsMINIDerivatives(&dbasismini[0][0],xyz_list, gauss);

	/*Build Bprime: */
	for(i=0;i<6;i++){
		Bprime[(3*7+6)*0+3*i+0] = 2.*dbasismini[0][i];
		Bprime[(3*7+6)*0+3*i+1] = dbasismini[1][i];
		Bprime[(3*7+6)*0+3*i+2] = 0.;
		Bprime[(3*7+6)*1+3*i+0] = dbasismini[0][i];
		Bprime[(3*7+6)*1+3*i+1] = 2.*dbasismini[1][i];
		Bprime[(3*7+6)*1+3*i+2] = 0.;
		Bprime[(3*7+6)*2+3*i+0] = dbasismini[1][i];
		Bprime[(3*7+6)*2+3*i+1] = dbasismini[0][i];
		Bprime[(3*7+6)*2+3*i+2] = 0.;
	}

	for(i=0;i<1;i++){ //Add zeros for the bubble function
		Bprime[(3*7+6)*0+3*(6+i)+0] = 0.;
		Bprime[(3*7+6)*0+3*(6+i)+1] = 0.;
		Bprime[(3*7+6)*0+3*(6+i)+2] = 0.;
		Bprime[(3*7+6)*1+3*(6+i)+0] = 0.;
		Bprime[(3*7+6)*1+3*(6+i)+1] = 0.;
		Bprime[(3*7+6)*1+3*(6+i)+2] = 0.;
		Bprime[(3*7+6)*2+3*(6+i)+0] = 0.;
		Bprime[(3*7+6)*2+3*(6+i)+1] = 0.;
		Bprime[(3*7+6)*2+3*(6+i)+2] = 0.;
	}

	for(i=0;i<6;i++){ //last column not for the bubble function
		Bprime[(3*7+6)*0+7*3+i] = 0.;
		Bprime[(3*7+6)*1+7*3+i] = 0.;
		Bprime[(3*7+6)*2+7*3+i] = 0.;
	}
}/*}}}*/
void           StressbalanceAnalysis::GetBprimeSSAFSTria(IssmDouble* Bprime,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute Bprime  matrix. Bprime=[Bprime1 Bprime2 Bprime3] where Bprimei is of size 3*2.
	 * For node i, Bprimei can be expressed in the actual coordinate system
	 * by:
	 *       Bprimei=[  dN/dx    0   ]
	 *               [    0    dN/dy ]
	 *               [  dN/dy  dN/dx ]
	 N               [  dN/dx  dN/dy ]
	 * where N is the finiteelement function for node i.
	 *
	 * We assume Bprime has been allocated already, of size: 3x(2*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions*/
	IssmDouble* dbasis=xNew<IssmDouble>(2*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build Bprime: */
	for(int i=0;i<numnodes;i++){
		Bprime[2*numnodes*0+2*i+0] = dbasis[0*numnodes+i];
		Bprime[2*numnodes*0+2*i+1] = 0.;
		Bprime[2*numnodes*1+2*i+0] = 0.;
		Bprime[2*numnodes*1+2*i+1] = dbasis[1*numnodes+i];
		Bprime[2*numnodes*2+2*i+0] = dbasis[1*numnodes+i];
		Bprime[2*numnodes*2+2*i+1] = dbasis[0*numnodes+i];
		Bprime[2*numnodes*3+2*i+0] = dbasis[0*numnodes+i];
		Bprime[2*numnodes*3+2*i+1] = dbasis[1*numnodes+i];
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetBSSAFS(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3 B4 B5 B6] where Bi is of size 5*2.
	 * For node i, Bi can be expressed in the actual coordinate system
	 * by:
	 *       Bi=[ dh/dx          0       0   0 ]
	 *          [   0           dh/dy    0   0 ]
	 *          [ 1/2*dh/dy  1/2*dh/dx   0   0 ]
	 *          [   0            0       0   h ]
	 * where h is the interpolation function for node i.
	 *
	 * We assume B has been allocated already, of size: 5x(2*NUMNODESP1)
	 */

	int i;
	IssmDouble dbasismini[3][7];
	IssmDouble basis[6];

	/*Get dbasis in actual coordinate system: */
	element->NodalFunctionsMINIDerivatives(&dbasismini[0][0],xyz_list, gauss);
	element->NodalFunctionsP1(basis,gauss);

	/*Build B: */
	for(i=0;i<6;i++){
		B[(3*7+6)*0+3*i+0] = dbasismini[0][i];
		B[(3*7+6)*0+3*i+1] = 0.;
		B[(3*7+6)*0+3*i+2] = 0.;
		B[(3*7+6)*1+3*i+0] = 0.;
		B[(3*7+6)*1+3*i+1] = dbasismini[1][i];
		B[(3*7+6)*1+3*i+2] = 0.;
		B[(3*7+6)*2+3*i+0] = 0.5*dbasismini[1][i];
		B[(3*7+6)*2+3*i+1] = 0.5*dbasismini[0][i];
		B[(3*7+6)*2+3*i+2] = 0.;
		B[(3*7+6)*3+3*i+0] = 0.;
		B[(3*7+6)*3+3*i+1] = 0.;
		B[(3*7+6)*3+3*i+2] = 0.;
	}
	for(i=0;i<1;i++){
		B[(3*7+6)*0+3*(6+i)+0] = 0.;
		B[(3*7+6)*0+3*(6+i)+1] = 0.;
		B[(3*7+6)*0+3*(6+i)+2] = 0.;
		B[(3*7+6)*1+3*(6+i)+0] = 0.;
		B[(3*7+6)*1+3*(6+i)+1] = 0.;
		B[(3*7+6)*1+3*(6+i)+2] = 0.;
		B[(3*7+6)*2+3*(6+i)+0] = 0.;
		B[(3*7+6)*2+3*(6+i)+1] = 0.;
		B[(3*7+6)*2+3*(6+i)+2] = 0.;
		B[(3*7+6)*3+3*(6+i)+0] = 0.;
		B[(3*7+6)*3+3*(6+i)+1] = 0.;
		B[(3*7+6)*3+3*(6+i)+2] = 0.;
	}

	for(i=0;i<6;i++){ //last column not for the bubble function
		B[(3*7+6)*0+7*3+i] = 0;
		B[(3*7+6)*1+7*3+i] = 0;
		B[(3*7+6)*2+7*3+i] = 0;
		B[(3*7+6)*3+7*3+i] = basis[i];
	}
}/*}}}*/
void           StressbalanceAnalysis::GetBSSAFSTria(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3] where Bi is of size 3*2.
	 * For node i, Bi can be expressed in the actual coordinate system
	 * by:
	 *       Bi=[   dN/dx         0     ]
	 *          [       0       dN/dy   ]
	 *          [  1/2*dN/dy  1/2*dN/dx ]
	 * where N is the finiteelement function for node i.
	 *
	 * We assume B has been allocated already, of size: 3x(2*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* dbasis=xNew<IssmDouble>(2*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B': */
	for(int i=0;i<numnodes;i++){
		B[2*numnodes*0+2*i+0] = dbasis[0*numnodes+i];
		B[2*numnodes*0+2*i+1] = 0.;
		B[2*numnodes*1+2*i+0] = 0.;
		B[2*numnodes*1+2*i+1] = dbasis[1*numnodes+i];
		B[2*numnodes*2+2*i+0] = 0.5*dbasis[1*numnodes+i];
		B[2*numnodes*2+2*i+1] = 0.5*dbasis[0*numnodes+i];
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetBSSAHO(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3 B4 B5 B6] where Bi is of size 3*2.
	 * For node i, Bi can be expressed in the actual coordinate system
	 * by:
	 *       Bi=[ dh/dx          0      ]
	 *          [   0           dh/dy   ]
	 *          [ 1/2*dh/dy  1/2*dh/dx  ]
	 * where h is the interpolation function for node i.
	 *
	 * We assume B has been allocated already, of size: 5x(2*NUMNODESP1)
	 */

	int numnodes = element->GetNumberOfNodes();
	IssmDouble* dbasis=xNew<IssmDouble>(3*numnodes);

	/*Get dbasis in actual coordinate system: */
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B: */
	for(int i=0;i<numnodes;i++){
		B[2*numnodes*0+2*i+0] = dbasis[0*numnodes+i];
		B[2*numnodes*0+2*i+1] = 0.;
		B[2*numnodes*1+2*i+0] = 0.;
		B[2*numnodes*1+2*i+1] = dbasis[1*numnodes+i];
		B[2*numnodes*2+2*i+0] = .5*dbasis[1*numnodes+i];
		B[2*numnodes*2+2*i+1] = .5*dbasis[0*numnodes+i];
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);
}/*}}}*/
void           StressbalanceAnalysis::GetLprimeFSSSA(IssmDouble* LprimeFS,Element* element,IssmDouble* xyz_list,Gauss* gauss_in){/*{{{*/
	/* Compute Lprime  matrix. Lprime=[Lp1 Lp2 Lp3] where Lpi is square and of size numdof.
	 * For node i, Lpi can be expressed in the actual coordinate system
	 * by:
	 *       Lpi=[ h    0 ]
	 *		       [ 0    h ]
	 *		       [ h    0 ]
	 *		       [ 0    h ]
	 * where h is the interpolation function for node i.
	 */
	int num_dof=2;
	IssmDouble basis[3];

	/*Cast gauss to GaussPenta*/
	_assert_(gauss_in->Enum()==GaussPentaEnum);
	GaussPenta* gauss = xDynamicCast<GaussPenta*>(gauss_in);

	/*Get basis in actual coordinate system: */
	basis[0]=gauss->coord1*(1-gauss->coord4)/2.0;
	basis[1]=gauss->coord2*(1-gauss->coord4)/2.0;
	basis[2]=gauss->coord3*(1-gauss->coord4)/2.0;

	/*Build LprimeFS: */
	for(int i=0;i<3;i++){
		LprimeFS[num_dof*3*0+num_dof*i+0] = basis[i];
		LprimeFS[num_dof*3*0+num_dof*i+1] = 0.;
		LprimeFS[num_dof*3*1+num_dof*i+0] = 0.;
		LprimeFS[num_dof*3*1+num_dof*i+1] = basis[i];
		LprimeFS[num_dof*3*2+num_dof*i+0] = basis[i];
		LprimeFS[num_dof*3*2+num_dof*i+1] = 0.;
		LprimeFS[num_dof*3*3+num_dof*i+0] = 0.;
		LprimeFS[num_dof*3*3+num_dof*i+1] = basis[i];
	}
}/*}}}*/
void           StressbalanceAnalysis::GetLprimeSSAFS(IssmDouble* LprimeFS,Element* element,IssmDouble* xyz_list,Gauss* gauss_in){/*{{{*/
	/* Compute Lprime  matrix. Lprime=[Lp1 Lp2 Lp3] where Lpi is square and of size numdof.
	 * For node i, Lpi can be expressed in the actual coordinate system
	 * by:
	 *       Lpi=[ h    0    0   0]
	 *		       [ 0    h    0   0]
	 *		       [ 0    0    h   0]
	 *		       [ 0    0    h   0]
	 *		       [ 0    0  dh/dz 0]
	 *		       [ 0    0  dh/dz 0]
	 *           [ 0    0    0   h]
	 *           [ 0    0    0   h]
	 * where h is the interpolation function for node i.
	 */
	int num_dof=3;
	int num_dof_vel=3*7;
	int num_dof_total=3*7+1*6;
	IssmDouble basis[3];
	IssmDouble dbasis[3][6];

	/*Cast gauss to GaussPenta*/
	_assert_(gauss_in->Enum()==GaussPentaEnum);
	GaussPenta* gauss = xDynamicCast<GaussPenta*>(gauss_in);

	/*Get basis in actual coordinate system: */
	basis[0]=gauss->coord1*(1-gauss->coord4)/2.0;
	basis[1]=gauss->coord2*(1-gauss->coord4)/2.0;
	basis[2]=gauss->coord3*(1-gauss->coord4)/2.0;

	element->NodalFunctionsP1Derivatives(&dbasis[0][0],xyz_list,gauss);

	/*Build LprimeFS: */
	for(int i=0;i<3;i++){
		LprimeFS[num_dof_total*0+num_dof*i+0] = basis[i];
		LprimeFS[num_dof_total*0+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*0+num_dof*i+2] = 0.;
		LprimeFS[num_dof_total*1+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*1+num_dof*i+1] = basis[i];
		LprimeFS[num_dof_total*1+num_dof*i+2] = 0.;
		LprimeFS[num_dof_total*2+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*2+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*2+num_dof*i+2] = basis[i];
		LprimeFS[num_dof_total*3+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*3+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*3+num_dof*i+2] = basis[i];
		LprimeFS[num_dof_total*4+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*4+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*4+num_dof*i+2] = dbasis[2][i];
		LprimeFS[num_dof_total*5+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*5+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*5+num_dof*i+2] = dbasis[2][i];
		LprimeFS[num_dof_total*6+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*6+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*6+num_dof*i+2] = 0.;
		LprimeFS[num_dof_total*7+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*7+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*7+num_dof*i+2] = 0.;
	}
	for(int i=3;i<7;i++){
		LprimeFS[num_dof_total*0+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*0+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*0+num_dof*i+2] = 0.;
		LprimeFS[num_dof_total*1+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*1+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*1+num_dof*i+2] = 0.;
		LprimeFS[num_dof_total*2+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*2+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*2+num_dof*i+2] = 0.;
		LprimeFS[num_dof_total*3+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*3+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*3+num_dof*i+2] = 0.;
		LprimeFS[num_dof_total*4+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*4+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*4+num_dof*i+2] = 0.;
		LprimeFS[num_dof_total*5+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*5+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*5+num_dof*i+2] = 0.;
		LprimeFS[num_dof_total*6+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*6+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*6+num_dof*i+2] = 0.;
		LprimeFS[num_dof_total*7+num_dof*i+0] = 0.;
		LprimeFS[num_dof_total*7+num_dof*i+1] = 0.;
		LprimeFS[num_dof_total*7+num_dof*i+2] = 0.;
	}
	for(int i=0;i<3;i++){
		LprimeFS[num_dof_total*0+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*1+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*2+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*3+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*4+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*5+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*6+num_dof_vel+i] = basis[i];
		LprimeFS[num_dof_total*7+num_dof_vel+i] = basis[i];
	}
	for(int i=3;i<6;i++){
		LprimeFS[num_dof_total*0+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*1+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*2+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*3+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*4+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*5+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*6+num_dof_vel+i] = 0.;
		LprimeFS[num_dof_total*7+num_dof_vel+i] = 0.;
	}
}/*}}}*/
void           StressbalanceAnalysis::GetLFSSSA(IssmDouble* LFS,Element* element,Gauss* gauss_in){/*{{{*/
	/* Compute L  matrix. L=[L1 L2 L3] where Li is square and of size numdof.
	 * For node i, Li can be expressed in the actual coordinate system
	 * by:
	 *       Li=[ h    0    0 ]
	 *	 	      [ 0    h    0 ]
	 *		      [ 0    0    h ]
	 *		      [ 0    0    h ]
	 * where h is the interpolation function for node i.
	 */

	int num_dof=3;
	IssmDouble basis[3];

	/*Cast gauss to GaussPenta*/
	_assert_(gauss_in->Enum()==GaussPentaEnum);
	GaussPenta* gauss = xDynamicCast<GaussPenta*>(gauss_in);

	/*Get basis in actual coordinate system: */
	basis[0]=gauss->coord1*(1-gauss->coord4)/2.0;
	basis[1]=gauss->coord2*(1-gauss->coord4)/2.0;
	basis[2]=gauss->coord3*(1-gauss->coord4)/2.0;

	/*Build LFS: */
	for(int i=0;i<3;i++){
		LFS[num_dof*3*0+num_dof*i+0] = basis[i];
		LFS[num_dof*3*0+num_dof*i+1] = 0.;
		LFS[num_dof*3*0+num_dof*i+2] = 0.;
		LFS[num_dof*3*1+num_dof*i+0] = 0.;
		LFS[num_dof*3*1+num_dof*i+1] = basis[i];
		LFS[num_dof*3*1+num_dof*i+2] = 0.;
		LFS[num_dof*3*2+num_dof*i+0] = 0.;
		LFS[num_dof*3*2+num_dof*i+1] = 0.;
		LFS[num_dof*3*2+num_dof*i+2] = basis[i];
		LFS[num_dof*3*3+num_dof*i+0] = 0.;
		LFS[num_dof*3*3+num_dof*i+1] = 0.;
		LFS[num_dof*3*3+num_dof*i+2] = basis[i];
	}
}/*}}}*/
void           StressbalanceAnalysis::GetLSSAFS(IssmDouble* LFS,Element* element,Gauss* gauss_in){/*{{{*/
	/*
	 * Compute L  matrix. L=[L1 L2 L3] where Li is square and of size numdof.
	 * For node i, Li can be expressed in the actual coordinate system
	 * by:
	 *       Li=[ h    0 ]
	 *	 	      [ 0    h ]
	 *	 	      [ h    0 ]
	 *	 	      [ 0    h ]
	 *	 	      [ h    0 ]
	 *	 	      [ 0    h ]
	 *	 	      [ h    0 ]
	 *	 	      [ 0    h ]
	 * where h is the interpolation function for node i.
	 */

	int num_dof=2;
	IssmDouble basis[3];

	/*Cast gauss to GaussPenta*/
	_assert_(gauss_in->Enum()==GaussPentaEnum);
	GaussPenta* gauss = xDynamicCast<GaussPenta*>(gauss_in);

	/*Get basis in actual coordinate system: */
	basis[0]=gauss->coord1*(1-gauss->coord4)/2.0;
	basis[1]=gauss->coord2*(1-gauss->coord4)/2.0;
	basis[2]=gauss->coord3*(1-gauss->coord4)/2.0;

	/*Build LFS: */
	for(int i=0;i<3;i++){
		LFS[num_dof*3*0+num_dof*i+0] = basis[i];
		LFS[num_dof*3*0+num_dof*i+1] = 0;
		LFS[num_dof*3*1+num_dof*i+0] = 0;
		LFS[num_dof*3*1+num_dof*i+1] = basis[i];
		LFS[num_dof*3*2+num_dof*i+0] = basis[i];
		LFS[num_dof*3*2+num_dof*i+1] = 0;
		LFS[num_dof*3*3+num_dof*i+0] = 0;
		LFS[num_dof*3*3+num_dof*i+1] = basis[i];
		LFS[num_dof*3*4+num_dof*i+0] = basis[i];
		LFS[num_dof*3*4+num_dof*i+1] = 0;
		LFS[num_dof*3*5+num_dof*i+0] = 0;
		LFS[num_dof*3*5+num_dof*i+1] = basis[i];
		LFS[num_dof*3*6+num_dof*i+0] = basis[i];
		LFS[num_dof*3*6+num_dof*i+1] = 0;
		LFS[num_dof*3*7+num_dof*i+0] = 0;
		LFS[num_dof*3*7+num_dof*i+1] = basis[i];
	}
}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolutionHOFS(IssmDouble* solution,Element* element){/*{{{*/

	int         i;
	IssmDouble  rho_ice,g,FSreconditioning;
	int*        doflistHO  = NULL;
	int*        doflistFSv = NULL;
	int*        doflistFSp = NULL;

	/*Only works with Penta for now*/
	if(element->ObjectEnum()!=PentaEnum) _error_("Coupling not supported for "<<EnumToStringx(element->ObjectEnum()));

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes  = 6;
	int numdofHO  = 6*2;
	int numdofFSv = 6*3;
	int numdofFSp = 6;

	/*Fetch dof list and allocate solution vectors*/
	element->GetDofListLocal(&doflistFSv,FSvelocityEnum,GsetEnum);
	element->GetDofListLocal(&doflistHO, HOApproximationEnum, GsetEnum);
	element->GetDofListLocalPressure(&doflistFSp,GsetEnum);
	IssmDouble* HOvalues  = xNew<IssmDouble>(numdofHO);
	IssmDouble* FSvalues  = xNew<IssmDouble>(numdofFSv+numdofFSp);
	IssmDouble* vx        = xNew<IssmDouble>(numnodes);
	IssmDouble* vy        = xNew<IssmDouble>(numnodes);
	IssmDouble* vz        = xNew<IssmDouble>(numnodes);
	IssmDouble* vzHO      = xNew<IssmDouble>(numnodes);
	IssmDouble* vzFS      = xNew<IssmDouble>(numnodes);
	IssmDouble* vel       = xNew<IssmDouble>(numnodes);
	IssmDouble* pressure  = xNew<IssmDouble>(numnodes);

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(2*numnodes);
	for(i=0;i<numnodes;i++) cs_list[i]          = XYZEnum;
	for(i=0;i<numnodes;i++) cs_list[numnodes+i] = PressureEnum;

	/*Use the dof list to index into the solution vector: */
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	for(i=0;i<numdofHO ;i++) HOvalues[i]=solution[doflistHO[i]];
	for(i=0;i<numdofFSv;i++) FSvalues[i]=solution[doflistFSv[i]];
	for(i=0;i<numdofFSp;i++) FSvalues[numdofFSv+i]=solution[doflistFSp[i]];

	/*Transform solution in Cartesian Space*/
	element->TransformSolutionCoord(FSvalues,2*numnodes,cs_list);
	element->TransformSolutionCoord(HOvalues,numnodes,XYEnum);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	for(i=0;i<numnodes;i++){
		vx[i]       = FSvalues[i*3+0]+HOvalues[i*2+0];
		vy[i]       = FSvalues[i*3+1]+HOvalues[i*2+1];
		vzFS[i]     = FSvalues[i*3+2];
		pressure[i] = FSvalues[numnodes*3+i]*FSreconditioning;

		/*Check solution*/
		if(xIsNan<IssmDouble>(vx[i]))       _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vx[i]))       _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(vy[i]))       _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vy[i]))       _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(vzFS[i]))     _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vzFS[i]))     _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(pressure[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(pressure[i])) _error_("Inf found in solution vector");
	}

	/*Get Vz and compute vel*/
	element->GetInputListOnVertices(vzHO,VzHOEnum);
	for(i=0;i<numnodes;i++){
		vz[i] = vzHO[i]+vzFS[i];
		vel[i]= sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	}

	/*Add vx and vy as inputs to element: */
	element->AddInput(VxEnum,vx,P1Enum);
	element->AddInput(VyEnum,vy,P1Enum);
	element->AddInput(VzEnum,vz,P1Enum);
	element->AddInput(VzFSEnum,vzFS,P1Enum);
	element->AddInput(VelEnum,vel,P1Enum);
	element->AddInput(PressureEnum,pressure,P1Enum);

	/*Free resources:*/
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vz);
	xDelete<IssmDouble>(vzHO);
	xDelete<IssmDouble>(vzFS);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(FSvalues);
	xDelete<IssmDouble>(HOvalues);
	xDelete<int>(doflistFSp);
	xDelete<int>(doflistFSv);
	xDelete<int>(doflistHO);
	xDelete<int>(cs_list);
}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolutionSSAFS(IssmDouble* solution,Element* element){/*{{{*/

	int         i;
	IssmDouble  rho_ice,g,FSreconditioning;
	int*        doflistSSA  = NULL;
	int*        doflistFSv = NULL;
	int*        doflistFSp = NULL;

	/*we have to add results of this element for FS and results from the element
	 * at base for SSA, so we need to have the pointer toward the basal element*/
	Element* basalelement=element->GetBasalElement();
	if(basalelement->ObjectEnum()!=PentaEnum){
		_error_("Coupling not supported for "<<EnumToStringx(basalelement->ObjectEnum()));
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes  = 6;
	int numdof2d  = numnodes;
	int numdofSSA = 6*2;
	int numdofFSv = 6*3;
	int numdofFSp = 6;

	/*Fetch dof list and allocate solution vectors*/
	element->GetDofListLocal(&doflistFSv,FSvelocityEnum,GsetEnum);
	element->GetDofListLocalPressure(&doflistFSp,GsetEnum);
	basalelement->GetDofListLocal(&doflistSSA, SSAApproximationEnum, GsetEnum);
	IssmDouble* SSAvalues  = xNew<IssmDouble>(numdofSSA);
	IssmDouble* FSvalues  = xNew<IssmDouble>(numdofFSv+numdofFSp);
	IssmDouble* vx        = xNew<IssmDouble>(numnodes);
	IssmDouble* vy        = xNew<IssmDouble>(numnodes);
	IssmDouble* vz        = xNew<IssmDouble>(numnodes);
	IssmDouble* vzSSA      = xNew<IssmDouble>(numnodes);
	IssmDouble* vzFS      = xNew<IssmDouble>(numnodes);
	IssmDouble* vel       = xNew<IssmDouble>(numnodes);
	IssmDouble* pressure  = xNew<IssmDouble>(numnodes);

	/*Prepare coordinate system list*/
	int* cs_list = xNew<int>(2*numnodes);
	for(i=0;i<numnodes;i++) cs_list[i]          = XYZEnum;
	for(i=0;i<numnodes;i++) cs_list[numnodes+i] = PressureEnum;

	/*Use the dof list to index into the solution vector: */
	element->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);
	for(i=0;i<numdof2d;i++){
		SSAvalues[i]          = solution[doflistSSA[i]];
		SSAvalues[i+numdof2d] = solution[doflistSSA[i]];
	}
	for(i=0;i<numdofFSv;i++) FSvalues[i]=solution[doflistFSv[i]];
	for(i=0;i<numdofFSp;i++) FSvalues[numdofFSv+i]=solution[doflistFSp[i]];

	/*Transform solution in Cartesian Space*/
	element->TransformSolutionCoord(FSvalues,2*numnodes,cs_list);
	element->TransformSolutionCoord(SSAvalues,numnodes,XYEnum);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */

	for(i=0;i<numnodes;i++){
		vx[i]       = FSvalues[i*3+0]+SSAvalues[i*2+0];
		vy[i]       = FSvalues[i*3+1]+SSAvalues[i*2+1];
		vzFS[i]     = FSvalues[i*3+2];
		pressure[i] = FSvalues[numnodes*3+i]*FSreconditioning;

		/*Check solution*/
		if(xIsNan<IssmDouble>(vx[i]))       _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vx[i]))       _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(vy[i]))       _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vy[i]))       _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(vzFS[i]))     _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vzFS[i]))     _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(pressure[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(pressure[i])) _error_("Inf found in solution vector");
	}

	/*Get Vz and compute vel*/
	element->GetInputListOnVertices(vzSSA,VzSSAEnum);
	for(i=0;i<numnodes;i++){
		vz[i] = vzSSA[i]+vzFS[i];
		vel[i]= sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	}

	/*Add vx and vy as inputs to element: */
	element->AddInput(VxEnum,vx,P1Enum);
	element->AddInput(VyEnum,vy,P1Enum);
	element->AddInput(VzEnum,vz,P1Enum);
	element->AddInput(VzFSEnum,vzFS,P1Enum);
	element->AddInput(VelEnum,vel,P1Enum);
	element->AddInput(PressureEnum,pressure,P1Enum);

	/*Free resources:*/
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vz);
	xDelete<IssmDouble>(vzSSA);
	xDelete<IssmDouble>(vzFS);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(FSvalues);
	xDelete<IssmDouble>(SSAvalues);
	xDelete<int>(doflistFSp);
	xDelete<int>(doflistFSv);
	xDelete<int>(doflistSSA);
	xDelete<int>(cs_list);
}/*}}}*/
void           StressbalanceAnalysis::InputUpdateFromSolutionSSAHO(IssmDouble* solution,Element* element){/*{{{*/

	int         i,domaintype;
	IssmDouble  rho_ice,g;
	int*        SSAdoflist = NULL;
	int*        HOdoflist  = NULL;
	IssmDouble* xyz_list   = NULL;

	/*we have to add results of this element for HO and results from the element
	 * at base for SSA, so we need to have the pointer toward the basal element*/
	Element* basalelement=element->GetBasalElement();
	if(basalelement->ObjectEnum()!=PentaEnum){
		_error_("Coupling not supported for "<<EnumToStringx(basalelement->ObjectEnum()));
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*2;
	int numdof2d = numnodes;

	/*Fetch dof list and allocate solution vectors*/
	basalelement->GetDofListLocal(&SSAdoflist,SSAApproximationEnum,GsetEnum);
	element     ->GetDofListLocal(&HOdoflist, HOApproximationEnum, GsetEnum);
	IssmDouble* HOvalues  = xNew<IssmDouble>(numdof);
	IssmDouble* SSAvalues = xNew<IssmDouble>(numdof);
	IssmDouble* vx        = xNew<IssmDouble>(numnodes);
	IssmDouble* vy        = xNew<IssmDouble>(numnodes);
	IssmDouble* vz        = xNew<IssmDouble>(numnodes);
	IssmDouble* vel       = xNew<IssmDouble>(numnodes);
	IssmDouble* pressure  = xNew<IssmDouble>(numnodes);
	IssmDouble* surface   = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numdof2d;i++){
		HOvalues[i]  = solution[HOdoflist[i]];
		SSAvalues[i] = solution[SSAdoflist[i]];
	}
	for(i=numdof2d;i<numdof;i++){
		HOvalues[i]  = solution[HOdoflist[i]];
		SSAvalues[i] = SSAvalues[i-numdof2d];
	}

	/*Transform solution in Cartesian Space*/
	basalelement->TransformSolutionCoord(SSAvalues,XYEnum);
	element->TransformSolutionCoord(HOvalues,XYEnum);

	/*Ok, we have vx and vy in values, fill in vx and vy arrays: */
	for(i=0;i<numnodes;i++){
		vx[i]=SSAvalues[i*2+0]+HOvalues[i*2+0];
		vy[i]=SSAvalues[i*2+1]+HOvalues[i*2+1];

		/*Check solution*/
		if(xIsNan<IssmDouble>(vx[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vx[i])) _error_("Inf found in solution vector");
		if(xIsNan<IssmDouble>(vy[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(vy[i])) _error_("Inf found in solution vector");
	}

	/*Get Vz and compute vel*/
	element->GetInputListOnNodes(&vz[0],VzEnum,0.);
	for(i=0;i<numnodes;i++) vel[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

	/*For pressure: we have not computed pressure in this analysis, for this element. We are in 2D,
	 *so the pressure is just the pressure at the bedrock: */
	rho_ice = element->FindParam(MaterialsRhoIceEnum);
	g       = element->FindParam(ConstantsGEnum);
	element->GetVerticesCoordinates(&xyz_list);
	element->GetInputListOnNodes(&surface[0],SurfaceEnum);
	for(i=0;i<numnodes;i++) pressure[i]=rho_ice*g*(surface[i]-xyz_list[i*3+2]);

	/*Add vx and vy as inputs to element: */
	element->AddInput(VxEnum,vx,P1Enum);
	element->AddInput(VyEnum,vy,P1Enum);
	element->AddInput(VelEnum,vel,P1Enum);
	element->AddInput(PressureEnum,pressure,P1Enum);

	/*Free resources:*/
	xDelete<IssmDouble>(surface);
	xDelete<IssmDouble>(pressure);
	xDelete<IssmDouble>(vel);
	xDelete<IssmDouble>(vz);
	xDelete<IssmDouble>(vy);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(SSAvalues);
	xDelete<IssmDouble>(HOvalues);
	xDelete<int>(SSAdoflist);
	xDelete<int>(HOdoflist);
}/*}}}*/
