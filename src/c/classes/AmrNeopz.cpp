/*!\file AdaptiveMeshrefinement.cpp
 * \brief: implementation of the adaptive mesh refinement tool based on NeoPZ library: github.com/labmec/neopz
 */

#ifdef HAVE_CONFIG_H
    #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./AmrNeopz.h"

/*Includes*/
/*{{{*/
/*Common includes*/
#include <iostream>
#include <fstream>
#include <string>
#include <climits>
#include <cfloat>

/*NeoPZ includes*/
#include <pz_config.h>
#include <pzreal.h>
#include <pzvec.h>
#include <pzeltype.h>

#include <TPZRefPatternTools.h>
#include <TPZRefPatternDataBase.h>
#include <TPZRefPattern.h>

#include <tpzchangeel.h>
#include <TPZGeoElement.h>
#include <pzreftriangle.h>
#include <pzgeotriangle.h>
#include <tpzgeoelrefpattern.h>
#include <pzgraphmesh.h>
#include <TPZVTKGeoMesh.h>
/*}}}*/

using namespace pzgeom;

/*Constructor, copy, clean up and destructor*/
AmrNeopz::AmrNeopz(){/*{{{*/

	/*Set pointers to NULL*/
	this->fathermesh						= NULL;
	this->previousmesh					= NULL;
	this->refinement_type				= -1;
	this->level_max						= -1;
	this->gradation						= -1;
	this->lag								= -1;
   this->groundingline_distance		= -1;
	this->icefront_distance				= -1;
	this->thicknesserror_threshold	= -1;
	this->deviatoricerror_threshold	= -1;
	this->deviatoricerror_maximum		= -1;
	this->thicknesserror_maximum		= -1;
	this->sid2index.clear();
	this->index2sid.clear();
	this->specialelementsindex.clear();
	this->x									= NULL;
	this->y									= NULL;
	this->elementslist					= NULL;
	this->numberofvertices				= -1;
	this->numberofelements				= -1;
}
/*}}}*/
AmrNeopz::AmrNeopz(const AmrNeopz &cp){/*{{{*/
	this->Initialize(); 
	this->operator =(cp);
}
/*}}}*/
AmrNeopz & AmrNeopz::operator =(const AmrNeopz &cp){/*{{{*/

	/*Clean all attributes*/
	this->CleanUp();
	/*Copy all data*/
	this->fathermesh						= new TPZGeoMesh(*cp.fathermesh);
	this->previousmesh					= new TPZGeoMesh(*cp.previousmesh);
	this->refinement_type				= cp.refinement_type;
	this->level_max						= cp.level_max;
	this->gradation						= cp.gradation;
	this->lag								= cp.lag;
   this->groundingline_distance		= cp.groundingline_distance;
	this->icefront_distance				= cp.icefront_distance;
	this->thicknesserror_threshold	= cp.thicknesserror_threshold;
	this->deviatoricerror_threshold	= cp.deviatoricerror_threshold;
	this->deviatoricerror_maximum		= cp.deviatoricerror_maximum;
	this->thicknesserror_maximum		= cp.thicknesserror_maximum;
	this->sid2index.clear();
	this->sid2index.resize(cp.sid2index.size());
	for(int i=0;i<cp.sid2index.size();i++) this->sid2index[i]=cp.sid2index[i];
	this->index2sid.clear();
	this->index2sid.resize(cp.index2sid.size());
	for(int i=0;i<cp.index2sid.size();i++) this->index2sid[i]=cp.index2sid[i];
	this->specialelementsindex.clear();
	this->specialelementsindex.resize(cp.specialelementsindex.size());
	for(int i=0;i<cp.specialelementsindex.size();i++) this->specialelementsindex[i]=cp.specialelementsindex[i];

	return *this;
}
/*}}}*/
AmrNeopz::~AmrNeopz(){/*{{{*/
	int writemesh = 0;//only to restart
	if(writemesh) this->WriteMesh();
	this->CleanUp();
	gRefDBase.clear();
}
/*}}}*/
void AmrNeopz::CleanUp(){/*{{{*/

	/*Verify and delete all data*/
	if(this->fathermesh)    delete this->fathermesh;
	if(this->previousmesh)  delete this->previousmesh;
	if(this->x)					xDelete<IssmDouble>(this->x);
	if(this->y)					xDelete<IssmDouble>(this->y);
	if(this->elementslist)	xDelete<int>(this->elementslist);
	this->refinement_type				= -1;
	this->level_max						= -1;
	this->gradation						= -1;
	this->lag								= -1;
   this->groundingline_distance		= -1;
	this->icefront_distance				= -1;
	this->thicknesserror_threshold	= -1;
	this->deviatoricerror_threshold	= -1;
	this->deviatoricerror_maximum		= -1;
	this->thicknesserror_maximum		= -1;
	this->numberofvertices				= -1;
	this->numberofelements				= -1;
	this->sid2index.clear();
	this->index2sid.clear();
	this->specialelementsindex.clear();
}
/*}}}*/

/*Mesh refinement methods*/
void AmrNeopz::SetMesh(int** elementslist_in,IssmDouble** x_in,IssmDouble** y_in,int* numberofvertices_in,int* numberofelements_in){/*{{{*/

   /*Delete previous mesh and keep the entire mesh*/
   if(this->elementslist) xDelete<int>(this->elementslist);
   if(this->x) xDelete<IssmDouble>(this->x);
   if(this->y) xDelete<IssmDouble>(this->y);

   this->elementslist      = *elementslist_in;
   this->x                 = *x_in;
   this->y                 = *y_in;
   this->numberofvertices  = *numberofvertices_in;
   this->numberofelements  = *numberofelements_in;
}/*}}}*/
void AmrNeopz::GetMesh(int** elementslist_out,IssmDouble** x_out,IssmDouble** y_out,int* numberofvertices_out,int* numberofelements_out){/*{{{*/

   /*Get the entire mesh*/
   *elementslist_out    = this->elementslist;
   *x_out               = this->x;
   *y_out               = this->y;
   *numberofvertices_out= this->numberofvertices;
   *numberofelements_out= this->numberofelements;
}/*}}}*/
void AmrNeopz::ExecuteRefinement(double* gl_distance,double* if_distance,double* deviatoricerror,double* thicknesserror,int** pdatalist,double** pxylist,int** pelementslist){/*{{{*/

	/*IMPORTANT! pelementslist are in Matlab indexing*/
	/*NEOPZ works only in C indexing*/
	if(!this->fathermesh || !this->previousmesh) _error_("Impossible to execute refinement: fathermesh or previousmesh is NULL!\n");
	if(this->refinement_type!=0 && this->refinement_type!=1) _error_("Impossible to execute refinement: refinement type is not defined!\n");

	/*Input verifications*/
	if(this->deviatoricerror_threshold>0	&& !deviatoricerror) _error_("deviatoricerror is NULL!\n");
	if(this->thicknesserror_threshold>0		&& !thicknesserror)	_error_("thicknesserror is NULL!\n");
	if(this->groundingline_distance>0		&& !gl_distance)		_error_("gl_distance is NULL!\n");
	if(this->icefront_distance>0				&& !if_distance)		_error_("if_distance is NULL!\n");
	/*Attributes verifications*/
	if(this->deviatoricerror_threshold>0	&& this->deviatoricerror_groupthreshold<DBL_EPSILON)	_error_("group threshold is too small!");
	if(this->thicknesserror_threshold>0		&& this->thicknesserror_groupthreshold<DBL_EPSILON)	_error_("group threshold is too small!");

	/*Intermediaries*/
	bool verbose=VerboseSolution();

	/*Execute refinement*/
	this->RefineMeshOneLevel(verbose,gl_distance,if_distance,deviatoricerror,thicknesserror);

	/*Get new geometric mesh in ISSM data structure*/
	this->GetMesh(this->previousmesh,pdatalist,pxylist,pelementslist);

	/*Verify the new geometry*/
	this->CheckMesh(pdatalist,pxylist,pelementslist);
}
/*}}}*/
void AmrNeopz::RefineMeshOneLevel(bool &verbose,double* gl_distance,double* if_distance,double* deviatoricerror,double* thicknesserror){/*{{{*/

	/*Intermediaries*/
	int nelem							=-1;
	int side2D							= 6;
	int sid								=-1;
	int count							=-1;
	int criteria						=-1;
	int numberofcriteria				=-1;
	int nconformelements				= this->sid2index.size();
	double gl_distance_h				=-1;
	double gl_distance_hmax			= this->groundingline_distance;
	double if_distance_h				=-1;
	double if_distance_hmax			= this->icefront_distance;
	double gl_groupdistance			=-1;
	double if_groupdistance			=-1;
	double d_maxerror					=-1;
	double t_maxerror					=-1;
	double deviatoric_grouperror	=-1;
	double thickness_grouperror	=-1;
	TPZGeoMesh* gmesh					= NULL; 
	TPZVec<REAL> qsi(2,0.),cp(3,0.);
	TPZVec<TPZGeoEl*> sons;
	std::vector<int> index;

	/*Calculate the number of criteria{{{*/
	numberofcriteria=0;
	if(this->deviatoricerror_threshold>0)	numberofcriteria++;
	if(this->thicknesserror_threshold>0)	numberofcriteria++;
	if(this->groundingline_distance>0)		numberofcriteria++;
	if(this->icefront_distance>0)				numberofcriteria++;
	/*}}}*/

	/*Calculate the maximum of the estimators, if requested{{{*/
	if(this->deviatoricerror_threshold>0 && this->deviatoricerror_maximum<DBL_EPSILON){ 
		for(int i=0;i<nconformelements;i++) this->deviatoricerror_maximum=max(this->deviatoricerror_maximum,deviatoricerror[i]);
	}
	if(this->thicknesserror_threshold>0 && this->thicknesserror_maximum<DBL_EPSILON){
		for(int i=0;i<nconformelements;i++) this->thicknesserror_maximum=max(this->thicknesserror_maximum,thicknesserror[i]);
	}
	/*}}}*/

	/*First, verify if special elements have min distance or high errors{{{*/
	gmesh=this->previousmesh;
	for(int i=0;i<this->specialelementsindex.size();i++){
		if(this->specialelementsindex[i]==-1) _error_("index is -1!\n");
		if(!gmesh->Element(this->specialelementsindex[i])) _error_("element is null!\n");
		if(!gmesh->Element(this->specialelementsindex[i])->Father()) _error_("father is null!\n");
		if(gmesh->Element(this->specialelementsindex[i])->HasSubElement()) _error_("special element has sub elements!\n");
		sons.clear();
		gmesh->Element(this->specialelementsindex[i])->Father()->GetHigherSubElements(sons);
		/*Limits*/
		gl_distance_h	= gl_distance_hmax*std::pow(this->gradation,this->level_max-gmesh->Element(this->specialelementsindex[i])->Level());
		if_distance_h	= if_distance_hmax*std::pow(this->gradation,this->level_max-gmesh->Element(this->specialelementsindex[i])->Level());
		d_maxerror		= this->deviatoricerror_threshold*this->deviatoricerror_maximum;
		t_maxerror		= this->thicknesserror_threshold*this->thicknesserror_maximum;
		/*Calculate the distance and error of the group (sons)*/
		gl_groupdistance=INFINITY;if_groupdistance=INFINITY;deviatoric_grouperror=0;thickness_grouperror=0;
		for(int s=0;s<sons.size();s++){
			sid=this->index2sid[sons[s]->Index()];
			if(sid<0) continue;
			if(this->groundingline_distance>0)		gl_groupdistance=std::min(gl_groupdistance,gl_distance[sid]); 
			if(this->icefront_distance>0)				if_groupdistance=std::min(if_groupdistance,if_distance[sid]); 
			if(this->deviatoricerror_threshold>0)	deviatoric_grouperror+=deviatoricerror[sid];
			if(this->thicknesserror_threshold>0)	thickness_grouperror+=thicknesserror[sid];
		}	
		criteria=0;
		if(this->groundingline_distance>0		&& gl_groupdistance<gl_distance_h+DBL_EPSILON)		criteria++;
		if(this->icefront_distance>0				&& if_groupdistance<if_distance_h+DBL_EPSILON)		criteria++;
		if(this->deviatoricerror_threshold>0	&& deviatoric_grouperror>d_maxerror-DBL_EPSILON)	criteria++;
		if(this->thicknesserror_threshold>0		&& thickness_grouperror>t_maxerror-DBL_EPSILON)		criteria++;
		/*Finally, it keeps the father index if it must be refine*/
		if(criteria) index.push_back(gmesh->Element(this->specialelementsindex[i])->FatherIndex());
	}
	/*}}}*/

	/*Now, detele the special elements{{{*/
	if(this->refinement_type==1) this->DeleteSpecialElements(verbose,gmesh);
	else this->specialelementsindex.clear();
	/*}}}*/

	/*Set the mesh and delete previousmesh if refinement type is 0{{{*/
	if(this->refinement_type==0){
		delete this->previousmesh;	
		gmesh=this->fathermesh;
	}
	/*}}}*/

	/*Unrefinement process: loop over level of refinements{{{*/
	if(verbose) _printf_("\tunrefinement process...\n");
	if(verbose) _printf_("\ttotal: ");
	count=0;

	nelem=gmesh->NElements();//must keep here
	for(int i=0;i<nelem;i++){
		if(!gmesh->Element(i)) continue;
		if(gmesh->Element(i)->MaterialId()!=this->GetElemMaterialID()) continue;
		if(gmesh->Element(i)->HasSubElement()) continue;
		if(gmesh->Element(i)->Level()==0) continue;
		if(!gmesh->Element(i)->Father()) _error_("father is NULL!\n");
		/*Limits with lag*/
		gl_distance_h = this->lag*gl_distance_hmax*std::pow(this->gradation,this->level_max-gmesh->Element(i)->Level());
		if_distance_h = this->lag*if_distance_hmax*std::pow(this->gradation,this->level_max-gmesh->Element(i)->Level());
		d_maxerror	  = this->deviatoricerror_groupthreshold*this->deviatoricerror_maximum;
		t_maxerror	  = this->thicknesserror_groupthreshold*this->thicknesserror_maximum;
		/*Get the sons of the father (sibilings)*/	
		sons.clear();
		gmesh->Element(i)->Father()->GetHigherSubElements(sons);
		if(sons.size()!=4) continue;//delete just group of 4 elements. This avoids big holes in the mesh
		/*Find the minimal distance and the error of the group*/	
		gl_groupdistance=INFINITY;if_groupdistance=INFINITY;deviatoric_grouperror=0;thickness_grouperror=0;
		for(int s=0;s<sons.size();s++){
			sid=this->index2sid[sons[s]->Index()];
			/*Verify if this group have solutions*/
			if(sid<0){gl_groupdistance=INFINITY;if_groupdistance=INFINITY;deviatoric_grouperror=INFINITY;thickness_grouperror=INFINITY;continue;} 
			/*Distance and error of the group*/
			if(this->groundingline_distance>0)		gl_groupdistance=std::min(gl_groupdistance,gl_distance[sid]); 
			if(this->icefront_distance>0)				if_groupdistance=std::min(if_groupdistance,if_distance[sid]); 
			if(this->deviatoricerror_threshold>0)	deviatoric_grouperror+=deviatoricerror[sid]; 
			if(this->thicknesserror_threshold>0)	thickness_grouperror+=thicknesserror[sid]; 
		}
		/*Verify the criteria*/
		criteria=0;
		if(this->groundingline_distance>0		&& gl_groupdistance>gl_distance_h-DBL_EPSILON)		criteria++;
		if(this->icefront_distance>0				&& if_groupdistance>if_distance_h-DBL_EPSILON)		criteria++;
		if(this->deviatoricerror_threshold>0	&& deviatoric_grouperror<d_maxerror+DBL_EPSILON)	criteria++;
		if(this->thicknesserror_threshold>0		&& thickness_grouperror<t_maxerror+DBL_EPSILON)		criteria++;
		/*Now, if the group attends the criteria, unrefine it*/
		if(criteria==numberofcriteria){ 
			gmesh->Element(i)->Father()->ResetSubElements(); count++;
			for(int s=0;s<sons.size();s++){this->index2sid[sons[s]->Index()]=-1;gmesh->DeleteElement(sons[s],sons[s]->Index());}
		}
	}
	if(verbose) _printf_(""<<count<<"\n");
	/*Adjust the connectivities before continue*/
	//gmesh->BuildConnectivity(); this is not necessary
	/*}}}*/

	/*Refinement process: loop over level of refinements{{{*/
	if(verbose) _printf_("\trefinement process (level max = "<<this->level_max<<")\n");
	if(verbose) _printf_("\ttotal: ");
	count=0;
	nelem=gmesh->NElements();//must keep here
	for(int i=0;i<nelem;i++){
		if(!gmesh->Element(i)) continue;
		if(gmesh->Element(i)->MaterialId()!=this->GetElemMaterialID()) continue;
		if(gmesh->Element(i)->HasSubElement()) continue;
		if(gmesh->Element(i)->Level()==this->level_max) continue;
		/*Verify if this element has solutions*/
		sid=this->index2sid[gmesh->Element(i)->Index()];
		if(sid<0) continue;
		/*Set the distance for level h*/
		gl_distance_h	= gl_distance_hmax*std::pow(this->gradation,this->level_max-(gmesh->Element(i)->Level()+1));//+1: current element level is <level_max
		if_distance_h	= if_distance_hmax*std::pow(this->gradation,this->level_max-(gmesh->Element(i)->Level()+1));//+1: current element level is <level_max
		d_maxerror		= this->deviatoricerror_threshold*this->deviatoricerror_maximum;
		t_maxerror		= this->thicknesserror_threshold*this->thicknesserror_maximum;
		/*Verify distance and error of the element, if requested*/
		criteria=0;
		if(this->groundingline_distance>0		&& gl_distance[sid]<gl_distance_h+DBL_EPSILON)	criteria++; 
		if(this->icefront_distance>0				&& if_distance[sid]<if_distance_h+DBL_EPSILON)	criteria++; 
		if(this->deviatoricerror_threshold>0	&& deviatoricerror[sid]>d_maxerror-DBL_EPSILON)	criteria++; 
		if(this->thicknesserror_threshold>0		&& thicknesserror[sid]>t_maxerror-DBL_EPSILON)	criteria++; 
		/*Now, if it attends any criterion, keep the element index to refine in next step*/
		if(criteria)index.push_back(i);
	}
	/*Now, refine the elements*/
	for(int i=0;i<index.size();i++){ 
		if(!gmesh->Element(index[i])) DebugStop();
		if(!gmesh->Element(index[i])->HasSubElement()){gmesh->Element(index[i])->Divide(sons);count++;}
	}
	if(verbose) _printf_(""<<count<<"\n");
	/*Adjust the connectivities before continue*/
	//gmesh->BuildConnectivity();//this is not necessary
	/*}}}*/

	/*Now, apply smoothing and insert special elements to avoid hanging nodes{{{*/
	this->RefineMeshWithSmoothing(verbose,gmesh);
	if(this->refinement_type==0) this->previousmesh=this->CreateRefPatternMesh(gmesh);//in this case, gmesh==this->fathermesh
	gmesh=this->previousmesh;//previous mesh is always refined to avoid hanging nodes
	this->RefineMeshToAvoidHangingNodes(verbose,gmesh);
	/*}}}*/
}
/*}}}*/
int AmrNeopz::VerifyRefinementType(TPZGeoEl* geoel,TPZGeoMesh* gmesh){/*{{{*/

	/*
	 * 0 : no refinement
	 * 1 : special refinement (to avoid hanging nodes)
	 * 2 : uniform refinment
	 * */
	if(!geoel) _error_("geoel is NULL!\n");

	/*Output*/
	int type=0;

	/*Intermediaries*/
	TPZVec<TPZGeoEl*> sons;

	/*Loop over neighboors (sides 3, 4 and 5)*/
	for(int j=3;j<6;j++){
		if(!gmesh->Element(geoel->NeighbourIndex(j))->HasSubElement()) continue;
		sons.clear();
		gmesh->Element(geoel->NeighbourIndex(j))->GetHigherSubElements(sons);
		if(sons.size()) type++; //if neighbour was refined
		if(sons.size()>4) type++; //if neighbour's level is > element level+1
		if(type>1) break;
	}

	/*Verify and return*/
	if(type>1) type=2;

	return type;
}
/*}}}*/
void AmrNeopz::RefineMeshWithSmoothing(bool &verbose,TPZGeoMesh* gmesh){/*{{{*/

	/*Intermediaries*/
	int nelem		=-1;
	int count		=-1;
	int type			=-1;
	int typecount	=-1;

	TPZVec<TPZGeoEl*> sons;

	/*Refinement process: loop over level of refinements*/
	if(verbose) _printf_("\tsmoothing process (level max = "<<this->level_max<<")\n");
	if(verbose) _printf_("\ttotal: ");

	count=1;

	while(count>0){
		count=0;
		nelem=gmesh->NElements();//must keep here
		for(int i=0;i<nelem;i++){
			if(!gmesh->Element(i)) continue;
			if(gmesh->Element(i)->MaterialId()!=this->GetElemMaterialID()) continue;
			if(gmesh->Element(i)->HasSubElement()) continue;
			if(gmesh->Element(i)->Level()==this->level_max) continue;
			/*loop over neighboors (sides 3, 4 and 5). Important: neighbours has the same dimension of the element*/
			type=this->VerifyRefinementType(gmesh->Element(i),gmesh);
			if(type<2){
				typecount=0;
				for(int j=3;j<6;j++){
					if(gmesh->Element(gmesh->Element(i)->NeighbourIndex(j))->HasSubElement()) continue;
					if(gmesh->Element(i)->NeighbourIndex(j)==i) typecount++;//neighbour==this element, element at the border
					if(this->VerifyRefinementType(gmesh->Element(gmesh->Element(i)->NeighbourIndex(j)),gmesh)==1) typecount++;
					if(typecount>1 && type==1) type=2;
					else if(typecount>2 && type==0) type=2;
					if(type==2) break;
				}
			}

			/*refine the element if requested*/
			if(type==2){gmesh->Element(i)->Divide(sons);	count++;}
		}
		if(verbose){
			if(count==0) _printf_(""<<count<<"\n");
			else _printf_(""<<count<<", ");
		}
		/*Adjust the connectivities before continue*/
		//gmesh->BuildConnectivity();//this is not necessary
	}
}
/*}}}*/
void AmrNeopz::RefineMeshToAvoidHangingNodes(bool &verbose,TPZGeoMesh* gmesh){/*{{{*/

	/*Now, insert special elements to avoid hanging nodes*/
	if(verbose) _printf_("\trefining to avoid hanging nodes (total: ");

	/*Intermediaries*/
	int nelem=-1;
	int count= 1;

	while(count>0){
		nelem=gmesh->NElements();//must keep here
		count=0;
		for(int i=0;i<nelem;i++){
			/*Get geometric element and verify if it has already been refined. Geoel may not have been previously refined*/
			TPZGeoEl * geoel=gmesh->Element(i);
			if(!geoel) continue;
			if(geoel->HasSubElement()) continue;
			if(geoel->MaterialId() != this->GetElemMaterialID()) continue;
			/*Get the refinement pattern for this element and refine it*/
			TPZAutoPointer<TPZRefPattern> refp=TPZRefPatternTools::PerfectMatchRefPattern(geoel);
			if(refp){
				/*Non-uniform refinement*/
				TPZVec<TPZGeoEl *> sons;
				geoel->SetRefPattern(refp);
				geoel->Divide(sons);
				count++;
				/*Keep the index of the special elements*/
				for(int j=0;j<sons.size();j++) this->specialelementsindex.push_back(sons[j]->Index());
			}
		}
		if(verbose){
			if(count==0) _printf_(""<<count<<")\n");
			else _printf_(""<<count<<", ");
		}
		//gmesh->BuildConnectivity();//this is not necessary
	}
}
/*}}}*/
void AmrNeopz::DeleteSpecialElements(bool &verbose,TPZGeoMesh* gmesh){/*{{{*/

	/*Intermediaries*/
	int count=0;

	if(verbose) _printf_("\tdelete "<<this->specialelementsindex.size()<<" special elements (total: ");
	for(int i=0;i<this->specialelementsindex.size();i++){
		if(this->specialelementsindex[i]==-1) continue;
		if(!gmesh->Element(this->specialelementsindex[i])) continue;
		if(gmesh->Element(this->specialelementsindex[i])->Father()) gmesh->Element(this->specialelementsindex[i])->Father()->ResetSubElements();
		gmesh->DeleteElement(gmesh->Element(this->specialelementsindex[i]),this->specialelementsindex[i]);
      this->index2sid[this->specialelementsindex[i]]=-1;
		count++;
	}
	if(verbose) _printf_(""<<count<<")\n");
	/*Cleanup*/
	this->specialelementsindex.clear();
	/*Adjust connectivities*/
	//gmesh->BuildConnectivity();//this is not necessary
}
/*}}}*/
void AmrNeopz::GetMesh(TPZGeoMesh* gmesh,int** pdata,double** pxy, int** pelements){/*{{{*/

	/* IMPORTANT! pelements are in Matlab indexing
	   NEOPZ works only in C indexing.
		This method cleans up and updated the this->sid2index
		and this->index2sid and fills in it with the new mesh.
		Avoid to call this method before Refinement Process.*/

	/*Intermediaries */
	long sid,nodeindex;
	int nconformelements,nconformvertices;
	int ntotalvertices		= gmesh->NNodes();
	int* newelements			= NULL;
	int* newdata				= NULL;
	double* newmeshXY			= NULL;
	TPZGeoEl* geoel			= NULL;
	long* vertex_index2sid 	= xNew<long>(ntotalvertices);
	this->index2sid.clear(); this->index2sid.resize(gmesh->NElements());
	this->sid2index.clear();

	/*Fill in the vertex_index2sid vector with non usual index value*/
	for(int i=0;i<gmesh->NNodes();i++) vertex_index2sid[i]=-1;

	/*Fill in the this->index2sid vector with non usual index value*/
	for(int i=0;i<gmesh->NElements();i++) this->index2sid[i]=-1;

	/*Get elements without sons and fill in the vertex_index2sid with used vertices (indexes) */
	sid=0;
	for(int i=0;i<gmesh->NElements();i++){//over gmesh elements index 
		geoel=gmesh->ElementVec()[i];
		if(!geoel) continue;
		if(geoel->HasSubElement()) continue;
		if(geoel->MaterialId() != this->GetElemMaterialID()) continue;
		this->sid2index.push_back(geoel->Index());//keep the element index
		this->index2sid[geoel->Index()]=this->sid2index.size()-1;//keep the element sid
		for(int j=0;j<this->GetNumberOfNodes();j++){
      	nodeindex=geoel->NodeIndex(j);
      	if(nodeindex<0) _error_("nodeindex is <0\n");
			if(vertex_index2sid[nodeindex]==-1){
      		vertex_index2sid[nodeindex]=sid; 
				sid++;
			}
      }	
	}

	/* Create new mesh structure and fill it */
	nconformelements	= (int)this->sid2index.size();
	nconformvertices	= (int)sid;
	newelements			= xNew<int>(nconformelements*this->GetNumberOfNodes());
	newmeshXY			= xNew<double>(nconformvertices*2);
	newdata				= xNew<int>(2);
	newdata[0]			= nconformvertices;
	newdata[1]			= nconformelements;

	for(int i=0;i<ntotalvertices;i++){//over the TPZNode index (fill in the ISSM vertices coords)
		sid=vertex_index2sid[i];
		if(sid==-1) continue;//skip this index (node no used)
		TPZVec<REAL> coords(3,0.);
		gmesh->NodeVec()[i].GetCoordinates(coords);
		newmeshXY[2*sid]		= coords[0]; // X
		newmeshXY[2*sid+1]	= coords[1]; // Y
	}

	for(int i=0;i<this->sid2index.size();i++){//over the sid (fill the ISSM elements)
		for(int j=0;j<this->GetNumberOfNodes();j++) {
			geoel	= gmesh->ElementVec()[this->sid2index[i]];
			sid	= vertex_index2sid[geoel->NodeIndex(j)];
			newelements[i*this->GetNumberOfNodes()+j]=(int)sid+1;//C to Matlab indexing
		}
		/*Verify the Jacobian determinant. If detJ<0, swap the 2 first postions:
		  a -> b
		  b -> a */
		double detJ,xa,xb,xc,ya,yb,yc;
		int a,b,c;

		a=newelements[i*this->GetNumberOfNodes()+0]-1;
		b=newelements[i*this->GetNumberOfNodes()+1]-1;
		c=newelements[i*this->GetNumberOfNodes()+2]-1;

		xa=newmeshXY[2*a]; ya=newmeshXY[2*a+1];
		xb=newmeshXY[2*b]; yb=newmeshXY[2*b+1];
		xc=newmeshXY[2*c]; yc=newmeshXY[2*c+1];

		detJ=(xb-xa)*(yc-ya)-(xc-xa)*(yb-ya);

		/*verify and swap, if necessary*/
		if(detJ<0) {
			newelements[i*this->GetNumberOfNodes()+0]=b+1;//a->b
			newelements[i*this->GetNumberOfNodes()+1]=a+1;//b->a
		}
	}

	/*Setting outputs*/
	*pdata		= newdata;
	*pxy		   = newmeshXY;
	*pelements	= newelements;

	/*Cleanup*/
	xDelete<long>(vertex_index2sid);
}
/*}}}*/
void AmrNeopz::Initialize(){/*{{{*/

	/* IMPORTANT! elements come in Matlab indexing
		NEOPZ works only in C indexing*/

	if(this->numberofvertices<=0) _error_("Impossible to create initial mesh: nvertices is <= 0!\n");
   if(this->numberofelements<=0) _error_("Impossible to create initial mesh: nelements is <= 0!\n");
	if(this->refinement_type!=0 && this->refinement_type!=1) _error_("Impossible to create initial mesh: refinement type is not defined!\n");

    /*Verify and creating initial mesh*/
   if(!this->x || !this->y || !this->elementslist) _error_("Mesh data structure is NULL!\n");
	if(this->fathermesh || this->previousmesh) _error_("Initial mesh already exists!\n");

   this->fathermesh = new TPZGeoMesh();
	this->fathermesh->NodeVec().Resize(this->numberofvertices);

	/*Set the vertices (geometric nodes in NeoPZ context)*/
	for(int i=0;i<this->numberofvertices;i++){  
      /*x,y,z coords*/
		TPZManVector<REAL,3> coord(3,0.);
      coord[0]= this->x[i];
      coord[1]= this->y[i];
      coord[2]= 0.;
      /*Insert in the mesh*/
      this->fathermesh->NodeVec()[i].SetCoord(coord);
		this->fathermesh->NodeVec()[i].SetNodeId(i);
	}

	/*Generate the elements*/
	int64_t index;
   const int mat = this->GetElemMaterialID();
   TPZManVector<int64_t> elem(this->GetNumberOfNodes(),0);
	this->index2sid.clear(); this->index2sid.resize(this->numberofelements);
   this->sid2index.clear();

	for(int i=0;i<this->numberofelements;i++){
		for(int j=0;j<this->GetNumberOfNodes();j++) elem[j]=this->elementslist[i*this->GetNumberOfNodes()+j]-1;//Convert Matlab to C indexing
      switch(this->GetNumberOfNodes()){
			case 3: this->fathermesh->CreateGeoElement(ETriangle,elem,mat,index,this->refinement_type);	break;
         default:	_error_("mesh not supported yet");
		}
      /*Define the element ID*/        
      this->fathermesh->ElementVec()[index]->SetId(i);
		/*Initialize sid2index and index2sid*/
		this->sid2index.push_back((int)index);
		this->index2sid[(int)index]=this->sid2index.size()-1;//keep the element sid
	}
   /*Build element and node connectivities*/
   this->fathermesh->BuildConnectivity();
	/*Set previous mesh*/
	if(this->refinement_type==1) this->previousmesh=new TPZGeoMesh(*this->fathermesh);
	else this->previousmesh=this->CreateRefPatternMesh(this->fathermesh); 
}
/*}}}*/
TPZGeoMesh* AmrNeopz::CreateRefPatternMesh(TPZGeoMesh* gmesh){/*{{{*/

	TPZGeoMesh *newgmesh = new TPZGeoMesh();
   newgmesh->CleanUp();

   int nnodes  = gmesh->NNodes();
	int nelem   = gmesh->NElements();
   int mat     = this->GetElemMaterialID();;
   int reftype = 1;
   int64_t index; 

	//nodes
	newgmesh->NodeVec().Resize(nnodes);
   for(int i=0;i<nnodes;i++) newgmesh->NodeVec()[i] = gmesh->NodeVec()[i];

   //elements
   for(int i=0;i<nelem;i++){
   	TPZGeoEl * geoel = gmesh->Element(i);

		if(!geoel){
			index=newgmesh->ElementVec().AllocateNewElement();
			newgmesh->ElementVec()[index] = NULL;
			continue;
		}

		TPZManVector<int64_t> elem(3,0);
      for(int j=0;j<3;j++) elem[j] = geoel->NodeIndex(j);

      newgmesh->CreateGeoElement(ETriangle,elem,mat,index,reftype);
		newgmesh->ElementVec()[index]->SetId(geoel->Id());

      TPZGeoElRefPattern<TPZGeoTriangle>* newgeoel = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle>*>(newgmesh->ElementVec()[index]);

      //old neighbourhood
      const int nsides = TPZGeoTriangle::NSides;
      TPZVec< std::vector<TPZGeoElSide> > neighbourhood(nsides);
      TPZVec<long> NodesSequence(0);
      for(int s = 0; s < nsides; s++){
      	neighbourhood[s].resize(0);
      	TPZGeoElSide mySide(geoel,s);
      	TPZGeoElSide neighS = mySide.Neighbour();
         if(mySide.Dimension() == 0){
         	long oldSz = NodesSequence.NElements();
            NodesSequence.resize(oldSz+1);
            NodesSequence[oldSz] = geoel->NodeIndex(s);
         }
      	while(mySide != neighS){
         	neighbourhood[s].push_back(neighS);
            neighS = neighS.Neighbour();
         }
      }

      //inserting in new element
      for(int s = 0; s < nsides; s++){
      	TPZGeoEl * tempEl = newgeoel;
         TPZGeoElSide tempSide(newgeoel,s);
         int byside = s;
         for(unsigned long n = 0; n < neighbourhood[s].size(); n++){
         	TPZGeoElSide neighS = neighbourhood[s][n];
            tempEl->SetNeighbour(byside, neighS);
            tempEl = neighS.Element();
            byside = neighS.Side();
         }
         tempEl->SetNeighbour(byside, tempSide);
      }

      long fatherindex = geoel->FatherIndex();
      if(fatherindex>-1) newgeoel->SetFather(fatherindex);

      if(!geoel->HasSubElement()) continue;

      int nsons = geoel->NSubElements();

      TPZAutoPointer<TPZRefPattern> ref = gRefDBase.GetUniformRefPattern(ETriangle);
      newgeoel->SetRefPattern(ref);

      for(int j=0;j<nsons;j++){
      	TPZGeoEl* son = geoel->SubElement(j);
         if(!son){
             DebugStop();
         }
         newgeoel->SetSubElement(j,son);
      }
   }

	/*Now, build connectivities*/
	newgmesh->BuildConnectivity();

	return newgmesh;
}
/*}}}*/
void AmrNeopz::CheckMesh(int** pdata,double** pxy,int** pelements){/*{{{*/

	/*Basic verification*/
	if(!pdata) _error_("Impossible to continue: pdata is NULL!\n");
	if(**pdata<=0) _error_("Impossible to continue: nvertices <=0!\n");
	if(*(*pdata+1)<=0) _error_("Impossible to continue: nelements <=0!\n");
	if(!pxy) _error_("Impossible to continue: pxy is NULL!\n");
	if(!pelements) _error_("Impossible to continue: pelements is NULL!\n");
}
/*}}}*/
void AmrNeopz::PrintGMeshVTK(TPZGeoMesh* gmesh,std::ofstream &file,bool matColor){/*{{{*/

	file.clear();
	long nelements = gmesh->NElements();
	TPZGeoEl *gel;
	std::stringstream node, connectivity, type, material;

	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";

	long actualNode = -1, size = 0, nVALIDelements = 0;
	for(long el = 0; el < nelements; el++){
	  gel = gmesh->ElementVec()[el];
	  if(!gel )//|| (gel->Type() == EOned && !gel->IsLinearMapping()))//Exclude Arc3D and Ellipse3D
	  {
			continue;
	  }
	  if(gel->HasSubElement())
	  {
			continue;
	  }
	  MElementType elt = gel->Type();
	  int elNnodes = MElementType_NNodes(elt);

	  size += (1+elNnodes);
	  connectivity << elNnodes;

	  for(int t = 0; t < elNnodes; t++)
	  {
			for(int c = 0; c < 3; c++)
			{
				 double coord = gmesh->NodeVec()[gel->NodeIndex(t)].Coord(c);
				 node << coord << " ";
			}
			node << std::endl;

			actualNode++;
			connectivity << " " << actualNode;
	  }
	  connectivity << std::endl;

	  int elType = this->GetVTK_ElType(gel);
	  type << elType << std::endl;

	  if(matColor == true)
	  {
			material << gel->MaterialId() << std::endl;
	  }
	  else
	  {
			material << gel->Index() << std::endl;
	  }

	  nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();

	file << "CELLS " << nVALIDelements << " ";

	file << size << std::endl;
	file << connectivity.str() << std::endl;

	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str() << std::endl;

	file << "CELL_DATA" << " " << nVALIDelements << std::endl;
	file << "FIELD FieldData 1" << std::endl;
	if(matColor == true)
	{
	  file << "material 1 " << nVALIDelements << " int" << std::endl;
	}
	else
	{
	  file << "ElementIndex 1 " << nVALIDelements << " int" << std::endl;
	}
	file << material.str();
	file.close();
}
/*}}}*/
int AmrNeopz::GetVTK_ElType(TPZGeoEl * gel){/*{{{*/

	MElementType pzElType = gel->Type();

    int elType = -1;
    switch (pzElType)
    {
        case(EPoint):
        {
            elType = 1;
            break;
        }
        case(EOned):
        {
            elType = 3;    
            break;
        }
        case (ETriangle):
        {
            elType = 5;
            break;                
        }
        case (EQuadrilateral):
        {
            elType = 9;
            break;                
        }
        case (ETetraedro):
        {
            elType = 10;
            break;                
        }
        case (EPiramide):
        {
            elType = 14;
            break;                
        }
        case (EPrisma):
        {
            elType = 13;
            break;                
        }
        case (ECube):
        {
            elType = 12;
            break;                
        }
        default:
        {
            std::cout << "Element type not found on " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
            break;    
        }
    }
    if(elType == -1)
    {
        std::cout << "Element type not found on " << __PRETTY_FUNCTION__ << std::endl;
        std::cout << "MIGHT BE CURVED ELEMENT (quadratic or quarter point)" << std::endl;
        DebugStop();
    }

    return elType;
}
/*}}}*/
void AmrNeopz::ReadMesh(){/*{{{*/

	std::string fathermeshfile		= "/home/santos/issm_fathermesh.txt";
	std::string previousmeshfile	= "/home/santos/issm_previousmesh.txt";
	std::string amrfile				= "/home/santos/issm_amr.txt";
	std::ifstream amrifstream(amrfile.c_str());				
	int size,value;

	TPZPersistenceManager::OpenRead(fathermeshfile);
	this->fathermesh = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::ReadFromFile());
	TPZPersistenceManager::CloseRead();

	TPZPersistenceManager::OpenRead(previousmeshfile);
	this->previousmesh = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::ReadFromFile());
   TPZPersistenceManager::CloseRead();

	if(!amrifstream.is_open()) _error_("amr ifstream is not open!");
	amrifstream.seekg(0);

	this->sid2index.clear();
	this->index2sid.clear();
	this->specialelementsindex.clear();

	amrifstream>>size;
	for(int i=0;i<size;i++){
		amrifstream>>value;
		this->sid2index.push_back(value);
	}

	amrifstream>>size;
	for(int i=0;i<size;i++){
		amrifstream>>value;
		this->index2sid.push_back(value);
	}

	amrifstream>>size;
	for(int i=0;i<size;i++){
		amrifstream>>value;
		this->specialelementsindex.push_back(value);
	}
}
/*}}}*/
void AmrNeopz::WriteMesh(){/*{{{*/

	std::string fathermeshfile		= "/home/santos/issm_fathermesh.txt";
	std::string previousmeshfile	= "/home/santos/issm_previousmesh.txt";
	std::string amrfile				= "/home/santos/issm_amr.txt";
	std::ofstream amrofstream(amrfile.c_str());				

	if(this->fathermesh){
		TPZPersistenceManager::OpenWrite(fathermeshfile);
		TPZPersistenceManager::WriteToFile(this->fathermesh);
		TPZPersistenceManager::CloseWrite();
	}

	if(this->previousmesh){
		TPZPersistenceManager::OpenWrite(previousmeshfile);
		TPZPersistenceManager::WriteToFile(this->previousmesh);
		TPZPersistenceManager::CloseWrite();
	}

	if(this->sid2index.size()>0){
		amrofstream << this->sid2index.size() << std::endl;
		for(int i=0;i<this->sid2index.size();i++) {
			amrofstream << this->sid2index[i] << std::endl;
		}
	}

	if(this->index2sid.size()>0){
		amrofstream << this->index2sid.size() << std::endl;
		for(int i=0;i<this->index2sid.size();i++) {
			amrofstream << this->index2sid[i] << std::endl;
		}
	}

	if(this->specialelementsindex.size()){
		amrofstream << this->specialelementsindex.size() << std::endl;
		for(int i=0;i<this->specialelementsindex.size();i++) {
			amrofstream << this->specialelementsindex[i] << std::endl;
		}
	}
	amrofstream.flush();
	amrofstream.close();
}
/*}}}*/
