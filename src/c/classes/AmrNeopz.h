#ifndef AMRNEOPZ
#define AMRNEOPZ

/*Includes*/
/*{{{*/
/*NeoPZ includes*/
#include <pzgmesh.h>
/*ISSM includes*/
#include "../shared/shared.h"
/*}}}*/

class AmrNeopz{

public:
	/*Public attributes{{{*/
	/* 
	 * to refine:
	 * distance_h = initial_distance * gradation ^ (level_max-h)
	 * to unrefine:
	 * distance_h = lag * initial_distance * gradation ^ (level_max-h)
	 */
	int refinement_type;						//0 uniform (faster); 1 refpattern  
	int level_max;								//max level of refinement
	double gradation;							//geometric progression ratio to calculate radius of level h
	double lag;									//lag used in the unrefine process
	/*Target and estimators*/
	double groundingline_distance;		//all elements with distance from grounding line <= groundingline_distance will be refined
	double icefront_distance;				//all elements with distance from ice front <= icefront_distance will be refined
	double thicknesserror_threshold;		//if ==0, it will not be used
	double thicknesserror_groupthreshold;//group threshold
	double thicknesserror_maximum;		//max value of the error estimator; in general, it is defined in the first time step. Attention with restart
	double deviatoricerror_threshold;	//if ==0, it will not be used
	double deviatoricerror_groupthreshold;//group threshold
	double deviatoricerror_maximum;		//max value of the error estimator; in general, it is defined in the first time step. Attention with restart
	/*}}}*/
	/*Public methods{{{*/
	/* Constructor, destructor etc*/
	AmrNeopz();															
	AmrNeopz(const AmrNeopz &cp); 					
	AmrNeopz & operator= (const AmrNeopz &cp);	
	virtual ~AmrNeopz();												
	/*General methods*/
	void CleanUp();
	void Initialize();
   void ExecuteRefinement(double* gl_distance,double* if_distance,double* deviatoricerror,double* thicknesserror,int** pdatalist,double** pxy,int** pelementslist);
	void SetMesh(int** elementslist_in,IssmDouble** x_in,IssmDouble** y_in,int* numberofvertices,int* numberofelements);
	void GetMesh(int** elementslist_out,IssmDouble** x_out,IssmDouble** y_out,int* numberofvertices,int* numberofelements);
	void CheckMesh(int** pdata,double** pxy,int** pelements);
	void ReadMesh();
	void WriteMesh();
	/*}}}*/
private:
	/*Private attributes{{{*/
	std::vector<int> sid2index;					// Vector that keeps index of PZGeoMesh elements used in the ISSM mesh (sid) 
	std::vector<int> index2sid;					// Vector that keeps sid of issm mesh elements used in the neopz mesh (index) 
	std::vector<int> specialelementsindex;		// Vector that keeps index of the special elements (created to avoid haning nodes) 
	TPZGeoMesh *fathermesh;							// Entire mesh without refinement if refinement_type==1; refined with hanging nodes if efinement_type==0
	TPZGeoMesh *previousmesh;						// Refined mesh without hanging nodes (it is always refpattern type), used to generate ISSM mesh
	IssmDouble* x; 									// Entire mesh
   IssmDouble* y;
   int* elementslist;
   int numberofvertices;
   int numberofelements;	
	/*}}}*/
	/*Private methods{{{*/
   void RefineMeshOneLevel(bool &verbose,double* gl_distance,double* if_distance,double* deviatoricerror,double* thicknesserror);
	void RefineMeshWithSmoothing(bool &verbose,TPZGeoMesh* gmesh);
	void RefineMeshToAvoidHangingNodes(bool &verbose,TPZGeoMesh* gmesh);
	void DeleteSpecialElements(bool &verbose,TPZGeoMesh* gmesh);
	void GetMesh(TPZGeoMesh* gmesh,int** pdata,double** pxy,int** pelements);
	TPZGeoMesh* CreateRefPatternMesh(TPZGeoMesh* gmesh);
   inline int GetElemMaterialID(){return 1;} 
	inline int GetNumberOfNodes(){return 3;}
	void PrintGMeshVTK(TPZGeoMesh *gmesh,std::ofstream &file,bool matColor=true);
	int GetVTK_ElType(TPZGeoEl* gel);
	int VerifyRefinementType(TPZGeoEl* geoel,TPZGeoMesh* gmesh);
	/*}}}*/
};

#endif
