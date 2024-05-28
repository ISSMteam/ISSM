#ifndef AMRBAMG
#define AMRBAMG

/*Includes*/
#include "../shared/shared.h"
#include "../toolkits/toolkits.h"

#include "../bamg/BamgMesh.h"
#include "../bamg/BamgGeom.h"
#include "../bamg/BamgOpts.h"

class AmrBamg{

	public:
		int fieldenum;
		int keepmetric;
		IssmDouble groundingline_resolution;
		IssmDouble groundingline_distance;
		IssmDouble icefront_resolution;
		IssmDouble icefront_distance;
		IssmDouble thicknesserror_resolution;
		IssmDouble thicknesserror_threshold;
		IssmDouble thicknesserror_groupthreshold;
		IssmDouble thicknesserror_maximum;
		IssmDouble deviatoricerror_resolution;
		IssmDouble deviatoricerror_threshold;
		IssmDouble deviatoricerror_groupthreshold;
		IssmDouble deviatoricerror_maximum;

		/* Constructor, destructor etc*/
		AmrBamg();

		~AmrBamg();

		/*General methods*/
		void Initialize();
		void SetMesh(int** elementslist_in,IssmDouble** x_in,IssmDouble** y_in,int* numberofvertices,int* numberofelements);
		void GetMesh(int** elementslist_out,IssmDouble** x_out,IssmDouble** y_out,int* numberofvertices,int* numberofelements);
		void ExecuteRefinementBamg(IssmDouble* field,IssmDouble* hmaxVertices,int** pdatalist,IssmDouble** pxylist,int** pelementslist);
		void SetBamgOpts(IssmDouble hmin_in,IssmDouble hmax_in,IssmDouble err_in,IssmDouble gradation_in);

		/*Access Method*/
		BamgOpts* GetBamgOpts(){return this->options;}

	private:
		BamgGeom* geometry;
		BamgMesh* fathermesh;
		BamgMesh* previousmesh;
		BamgOpts* options;
		/*entire mesh*/
		IssmDouble* x;
		IssmDouble* y;
		int* elementslist;
		int numberofvertices;
		int numberofelements;
};

#endif
