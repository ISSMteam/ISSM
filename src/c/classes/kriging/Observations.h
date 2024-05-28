#ifndef _CONTAINER_OBSERVATIONS_H_
#define  _CONTAINER_OBSERVATIONS_H_

class Quadtree;
class Covertree;
class Variogram;
class Options;
#include "../../datastructures/datastructures.h"
#include "../../shared/shared.h"

/*!\brief Declaration of Observations class.
 *
 * Declaration of Observations class.  Observations are vector lists (Containers) of Observation objects.
 */ 

class Observations: public DataSet{

	private:
		int        treetype;
		Quadtree*  quadtree;
		Covertree* covertree;

	public:

		/*constructors, destructors*/
		Observations();
		Observations(IssmPDouble* observations_list,IssmPDouble* x,IssmPDouble* y,int n,Options* options);
		~Observations();

		/*Initialize data structures*/
		void InitCovertree(IssmPDouble* observations_list,IssmPDouble* x,IssmPDouble* y,int n,Options* options);
		void InitQuadtree(IssmPDouble* observations_list,IssmPDouble* x,IssmPDouble* y,int n,Options* options);

		/*Methods*/
		void ClosestObservation(IssmPDouble *px,IssmPDouble *py,IssmPDouble *pobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius);
		void ClosestObservationCovertree(IssmPDouble *px,IssmPDouble *py,IssmPDouble *pobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius);
		void ClosestObservationQuadtree(IssmPDouble *px,IssmPDouble *py,IssmPDouble *pobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius);
		void Distances(IssmPDouble* distances,IssmPDouble *x,IssmPDouble *y,int n,IssmPDouble radius);
		void InterpolationIDW(IssmPDouble *pprediction,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int mindata,int maxdata,IssmPDouble power);
		void InterpolationV4(IssmPDouble *pprediction,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int mindata,int maxdata);
		void InterpolationKriging(IssmPDouble *pprediction,IssmPDouble *perror,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int mindata,int maxdata,Variogram* variogram);
		void InterpolationNearestNeighbor(IssmPDouble *pprediction,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius);
		void ObservationList(IssmPDouble **px,IssmPDouble **py,IssmPDouble **pobs,int* pnobs);
		void ObservationList(IssmPDouble **px,IssmPDouble **py,IssmPDouble **pobs,int* pnobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int maxdata);
		void ObservationListCovertree(IssmPDouble **px,IssmPDouble **py,IssmPDouble **pobs,int* pnobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int maxdata);
		void ObservationListQuadtree(IssmPDouble **px,IssmPDouble **py,IssmPDouble **pobs,int* pnobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int maxdata);
		void QuadtreeColoring(IssmPDouble* A,IssmPDouble *x,IssmPDouble *y,int n);
		void Variomap(IssmPDouble* gamma,IssmPDouble *x,int n);

};
#endif //ifndef _OBSERVATIONS_H_
