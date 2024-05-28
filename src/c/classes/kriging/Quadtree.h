
#ifndef _QUADTREE_H
#define _QUADTREE_H

class Observation;

class Quadtree{

	private:
		/* A quadtree box contains up to 4 points (locations). 4 new quadtree boxes are
		 * created if a fifth point is added to the same box. A Quadtree box is therefore
		 * composed of EITHER:
		 * - up to 4 vertices
		 * - 4 "sub" quadtree boxes*/

		class QuadtreeBox: public Object{ 
			public:
				int    nbitems; // number of current vertices in the box
				double xcenter; // x position of the center (double)
				double ycenter; // x position of the center (double)
				double length;  // width of the box
				union{
					QuadtreeBox *box[4];
					Observation *obs[4];
				};

				/*Object functions (Needed because the Quadtree uses a Container*/
				Object *copy();
				void    DeepEcho()  {_error_("not implemented yet"); };
				void    Echo();
				int     Id()        {_error_("not implemented yet"); };
				void    Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!");};
				int     ObjectEnum(){_error_("not implemented yet"); };

				/*Methods*/
				int          IsWithinRange(double  x,double y,double range);
				void         RangeSearch(int *indices,int *pnobs,double x,double y,double range);
				void         WriteObservations(int *indices,int *pnobs);

		};

		/*Quadtree private Fields*/
		DataSet* boxcontainer;

	public:
		int          MaxDepth;          // maximum number of subdivision
		QuadtreeBox *root;              // main box
		int          NbQuadtreeBox;     // total number of boxes
		int          NbObs;             // number of points

		Quadtree();
		Quadtree(double xmin,double xmax,double ymin,double ymax,int maxdepth_in);
		~Quadtree();
		void         Add(Observation *observation);
		void         AddAndAverage(double x,double y,double value);
		void         ClosestObs(int *pindex,double x,double y);
		void         DeepEcho(void);
		void         Echo(void);
		void         IntergerCoordinates(int *xi,int *yi,double x,double y);
		QuadtreeBox *NewQuadtreeBox(double xcenter,double ycenter,double length);
		QuadtreeBox *NewQuadtreeBox(QuadtreeBox* master,int index);
		void         QuadtreeDepth(int *A,int xi,int yi);
		void         QuadtreeDepth2(int *A,int xi,int yi);
		void         RangeSearch(int **pindices,int *pnobs,double x,double y,double range);
};
#endif //_QUADTREE_H
