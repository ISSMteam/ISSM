#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include "./bamgobjects.h"

namespace bamg {

#define ABS(i) ((i)<0 ?-(i) :(i))
#define MAXDEPTH  30
#define MAXISIZE  1073741824 //2^30
#define MAXICOORD 1073741823 //2^30 - 1 = =111...111 (29 times one)

	/*DOCUMENTATION What is a BamgQuadtree? {{{
	 * A Quadtree is a very simple way to group vertices according
	 * to their locations. A square that holds all the points of the mesh
	 * (or the geometry) is divided into 4 boxes. As soon as one box
	 * hold more than 4 vertices, it is divided into 4 new boxes, etc...
	 * There cannot be more than MAXDEEP (=30) subdivision.
	 * This process is like a Dichotomy in dimension 2
	 *
	 *  + - -  -    - -    -    - - + -   - + - + - + - -     - - +
	 *  |                           |       |   | X |             |
	 *                                      + - + - +
	 *  |                           |       |   |   |             |
	 *                              + -   - + - + - +             +
	 *  |                           |       |       |             |
	 *                         
	 *  |                           |       |       |             |
	 *  + - -  -    - -    -    - - + -   - + -   - + - -     - - +
	 *  |                           |               |             |
	 *                         
	 *  |                           |               |             |
	 *                         
	 *  |                           |               |             |
	 *  |                           |               |             |
	 *  + - -  -    - -    -    - - + -   -   -   - + - -     - - +
	 *  |                           |                             |
	 *                         
	 *  |                           |                             |
	 *                         
	 *  |                           |                             |
	 *                         
	 *  |                           |                             |
	 *  |                           |                             |
	 *  |                           |                             |
	 *  |                           |                             |
	 *  |                           |                             |
	 *  + - -  -    - -    -    - - + -   -   -   -   - -     - - +
	 *
	 * The coordinate system used in a quadtree are integers to avoid
	 * round-off errors. The vertex in the lower left box has the coordinates
	 * (0 0) 
	 * The upper right vertex has the follwing coordinates:
	 * 2^30 -1           2^30 -1        in decimal
	 * 0 1 1 1 .... 1    0 1 1 1 .... 1 in binary
	 *  \--   29  --/     \--   29  --/
	 * Using binaries is therefore very easy to locate a vertex in a box:
	 * we just need to look at the bits from the left to the right (See ::Add)
	 }}}*/

	/*Constructors/Destructors*/
	BamgQuadtree::BamgQuadtree(){/*{{{*/

		/*Initialize fields*/
		this->NbQuadtreeBox = 0;
		this->NbVertices    = 0;

		/*Create Root, pointer toward the main box*/
		this->root=NewBamgQuadtreeBox();

	} /*}}}*/
	BamgQuadtree::BamgQuadtree(Mesh * t,long nbv){ /*{{{*/

		/*Initialize fields*/
		this->NbQuadtreeBox = 0;
		this->NbVertices    = 0;

		/*Create Root, pointer toward the main box*/
		this->root=NewBamgQuadtreeBox();

		/*Check Sizes*/
		_assert_(MAXISIZE>MAXICOORD);

		/*Add all vertices of the mesh*/
		if(nbv==-1) nbv=t->nbv;
		for(int i=0;i<nbv;i++){
			this->Add(t->vertices[i]);
		}

	}
	/*}}}*/
	BamgQuadtree::~BamgQuadtree() {/*{{{*/

		vector<BamgQuadtreeBox*>::reverse_iterator object;
		for(object=boxcontainer.rbegin() ; object <boxcontainer.rend(); object++ ){
			delete (*object);
		}
		boxcontainer.clear();
		root=NULL;
	}
	/*}}}*/

	/*Methods*/
	void          BamgQuadtree::Add(BamgVertex &w){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, BamgQuadtree.cpp/Add)*/
		BamgQuadtreeBox** pb=NULL;
		BamgQuadtreeBox*  b=NULL;

		/*Get integer coodinate of current point w*/
		long i=w.i.x, j=w.i.y;

		/*Initialize level*/
		long level=MAXISIZE;

		/*Get inital box (the largest)*/
		pb = &root;

		/*Find the smallest box where w is located*/
		while((b=*pb) && (b->nbitems<0)){ 

			//shift b->nbitems by -1
			b->nbitems--;

			//shifted righ by one bit: level=00000010 -> 00000001
			level >>= 1;

			//Get next subbox according to the bit value (level)
			pb = &b->box[BoxNumber(i,j,level)];
		}

		/*OK, we have found b, a Subbox holding vertices (might be full)
		  check that the vertex is not already in the box*/
		if (b){      
			if (b->nbitems > 3 &&  b->v[3] == &w) return;
			if (b->nbitems > 2 &&  b->v[2] == &w) return;
			if (b->nbitems > 1 &&  b->v[1] == &w) return;
			if (b->nbitems > 0 &&  b->v[0] == &w) return;
		}

		/*check that l is not 0 (this should not happen as MAXDEPTH = 30)*/
		_assert_(level>0);

		/*Now, try to add the vertex, if the subbox is full (nbitems=4), we have to divide it
		  in 4 new subboxes*/
		while ((b= *pb) && (b->nbitems == 4)){ // the BamgQuadtreeBox is full

			/*Copy the 4 vertices in the current BamgQuadtreebox*/
			BamgVertex* v4[4];
			v4[0]= b->v[0];
			v4[1]= b->v[1];
			v4[2]= b->v[2];
			v4[3]= b->v[3];

			/*set nbitems as negative 
			 * (box full -> holds 4 pointers toward subboxes and not 4 vertices)*/
			b->nbitems = -b->nbitems;

			/*Initialize the 4 pointers toward the 4 subboxes*/
			b->box[0]=b->box[1]=b->box[2]=b->box[3]=NULL;

			/*level = 0010000 -> 0001000*/
			level >>= 1;

			/*Put the four vertices in the new boxes*/
			for (int k=0;k<4;k++){

				int          ij;
				/*bb is the new "sub"box of b where v4[k] is located*/
				BamgQuadtreeBox *bb = b->box[ij=BoxNumber(v4[k]->i.x,v4[k]->i.y,level)];

				// alloc the BamgQuadtreeBox
				if (!bb) bb=b->box[ij]=NewBamgQuadtreeBox(); 

				/*Copy the current vertex*/
				bb->v[bb->nbitems++] = v4[k];
			}

			/*Get the subbox where w (i,j) is located*/
			pb = &b->box[BoxNumber(i,j,level)];
		}

		/*alloc the BamgQuadtreeBox if necessary*/
		if (!(b=*pb)) b=*pb= NewBamgQuadtreeBox();

		/*Add w*/
		b->v[b->nbitems++]=&w;

		//Increase NbVertices by one (we have one new vertex)
		NbVertices++;
	}
	/*}}}*/
	int           BamgQuadtree::BoxNumber(int i,int j,int l){/*{{{*/
		/* 
		 * 
		 *    J    j
		 *    ^    ^
		 *    |    | +--------+--------+
		 *    |    | |        |        |
		 * 1X |    | |   2    |   3    |
		 *    |    | |        |        |
		 *    |    | +--------+--------+
		 *    |    | |        |        |
		 * 0X |    | |   0    |   1    |
		 *    |    | |        |        |
		 *    |    | +--------+--------+
		 *    |    +-----------------------> i
		 *    |         
		 *    |----------------------------> I
		 *              X0        X1  
		 *
		 * box 0 -> I=0 J=0 BoxNumber=00  = 0
		 * box 1 -> I=1 J=0 BoxNumber=01  = 1
		 * box 2 -> I=0 J=1 BoxNumber=10  = 2
		 * box 3 -> I=1 J=1 BoxNumber=11  = 3
		 */
		//BoxNumber(i,j,l) returns the box number of i and j with respect to l
		//if !j&l and !i&l -> 0 (box zero: lower left )
		//if !j&l and  i&l -> 1 (box one:  lower right)
		//if  j&l and !i&l -> 2 (box two:  upper left )
		//if  j&l and  i&l -> 3 (box three:upper right)
		return ((j&l) ? ((i&l) ? 3:2 ) :((i&l) ? 1:0 ));
	}/*}}}*/
	bool          BamgQuadtree::Intersection(int a,int b,int x,int y){/*{{{*/
		/*is [x y] in [a b]*/
		return ((y) > (a)) && ((x) <(b));
	}/*}}}*/
	BamgVertex*   BamgQuadtree::NearestVertex(int xi,int yi) {/*{{{*/

		/*initial output as NULL (no vertex found)*/
		BamgVertex*  nearest_v=NULL;

		/*if the tree is empty, return NULL pointer*/
		if(!this->root->nbitems) return nearest_v; 

		/*Project coordinates (xi,yi) onto [0,MAXICOORD-1] x [0,MAXICOORD-1]*/
		int xi2 = xi;
		int yi2 = yi;
		if(xi<0)        xi2 = 0;
		if(xi>MAXISIZE) xi2 = MAXICOORD;
		if(yi<0)        yi2 = 0;
		if(yi>MAXISIZE) yi2 = MAXICOORD;

		/*Get inital box (the largest)*/
		BamgQuadtreeBox* b = this->root;

		/*Initialize level and sizes for largest box*/
		int levelbin = (1L<<MAXDEPTH);// = 2^30
		int        h = MAXISIZE;
		int       hb = MAXISIZE;
		int       i0 = 0;
		int       j0 = 0;

		/*else, find the smallest non-empty BamgQuadtreeBox containing  the point (i,j)*/
		while(b->nbitems<0){

			int hb2 = hb >> 1;                  //size of the current box
			int k   = BoxNumber(xi2,yi2,hb2);   //box number (0,1,2 or 3)
			BamgQuadtreeBox* b0  = b->box[k];   //pointer toward current box

			/* break if box is empty (Keep previous box b)*/
			if((b0==NULL) || (b0->nbitems==0)) break;

			/*Get next Quadtree box*/
			b  = b0;
			hb = hb2;
			this->SubBoxCoords(&i0,&j0,k,hb2);
		}

		/*The box b, is the smallest box containing the point (i,j) and
		 * has the following properties:
		 * - n0: number of items (>0 if vertices, else boxes)
		 * - hb: box size (int)
		 * - i0: x coordinate of the lower left corner
		 * - j0: y coordinate of the lower left corner*/
		int n0 = b->nbitems;

		/* if the current subbox is holding vertices, we are almost done*/
		if(n0>0){  
			/*loop over the vertices of the box and find the closest vertex*/
			for(int k=0;k<n0;k++){
				int xiv = b->v[k]->i.x;
				int yiv = b->v[k]->i.y;

				int h0 = Norm(xi2,xiv,yi2,yiv);

				/*is it smaller than previous value*/
				if(h0<h){
					h = h0;
					nearest_v = b->v[k];
				}
			}
			/*return closest vertex*/
			return nearest_v;
		}

		/* general case: the current box is empty, we have to go backwards
			and find the closest not-empty box and find the closest vertex*/

		/*Initialize search variables*/
		BamgQuadtreeBox* pb[MAXDEPTH];
		int pi[MAXDEPTH];
		int ii[MAXDEPTH];
		int jj[MAXDEPTH];
		pb[0]=b;                             //pointer toward the box b
		pi[0]=b->nbitems>0?(int)b->nbitems:4;//number of boxes in b
		ii[0]=i0;                            //i coordinate of the box lowest left corner
		jj[0]=j0;                            //j coordinate of the box lowest left corner

		/*initialize h: smallest box size, containing a vertex close to w*/
		h=hb;

		/*Main loop*/
		int level=0;
		do{
			/*get current box*/
			b = pb[level];

			/*Loop over the items in current box (if not empty!)*/
			while(pi[level]){

				/*We are looping now over the items of b. k is the current index (in [0 3])*/
				pi[level]--;
				int k=pi[level];

				/*if the current subbox is holding vertices (b->nbitems<0 is subboxes)*/
				if (b->nbitems>0){
					int h0 = Norm(xi2,b->v[k]->i.x,yi2,b->v[k]->i.y);
					if (h0<h){
						h=h0;
						nearest_v=b->v[k];
					}
				}
				/*else: current box b is pointing toward 4 boxes
				 * test sub-box k and go deeper into the tree if it is non empty
				 * and contains the point w modulo a size h that is either the size of the smallest
				 * non empty box containing w, or the closest point to w (so far) */
				else{
					BamgQuadtreeBox* b0=b;

					/*if the next box exists:*/
					if((b=b->box[k])){

						/*Get size (hb) and coordinates of the current sub-box lowest left corner*/
						hb>>=1;
						int iii = ii[level];
						int jjj = jj[level];
						this->SubBoxCoords(&iii,&jjj,k,hb);

						/*if the current point (xi2,yi2) is in b (modulo h), this box is good:
						 * it is holding vertices that are close to w */
						if(this->Intersection(iii,iii+hb,xi2-h,xi2+h) && this->Intersection(jjj,jjj+hb,yi2-h,yi2+h)){
							level++;
							pb[level] = b;
							pi[level] = b->nbitems>0 ? b->nbitems:4  ;
							ii[level] = iii;
							jj[level] = jjj;
						}

						/*else go backwards*/
						else{
							/*shifted righ by one bit: hb=001000000 -> 01000000*/
							b=b0;
							hb<<=1;
						}
					}
					else{
						/*Current box is NULL, go to next subbox of b (k=k-1)*/
						b=b0;
					}
				}
			}

			/*We have found a vertex, now, let's try the other boxes of the previous level
			 * in case there is a vertex closest to w that has not yet been tested*/
			hb <<= 1;
		}while(level--);

		/*return nearest vertex pointer*/
		return nearest_v;
	}
	/*}}}*/
	int           BamgQuadtree::Norm(int xi1,int xi2,int yi1,int yi2){/*{{{*/

		int deltax = xi2 - xi1;
		int deltay = yi2 - yi1;

		if(deltax<0) deltax = -deltax;
		if(deltay<0) deltay = -deltay;

		if(deltax> deltay){
			return deltax;
		}
		else{
			return deltay;
		}
	}/*}}}*/
	void          BamgQuadtree::SubBoxCoords(int* pi,int*pj,int boxnumber,int length){/*{{{*/
		/* 
		 *         j (first bit)
		 *         ^
		 *         | +--------+--------+
		 *         | |        |        |
		 *   1X    | |   2    |   3    |
		 *         | |        |        |
		 *         | +--------+--------+
		 *         | |        |        |
		 *   0X    | |   0    |   1    |
		 *         | |        |        |
		 *         | +--------+--------+
		 *         +-----------------------> i (second bit)
		 *               X0       X1
		 */

		/*Add sub-box coordinate to i and j*/
		//_assert_(!(*pi & length));
		//_assert_(!(*pj & length));

		/*length if first  bit of boxnumber is 1, elengthse 0*/
		*pi += ((boxnumber & 1) ? length:0);
		/*length if second bit of boxnumber is 1, elengthse 0*/
		*pj += ((boxnumber & 2) ? length:0);

	}/*}}}*/
	BamgVertex*   BamgQuadtree::TooClose(BamgVertex* v,double threshold,int hx,int hy){/*{{{*/

		const int i=v->i.x;
		const int j=v->i.y;
		const double Xx = v->r.x;
		const double Xy = v->r.y;
		Metric* Mx = new Metric(v->m);

		BamgQuadtreeBox *pb[MAXDEPTH];
		int  pi[MAXDEPTH];
		int  ii[MAXDEPTH], jj[MAXDEPTH];
		int l=0; // level
		BamgQuadtreeBox* b;
		int hb =  MAXISIZE;
		int i0=0,j0=0;

		// BamgVertex *vn=0;

		if(!root->nbitems) return 0; // empty tree 

		// general case
		pb[0]=root;
		pi[0]=root->nbitems>0 ?(int)root->nbitems:4;
		ii[0]=i0;
		jj[0]=j0;
		do{
			b= pb[l];
			while(pi[l]--){ 	      
				int k = pi[l];

				if(b->nbitems>0){ // BamgVertex BamgQuadtreeBox none empty
					int i2x = b->v[k]->i.x;
					int i2y = b->v[k]->i.y;
					if (ABS(i-i2x)<hx && ABS(j-i2y) <hy ){
						double XYx = b->v[k]->r.x - Xx;
						double XYy = b->v[k]->r.y - Xy;
						if(LengthInterpole(Mx->Length(XYx,XYy),b->v[k]->m.Length(XYx,XYy)) < threshold){
							delete Mx;
							return b->v[k]; 
						}
					}
				}
				else{ // Pointer BamgQuadtreeBox 
					BamgQuadtreeBox *b0=b;
					if ((b=b->box[k])){
						hb >>=1 ; // div by 2
						int iii = ii[l];
						int jjj = jj[l];
						this->SubBoxCoords(&iii,&jjj,k,hb);

						if(this->Intersection(iii,iii+hb,i-hx,i+hx) && this->Intersection(jjj,jjj+hb,j-hy,j+hy)){
							pb[++l]=  b;
							pi[l]= b->nbitems>0 ?(int)  b->nbitems : 4  ;
							ii[l]= iii;
							jj[l]= jjj;

						}
						else{
							b=b0;
							hb <<=1 ;
						}
					}
					else{
						b=b0;
					}
				}
			}
			hb <<= 1; // mul by 2 
		}while(l--);

		delete Mx;
		return 0;
	}
	/*}}}*/

	BamgQuadtree::BamgQuadtreeBox* BamgQuadtree::NewBamgQuadtreeBox(void){/*{{{*/

		/*Output*/
		BamgQuadtreeBox* newbox=NULL;

		/*Create and initialize a new box*/
		newbox=new BamgQuadtreeBox;
		newbox->nbitems=0;
		newbox->box[0]=NULL;
		newbox->box[1]=NULL;
		newbox->box[2]=NULL;
		newbox->box[3]=NULL;

		/*Add root to the container*/
		boxcontainer.push_back(newbox);
		this->NbQuadtreeBox++;

		/*currentbox now points toward next quadtree box*/
		return newbox;
	}/*}}}*/
}
