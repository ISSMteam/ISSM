#include "../classes.h"

/*DOCUMENTATION What is a Quadtree? {{{
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
/*MACROS {{{*/
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
 * box 0 -> I=0 J=0 IJ=00  = 0
 * box 1 -> I=1 J=0 IJ=01  = 1
 * box 2 -> I=0 J=1 IJ=10  = 2
 * box 3 -> I=1 J=1 IJ=11  = 3
 */
//IJ(i,j,l) returns the box number of i and j with respect to l
//if !j&l and !i&l -> 0 (box zero: lower left )
//if !j&l and  i&l -> 1 (box one:  lower right)
//if  j&l and !i&l -> 2 (box two:  upper left )
//if  j&l and  i&l -> 3 (box three:upper right)
#define IJ(i,j,l)  ((j&l) ? ((i&l) ? 3:2 ) :((i&l) ? 1:0 ))
/*}}}*/

	/*Constructors/Destructors*/
Quadtree::Quadtree(){/*{{{*/
	_error_("Constructor not supported");

}
/*}}}*/
Quadtree::Quadtree(double xmin,double xmax,double ymin,double ymax,int maxdepth){/*{{{*/

	/*Intermediaries*/
	double length;

	/*Initialize fields*/
	this->MaxDepth=maxdepth;
	this->NbQuadtreeBox=0;
	this->NbObs=0;

	/*Create container*/
	this->boxcontainer=new DataSet();

	/*Create Root, pointer toward the main box*/
	length=max(xmax-xmin,ymax-ymin);
	this->root=NewQuadtreeBox(xmin+length/2,ymin+length/2,length);
}
/*}}}*/
	Quadtree::~Quadtree(){/*{{{*/

		delete boxcontainer;
		root=NULL;

	}
	/*}}}*/

	/*Methods*/
void  Quadtree::Add(Observation* observation){/*{{{*/

	/*Intermediaries*/
	int          xi,yi,ij,level,levelbin;
	QuadtreeBox **pbox    = NULL; // pointer toward current box b
	QuadtreeBox **pmaster = NULL; // pointer toward master of b
	QuadtreeBox  *box     = NULL; // current box b
	QuadtreeBox  *slave   = NULL; // suslaveox of b (if necessary)
	Observation  *obs[4];

	/*Get integer coodinates*/
	xi = observation->xi;
	yi = observation->yi;

	/*Initialize levels*/
	level    = 0;
	levelbin = (1L<<this->MaxDepth);// = 2^30

	/*Get inital box (the largest)*/
	pmaster = &root;
	pbox    = &root;

	/*Find the smallest box where the observation is located*/
	while((box=*pbox) && (box->nbitems<0)){ 

		/*Go down one level (levelbin = 00100 -> 00010)*/
		levelbin>>=1; level+=1; _assert_(level<this->MaxDepth);

		/*Get next box according to the bit value (levelbin)*/
		pmaster = pbox;
		pbox    = &box->box[IJ(xi,yi,levelbin)];
	}
	_assert_(levelbin>0);

	/*Now, try to add the vertex, if the box is full (nbitems=4), we have to divide it in 4 new boxes*/
	while((box=*pbox) && (box->nbitems==4)){

		/*Copy the 4 observation in the current Quadtreebox*/
		obs[0] = box->obs[0];
		obs[1] = box->obs[1];
		obs[2] = box->obs[2];
		obs[3] = box->obs[3];

		/*set nbitems as -1 (now holding boxes instead of observations)*/
		box->nbitems = -1;
		box->box[0]  = NULL;
		box->box[1]  = NULL;
		box->box[2]  = NULL;
		box->box[3]  = NULL;

		/*Go down one level (levelbin = 00010 -> 00001)*/
		levelbin>>=1; level+=1; _assert_(level<this->MaxDepth);

		/*Put the four observations in the new boxes*/
		for (int k=0;k<4;k++){

			/*Get box for observation number k*/
			ij    = IJ(obs[k]->xi,obs[k]->yi,levelbin);
			slave = box->box[ij];
			if(!slave){
				box->box[ij] = NewQuadtreeBox(box,ij);
				slave        = box->box[ij];
			}
			slave->obs[slave->nbitems++] = obs[k];
		}

		/*Get the suslaveox where the current observation is located*/
		ij      = IJ(xi,yi,levelbin);
		pmaster = pbox;
		pbox    = &box->box[ij];
	}

	/*alloc the QuadtreeBox if necessary and add current observation*/
	box = *pbox;
	if(!box){
		ij  = IJ(xi,yi,levelbin);
		box = *pbox = NewQuadtreeBox(*pmaster,ij);
	}
	box->obs[box->nbitems++]=observation;
	NbObs++;

}/*}}}*/
void Quadtree::AddAndAverage(double x,double y,double value){/*{{{*/

	QuadtreeBox **pbox = NULL;
	QuadtreeBox  *box  = NULL;
	int           xi,yi;
	int           levelbin;
	int           index;
	double        length,length2;

	/*Get integer coodinates*/
	this->IntergerCoordinates(&xi,&yi,x,y);

	/*Initialize level*/
	levelbin = (1L<<this->MaxDepth);// = 2^30

	/*Get inital box (the largest)*/
	pbox=&root;

	/*Find the smallest box where this point is located*/
	while((box=*pbox) && (box->nbitems<0)){ 
		levelbin>>=1;
		pbox = &box->box[IJ(xi,yi,levelbin)];
	}

	/*Add obervation in this box (should be full)*/
	if(box && box->nbitems==4){
		index  = 0;
		length = pow(box->obs[0]->x - x,2) + pow(box->obs[0]->y - y,2);
		for(int i=1;i<4;i++){
			length2 = pow(box->obs[i]->x - x,2) + pow(box->obs[i]->y - y,2);
			if(length2<length){
				index  = i;
				length = length2;
			}
		}

		/*We found the closest observation, now average observation (do not change xi and yi to avoid round off errors*/
		box->obs[index]->x = (box->obs[index]->weight*box->obs[index]->x + x)/(box->obs[index]->weight+1.);
		box->obs[index]->y = (box->obs[index]->weight*box->obs[index]->y + y)/(box->obs[index]->weight+1.);
		box->obs[index]->xi= int((box->obs[index]->weight*double(box->obs[index]->xi) + double(xi))/(box->obs[index]->weight+1.));
		box->obs[index]->yi= int((box->obs[index]->weight*double(box->obs[index]->yi) + double(yi))/(box->obs[index]->weight+1.));
		box->obs[index]->value   = (box->obs[index]->weight*box->obs[index]->value + value)/(box->obs[index]->weight+1.);
		box->obs[index]->weight += 1.;
	}
	else{
		_error_("Box is not full");
	}
}/*}}}*/
void Quadtree::ClosestObs(int *pindex,double x,double y){/*{{{*/

	QuadtreeBox **pbox = NULL;
	QuadtreeBox  *box  = NULL;
	int           xi,yi;
	int           levelbin;
	int           index = -1;
	double        length,length2;

	/*Get integer coodinates*/
	this->IntergerCoordinates(&xi,&yi,x,y);

	/*Initialize level*/
	levelbin = (1L<<this->MaxDepth);// = 2^30

	/*Get inital box (the largest)*/
	pbox=&root;

	/*Find the smallest box where this point is located*/
	while((box=*pbox) && (box->nbitems<0)){ 
		levelbin>>=1;
		pbox = &box->box[IJ(xi,yi,levelbin)];
	}

	/*Add obervation in this box (should be full)*/
	if(box && box->nbitems>0){
		index  = box->obs[0]->index;
		length = pow(box->obs[0]->x - x,2) + pow(box->obs[0]->y - y,2);
		for(int i=1;i<box->nbitems;i++){
			length2 = pow(box->obs[i]->x - x,2) + pow(box->obs[i]->y - y,2);
			if(length2<length){
				index  = box->obs[i]->index;
				length = length2;
			}
		}
	}

	*pindex=index;
}/*}}}*/
void  Quadtree::DeepEcho(void){/*{{{*/

	_printf_("Quadtree:\n");
	_printf_("   MaxDepth      = " << this->MaxDepth << "\n");
	_printf_("   NbQuadtreeBox = " << this->NbQuadtreeBox << "\n");
	_printf_("   NbObs         = " << this->NbObs << "\n");
	_printf_("   root          = " << this->root << "\n");
	boxcontainer->Echo();

}/*}}}*/
void  Quadtree::Echo(void){/*{{{*/

	_printf_("Quadtree:\n");
	_printf_("   MaxDepth      = " << this->MaxDepth << "\n");
	_printf_("   NbQuadtreeBox = " << this->NbQuadtreeBox << "\n");
	_printf_("   NbObs         = " << this->NbObs << "\n");
	_printf_("   root          = " << this->root << "\n");

}/*}}}*/
void  Quadtree::IntergerCoordinates(int *xi,int *yi,double x,double y){/*{{{*/

	/*Intermediaries*/
	double coefficient;
	double xmin,ymin;

	/*Checks in debugging mode*/
	_assert_(xi && yi);
	_assert_(this->root);

	/*coeffIcoor is the coefficient used for integer coordinates:
	 *                (x-xmin)
	 * xi = (2^30 -1) --------- 
	 *                 length
	 * coefficient = (2^30 -1)/length
	 */
	coefficient = double((1L<<this->MaxDepth)-1)/(this->root->length);
	xmin        = this->root->xcenter - this->root->length/2;
	ymin        = this->root->ycenter - this->root->length/2;

	*xi=int(coefficient*(x - xmin));
	*yi=int(coefficient*(y - ymin));
}/*}}}*/
Quadtree::QuadtreeBox* Quadtree::NewQuadtreeBox(double xcenter,double ycenter,double length){/*{{{*/

	/*Output*/
	QuadtreeBox* newbox=NULL;

	/*Create and initialize a new box*/
	newbox=new QuadtreeBox();
	newbox->nbitems=0;
	newbox->xcenter=xcenter;
	newbox->ycenter=ycenter;
	newbox->length=length;
	newbox->box[0]=NULL;
	newbox->box[1]=NULL;
	newbox->box[2]=NULL;
	newbox->box[3]=NULL;

	/*Add to container*/
	this->boxcontainer->AddObject(newbox);
	NbQuadtreeBox++;

	/*currentbox now points toward next quadtree box*/
	return newbox;
}/*}}}*/
Quadtree::QuadtreeBox* Quadtree::NewQuadtreeBox(QuadtreeBox* master,int index){/*{{{*/

	/*Output*/
	QuadtreeBox* newbox=NULL;

	/*Checks in debugging mode*/
	_assert_(master);

	/*Create and initialize a new box*/
	newbox=new QuadtreeBox();
	newbox->nbitems=0;
	newbox->box[0]=NULL;
	newbox->box[1]=NULL;
	newbox->box[2]=NULL;
	newbox->box[3]=NULL;
	switch(index){
		case 0:
			newbox->xcenter=master->xcenter - master->length/4;
			newbox->ycenter=master->ycenter - master->length/4;
			break;
		case 1:
			newbox->xcenter=master->xcenter + master->length/4;
			newbox->ycenter=master->ycenter - master->length/4;
			break;
		case 2:
			newbox->xcenter=master->xcenter - master->length/4;
			newbox->ycenter=master->ycenter + master->length/4;
			break;
		case 3:
			newbox->xcenter=master->xcenter + master->length/4;
			newbox->ycenter=master->ycenter + master->length/4;
			break;
		default:
			_error_("Case " << index << " not supported");
	}
	newbox->length=master->length/2;

	/*Add to container*/
	this->boxcontainer->AddObject(newbox);
	NbQuadtreeBox++;

	/*currentbox now points toward next quadtree box*/
	return newbox;
}/*}}}*/
void Quadtree::RangeSearch(int **pindices,int *pnobs,double x,double y,double range){/*{{{*/

	/*Intermediaries*/
	int  nobs;
	int *indices = NULL;

	/*Allocate indices (maximum by default*/
	if(this->NbObs) indices = xNew<int>(this->NbObs);
	nobs = 0;

	if(this->root) this->root->RangeSearch(indices,&nobs,x,y,range);

	/*Clean-up and return*/
	*pnobs=nobs;
	*pindices=indices;

}/*}}}*/
void Quadtree::QuadtreeDepth(int* A,int xi,int yi){/*{{{*/

	QuadtreeBox **pbox = NULL;
	QuadtreeBox  *box  = NULL;
	int           level,levelbin;

	/*Initialize levels*/
	level    = 0;
	levelbin = (1L<<this->MaxDepth);// = 2^30

	/*Get inital box (the largest)*/
	pbox=&root;

	/*Find the smallest box where this point is located*/
	while((box=*pbox) && (box->nbitems<0)){ 

		levelbin>>=1; level+=1; _assert_(level<this->MaxDepth);

		pbox = &box->box[IJ(xi,yi,levelbin)];
	}
	if(box && box->nbitems>0){
		/*This box is not empty, add one level*/
		level+=1;
	}

	*A=level;
}/*}}}*/
void Quadtree::QuadtreeDepth2(int* A,int xi,int yi){/*{{{*/

	QuadtreeBox **pbox = NULL;
	QuadtreeBox  *box  = NULL;
	int           level,levelbin;

	/*Initialize levels*/
	level    = 0;
	levelbin = (1L<<this->MaxDepth);// = 2^30

	/*Get inital box (the largest)*/
	pbox=&root;

	/*Find the smallest box where this point is located*/
	while((box=*pbox) && (box->nbitems<0)){ 

		levelbin>>=1; level+=1; 

		pbox = &box->box[IJ(xi,yi,levelbin)];
	}
	if(box && box->nbitems>0){
		/*This box is not empty, add one level*/
		level+=1;
	}

	/*If we were to add the vertex, get level*/
	if(box && box->nbitems==4){
		int ij;
		bool flag=true;
		while(flag){

			levelbin>>=1; level+=1;
			if(level>this->MaxDepth){
				level+=1;
				break;
			}

			/*loop over the four observations*/
			ij=IJ(box->obs[0]->xi,box->obs[0]->yi,levelbin);
			for (int k=1;k<4;k++){
				if(IJ(box->obs[k]->xi,box->obs[k]->yi,levelbin) != ij){
					flag = false;
				}
			}
			if(IJ(xi,yi,levelbin)!=ij){
				flag = false;
			}
		}
	}

	*A=level;
}/*}}}*/

/*QuadtreeBox methos*/
Object* Quadtree::QuadtreeBox::copy(void){/*{{{*/

	   QuadtreeBox* qtreebox = new QuadtreeBox(*this);

		for (int i=0; i<4; ++i){
			if(this->box[i]) qtreebox->box[i] = reinterpret_cast<QuadtreeBox*>(this->box[i]->copy());
			else qtreebox->box[i] = NULL;
		}
		for (int i=0; i<4; ++i){
			if(this->obs[i]) qtreebox->obs[i] = reinterpret_cast<Observation*>(this->obs[i]->copy());
			else qtreebox->obs[i] = NULL;
		}

		return (Object*) qtreebox;
}
/*}}}*/
void  Quadtree::QuadtreeBox::Echo(void){/*{{{*/

	_printf_("QuadtreeBox:\n");
	_printf_("   nbitems = " << this->nbitems << "\n");
	_printf_("   xcenter = " << this->xcenter << "\n");
	_printf_("   ycenter = " << this->ycenter << "\n");
	_printf_("   length  = " << this->length << "\n");

}/*}}}*/
int Quadtree::QuadtreeBox::IsWithinRange(double x,double y,double range){/*{{{*/

	/*Return 0 if the 2 boxes do not overlap*/
	if(this->xcenter+this->length/2 < x-range) return 0;
	if(this->xcenter-this->length/2 > x+range) return 0;
	if(this->ycenter+this->length/2 < y-range) return 0;
	if(this->ycenter-this->length/2 > y+range) return 0;

	/*Return 2 if the this box is included in the range*/
	if(this->xcenter+this->length/2 <= x+range &&
		this->ycenter+this->length/2 <= y+range &&
		this->xcenter-this->length/2 >= x-range &&
		this->ycenter-this->length/2 >= y-range) return 2;

	/*This is a simple overlap*/
	return 1;

}/*}}}*/
void Quadtree::QuadtreeBox::RangeSearch(int* indices,int *pnobs,double x,double y,double range){/*{{{*/

	/*Intermediaries*/
	int i,nobs;

	/*Recover current number of observations*/
	nobs = *pnobs;

	switch(this->IsWithinRange(x,y,range)){
		case 0:
			/*If this box is not within range, return*/
			break;
		case 2:
			/*This box is included in range*/
			this->WriteObservations(indices,&nobs);
			break;
		case 1:
			/*This box is partly included*/
			if(this->nbitems>0){
				/*If this box has only observations, add indices that are within range*/
				for(i=0;i<this->nbitems;i++){
					if(fabs(this->obs[i]->x-x) <= range && fabs(this->obs[i]->y-y) <= range){
						indices[nobs++]=this->obs[i]->index;
					}
				}
			}
			else{
				/*This box points toward boxes*/
				if(this->box[0]) this->box[0]->RangeSearch(indices,&nobs,x,y,range);
				if(this->box[1]) this->box[1]->RangeSearch(indices,&nobs,x,y,range);
				if(this->box[2]) this->box[2]->RangeSearch(indices,&nobs,x,y,range);
				if(this->box[3]) this->box[3]->RangeSearch(indices,&nobs,x,y,range);
			}
			break;
		default:
			_error_("Case " << this->IsWithinRange(x,y,range) << " not supported");
	}

	/*Assign output pointers: */
	*pnobs=nobs;
}/*}}}*/
void Quadtree::QuadtreeBox::WriteObservations(int* indices,int *pnobs){/*{{{*/

	/*Intermediaries*/
	int i,nobs;

	/*Recover current number of observations*/
	nobs = *pnobs;

	if(this->nbitems>0){
		/*If this box has only observations, add all indices*/
		for(i=0;i<this->nbitems;i++){
			indices[nobs++]=this->obs[i]->index;
		}
	}
	else{
		/*This box points toward boxes, */
		if(this->box[0]) this->box[0]->WriteObservations(indices,&nobs);
		if(this->box[1]) this->box[1]->WriteObservations(indices,&nobs);
		if(this->box[2]) this->box[2]->WriteObservations(indices,&nobs);
		if(this->box[3]) this->box[3]->WriteObservations(indices,&nobs);
	}

	/*Assign output pointers: */
	*pnobs=nobs;
}/*}}}*/
