#include "./bamgobjects.h"

using namespace std;
namespace bamg {

	/*Constructor*/
	SetOfEdges4::SetOfEdges4(long mmx,long nnx){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, SetOfEdges4.cpp/SetOfEdges4)*/

		/*Intermediary*/
		int i;

		//initialize fields
		nx   =nnx;   //number of vertices
		nbax =mmx;   // 3 * number of triangles
		NbOfEdges=0;
		head = new long [nx];
		Edges= new IntEdge[nbax];

		//initialize head (-1 everywhere)
		i=nx;
		while(i--) head[i]=-1;
	}
	/*}}}*/

	/*Methods*/ 
	long SetOfEdges4::add(long ii,long jj) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, SetOfEdges4.cpp/add)*/

		/*Intermediary*/
		int h,n;

		//get n from h (usually h=ii)
		_assert_(head);
		n=head[h=Abs(ii)%nx];

		//go through the existing edges that holds h (=ii) and check that 
		//the edge ii jj is not already in Edge
		while (n >= 0){

			//if the edge ii jj is already in Edges, return n
			if (ii == Edges[n].i && jj == Edges[n].j) return n;

			//else go to next edge that holds ii
			else n = Edges[n].next;
		}

		//check that nbax <=NbOfEdges
		if (nbax <=NbOfEdges ) {
			_error_("SetOfEdges4::add overflow: NbOfEdges=" << NbOfEdges << " > nbax=" << nbax);
		}

		//update chain
		Edges[NbOfEdges].i=ii;
		Edges[NbOfEdges].j=jj;
		Edges[NbOfEdges].next= head[h];
		head[h] = NbOfEdges;
		return NbOfEdges ++;
	}
	/*}}}*/
	long SetOfEdges4::find(long ii,long jj) { /*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, SetOfEdges4.cpp/find)*/

		/*Intermediary*/
		int n;

		//check that head is not empty
		_assert_(head);

		//get n from h (usually h=ii)
		n=head[Abs(ii)%nx];

		//go through the existing edges that holds h (=ii) and return position in Edge
		while (n >= 0){

			//if the edge ii jj is already in Edges, return n
			if (ii == Edges[n].i && jj == Edges[n].j) return n;

			//else go to next edge that holds ii
			else n = Edges[n].next;
		}

		//if we reach this point, the edge does not exist return -1
		return -1;
	}
	/*}}}*/
	long SetOfEdges4::i(long k){/*{{{*/
		return Edges[k].i;
	}
	/*}}}*/
	long SetOfEdges4::j(long k){/*{{{*/
		return Edges[k].j;
	}
	/*}}}*/
	long SetOfEdges4::nb(){/*{{{*/
		return NbOfEdges;
	}
	/*}}}*/
	long SetOfEdges4::SortAndAdd (long ii,long jj) {/*{{{*/
		return ii <=jj ? add (ii,jj)  : add (jj,ii) ;
	}
	/*}}}*/
	long SetOfEdges4::SortAndFind (long ii,long jj) {/*{{{*/
		return ii <=jj ? find (ii,jj)  : find (jj,ii) ;
	}
	/*}}}*/
}
