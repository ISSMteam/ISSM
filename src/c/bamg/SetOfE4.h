#ifndef _SetOfEdge4_h
#define _SetOfEdge4_h

namespace bamg {

	class SetOfEdges4;

	class IntEdge{
		friend class SetOfEdges4;
		public:
		long i,j;
		long next; 
	};

	class SetOfEdges4 {

		private:
			long nx,nbax,NbOfEdges;
			long* head; 
			IntEdge* Edges;

		public:

			// Constructors
			SetOfEdges4(long ,long);// nb Edges mx , nb de sommet 
			~SetOfEdges4() {delete [] head; delete [] Edges;}

			//operators
			IntEdge & operator[](long k){return  Edges[k];}

			//Methods
			long add (long ii,long jj);
			long SortAndAdd (long ii,long jj);
			long nb();
			long find (long ii,long jj);
			long SortAndFind (long ii,long jj);
			long i(long k);
			long j(long k);
	};
}
#endif 
