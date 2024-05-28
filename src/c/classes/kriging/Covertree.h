
#ifndef _COVERTREE_H
#define _COVERTREE_H

#include <map>
class Observation;

class Covertree{

	/* Cover tree node. Consists of arbitrarily many points P, as long as they
	 * have distance 0 to each other. Keeps track of its children.  */
	class CoverTreeNode{
		private:
			//_childMap[i] is a vector of the node's children at level i
			std::map<int,std::vector<CoverTreeNode*> > _childMap;
			//_observations is all of the points with distance 0 which are not equal.
			std::vector<Observation> _observations;
		public:
			CoverTreeNode(const Observation& o);
			/**
			 * Returns the children of the node at level i. Note that this means
			 * the children exist in cover set i-1, not level i.
			 *
			 * Does not include the node itself, though technically every node
			 * has itself as a child in a cover tree.
			 */
			void addChild(int level, CoverTreeNode* p);
			void addObservation(const Observation& o);
			double distance(const CoverTreeNode& p) const;
			bool   isSingle() const;
			bool   hasObservation(const Observation& o) const;
			std::vector<CoverTreeNode*> getChildren(int level) const;
			const Observation& getObservation() const;
			const std::vector<Observation>& getObservations() { return _observations; }
			void removeChild(int level, CoverTreeNode* p);
			void removeObservation(const Observation& o);

			/**
			 * Return every child of the node from any level. This is handy for
			 * the destructor.
			 */
			std::vector<CoverTreeNode*> getAllChildren() const;
	  }; // CoverTreeNode class
	private:
	typedef std::pair<double, CoverTreeNode*> distNodePair;

	CoverTreeNode *_root;
	unsigned int   _numNodes;
	int            _maxLevel;   //base^_maxLevel should be the max distance
	//between any 2 points
	int            _minLevel;   //A level beneath which there are no more new nodes.

	/* Finds the node in Q with the minimum distance to p. Returns a pair
	 * consisting of this node and the distance.  */
	distNodePair distance(const Observation& p,const std::vector<CoverTreeNode*>& Q);
	/**
	 * Recursive implementation of the insert algorithm (see paper).
	 */
	bool insert_rec(const Observation& p, const std::vector<distNodePair>& Qi,const int& level);

	std::vector<CoverTreeNode*> kNearestNodes(const Observation& o, const unsigned int& k) const;
	void remove_rec(const Observation& p, std::map<int,std::vector<distNodePair> >& coverSets, int level, bool& multi);

	public:
	double base;

	/**
	 * Constructs a cover tree which begins with all points in points.
	 *
	 * maxDist should be the maximum distance that any two points
	 * can have between each other. IE p.distance(q) < maxDist for all
	 * p,q that you will ever try to insert. The cover tree may be invalid
	 * if an inaccurate maxDist is given.
	 */

	Covertree(int maxDist,const std::vector<Observation>& points=std::vector<Observation>()); 
	~Covertree();

	/**
	 * Insert newPoint into the cover tree. If newPoint is already present,
	 * (that is, newPoint==p for some p already in the tree), then the tree
	 * is unchanged. If p.distance(newPoint)==0.0 but newPoint!=p, then
	 * newPoint WILL be inserted and both points may be returned in k-nearest-
	 * neighbor searches.
	 */
	void insert(const Observation& newObservation);

	/**
	 * Just for testing/debugging. Returns true iff the cover tree satisfies the
	 * the covering tree invariants (every node in level i is greater than base^i
	 * distance from every other node, and every node in level i is less than
	 * or equal to base^i distance from its children). See the cover tree
	 * papers for details.
	 */
	bool isValidTree() const;

	/**
	 * Remove point p from the cover tree. If p is not present in the tree,
	 * it will remain unchanged. Otherwise, this will remove exactly one
	 * point q from the tree satisfying p==q.
	 */
	void remove(const Observation& p);

	/**
	 * Returns the k nearest points to p in order (the 0th element of the vector
	 * is closest to p, 1th is next, etc). It may return greater than k points
	 * if there is a tie for the kth place.
	 */
	std::vector<Observation> kNearestNeighbors(const Observation& p, const unsigned int& k) const;

	int get_numberofobs();

	CoverTreeNode* getRoot() const;

	/**
	 * Print the cover tree.
	 */
	void print() const;
};
#endif //_COVERTREE_H
