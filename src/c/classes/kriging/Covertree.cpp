#include "../classes.h"
#include <set>
#include <algorithm>

	/*Constructors/Destructors*/
Covertree::Covertree(int maxLevel,const std::vector<Observation>& points){/*{{{*/
	this->base = 2.;
	_root=NULL;
	_numNodes=0;
	_maxLevel=maxLevel;//ceilf(log(maxDist)/log(base));
	_minLevel=_maxLevel-1;
	std::vector<Observation>::const_iterator it;
	for(it=points.begin(); it!=points.end(); ++it) {
		this->insert(*it);//adds data to the covertree object
	}
}/*}}}*/
Covertree::~Covertree(){/*{{{*/
	if(_root==NULL) return;
	//Get all of the root's children (from any level),
	//delete the root, repeat for each of the children
	std::vector<CoverTreeNode*> nodes;
	nodes.push_back(_root);
	while(!nodes.empty()) {
		CoverTreeNode* byeNode = nodes[0];
		nodes.erase(nodes.begin());
		std::vector<CoverTreeNode*> children = byeNode->getAllChildren();
		nodes.insert(nodes.begin(),children.begin(),children.end());
		delete byeNode;
	}   
}/*}}}*/

	/*Methods*/
std::pair<double, Covertree::CoverTreeNode*>		Covertree::distance(const Observation& p, const std::vector<CoverTreeNode*>& Q){/*{{{*/
	double minDist = 1.e+50;
	CoverTreeNode* minNode;
	std::vector<CoverTreeNode*>::const_iterator it;
	for(it=Q.begin();it!=Q.end();++it) {
		double dist = p.distance((*it)->getObservation());
		if(dist < minDist) {
			minDist = dist;
			minNode = *it;
		}
	}
	return std::make_pair(minDist,minNode);  
}/*}}}*/
int Covertree::get_numberofobs(){/*{{{*/
	return _numNodes;
}/*}}}*/
void   Covertree::insert(const Observation& newObservation){/*{{{*/
	if(_root==NULL) {
		_root = new CoverTreeNode(newObservation);
		_numNodes=1;
		return;
	}
	//TODO: this is pretty inefficient, there may be a better way
	//to check if the node already exists...
	CoverTreeNode* n = kNearestNodes(newObservation,1)[0];
	if(newObservation.distance(n->getObservation())==0.0) {
		n->addObservation(newObservation);
	} else {
		//insert_rec acts under the assumption that there are no nodes with
		//distance 0 to newObservation in the cover tree (the previous lines check it)
		insert_rec(newObservation,
					std::vector<distNodePair>
					(1,std::make_pair(_root->distance(newObservation),_root)),
					_maxLevel);
	}
}/*}}}*/
bool		Covertree::insert_rec(const Observation& p, const std::vector<distNodePair>& Qi, const int& level){/*{{{*/
	std::vector<std::pair<double, CoverTreeNode*> > Qj;
	double sep = pow(base,level);
	double minDist = 1.e+50;
	std::pair<double,CoverTreeNode*> minQiDist(1.e+50,NULL);
	std::vector<std::pair<double, CoverTreeNode*> >::const_iterator it;
	for(it=Qi.begin(); it!=Qi.end(); ++it) {
		if(it->first<minQiDist.first) minQiDist = *it;
		if(it->first<minDist) minDist=it->first;
		if(it->first<=sep) Qj.push_back(*it);
		std::vector<CoverTreeNode*> children = it->second->getChildren(level);
		std::vector<CoverTreeNode*>::const_iterator it2;
		for(it2=children.begin();it2!=children.end();++it2) {
			double d = p.distance((*it2)->getObservation());
			if(d<minDist) minDist = d;
			if(d<=sep) {
				Qj.push_back(std::make_pair(d,*it2));
			}
		}
	}
	//std::cout << "level: " << level << ", sep: " << sep << ", dist: " << minQDist.first << "\n";
	if(minDist > sep) {
		return true;
	} else {
		bool found = insert_rec(p,Qj,level-1);
		//distNodePair minQiDist = distance(p,Qi);
		if(found && minQiDist.first <= sep) {
			if(level-1<_minLevel) _minLevel=level-1;
			minQiDist.second->addChild(level,
						new CoverTreeNode(p));
			//std::cout << "parent is ";
			//minQiDist.second->getObservation().print();
			_numNodes++;
			return false;
		} else {
			return found;
		}
	}
}/*}}}*/
std::vector<Covertree::CoverTreeNode*> Covertree::kNearestNodes(const Observation& p, const unsigned int& k) const{/*{{{*/
	if(_root==NULL) return std::vector<CoverTreeNode*>();
	//maxDist is the kth nearest known point to p, and also the farthest
	//point from p in the set minNodes defined below.
	double maxDist = p.distance(_root->getObservation());
	//minNodes stores the k nearest known points to p.
	std::set<distNodePair> minNodes;

	minNodes.insert(std::make_pair(maxDist,_root));
	std::vector<distNodePair> Qj(1,std::make_pair(maxDist,_root));
	for(int level = _maxLevel; level>=_minLevel;level--) {
		std::vector<distNodePair>::const_iterator it;
		int size = Qj.size();
		for(int i=0; i<size; i++) {
			std::vector<CoverTreeNode*> children =
			  Qj[i].second->getChildren(level);
			std::vector<CoverTreeNode*>::const_iterator it2;
			for(it2=children.begin(); it2!=children.end(); ++it2) {
				double d = p.distance((*it2)->getObservation());
				if(d < maxDist || minNodes.size() < k) {
					minNodes.insert(std::make_pair(d,*it2));
					//--minNodes.end() gives us an iterator to the greatest
					//element of minNodes.
					if(minNodes.size() > k) minNodes.erase(--minNodes.end());
					maxDist = (--minNodes.end())->first;
				}
				Qj.push_back(std::make_pair(d,*it2));
			}
		}
		double sep = maxDist + pow(base, level);
		size = Qj.size();
		for(int i=0; i<size; i++) {
			if(Qj[i].first > sep) {
				//quickly removes an element from a vector w/o preserving order.
				Qj[i]=Qj.back();
				Qj.pop_back();
				size--; i--;
			}
		}
	}
	std::vector<CoverTreeNode*> kNN;
	std::set<distNodePair>::const_iterator it;
	for(it=minNodes.begin();it!=minNodes.end();++it) {
		kNN.push_back(it->second);
	}
	return kNN;
}/*}}}*/
std::vector<Observation> Covertree::kNearestNeighbors(const Observation& p, const unsigned int& k) const{/*{{{*/
	if(_root==NULL) return std::vector<Observation>();
	std::vector<CoverTreeNode*> v = kNearestNodes(p, k);
	std::vector<Observation> kNN;
	std::vector<CoverTreeNode*>::const_iterator it;
	for(it=v.begin();it!=v.end();++it) {
		const std::vector<Observation>& p = (*it)->getObservations();
		kNN.insert(kNN.end(),p.begin(),p.end());
		if(kNN.size() >= k) break;
	}
	return kNN;
}/*}}}*/
void   Covertree::print() const{/*{{{*/
	int d = _maxLevel-_minLevel+1;
	std::vector<CoverTreeNode*> Q;
	Q.push_back(_root);
	for(int i=0;i<d;i++) {
		std::cout << "LEVEL " << _maxLevel-i << "\n";
		std::vector<CoverTreeNode*>::const_iterator it;
		for(it=Q.begin();it!=Q.end();++it) {
			(*it)->getObservation().print();
			std::vector<CoverTreeNode*>
			  children = (*it)->getChildren(_maxLevel-i);
			std::vector<CoverTreeNode*>::const_iterator it2;
			for(it2=children.begin();it2!=children.end();++it2) {
				std::cout << "  ";
				(*it2)->getObservation().print();
			}
		}
		std::vector<CoverTreeNode*> newQ;
		for(it=Q.begin();it!=Q.end();++it) {
			std::vector<CoverTreeNode*>
			  children = (*it)->getChildren(_maxLevel-i);
			newQ.insert(newQ.end(),children.begin(),children.end());
		}
		Q.insert(Q.end(),newQ.begin(),newQ.end());
		std::cout << "\n\n";
	}
}/*}}}*/
void   Covertree::remove(const Observation& p){/*{{{*/
	//Most of this function's code is for the special case of removing the root
	if(_root==NULL) return;
	bool removingRoot=_root->hasObservation(p);
	if(removingRoot && !_root->isSingle()) {
		_root->removeObservation(p);
		return;
	}
	CoverTreeNode* newRoot=NULL;
	if(removingRoot) {
		if(_numNodes==1) {
			//removing the last node...
			delete _root;
			_numNodes--;
			_root=NULL;
			return;
		} else {
			for(int i=_maxLevel;i>_minLevel;i--) {
				if(!(_root->getChildren(i).empty())) {
					newRoot = _root->getChildren(i).back();
					_root->removeChild(i,newRoot);
					break;
				}
			}
		}
	}
	std::map<int, std::vector<distNodePair> > coverSets;
	coverSets[_maxLevel].push_back(std::make_pair(_root->distance(p),_root));
	if(removingRoot)
	 coverSets[_maxLevel].push_back(std::make_pair(newRoot->distance(p),newRoot));
	bool multi = false;
	remove_rec(p,coverSets,_maxLevel,multi);
	if(removingRoot) {
		delete _root;
		_numNodes--;
		_root=newRoot;
	}
}/*}}}*/
void   Covertree::remove_rec(const Observation& p, std::map<int,std::vector<distNodePair> >& coverSets, int level, bool& multi){/*{{{*/
	std::vector<distNodePair>& Qi = coverSets[level];
	std::vector<distNodePair>& Qj = coverSets[level-1];
	double minDist = 1.e+50;
	CoverTreeNode* minNode = _root;
	CoverTreeNode* parent = 0;
	double sep = pow(base, level);
	std::vector<distNodePair>::const_iterator it;
	//set Qj to be all children q of Qi such that p.distance(q)<=sep
	//and also keep track of the minimum distance from p to a node in Qj
	//note that every node has itself as a child, but the
	//getChildren function only returns non-self-children.
	for(it=Qi.begin();it!=Qi.end();++it) {
		std::vector<CoverTreeNode*> children = it->second->getChildren(level);
		double dist = it->first;
		if(dist<minDist) {
			minDist = dist;
			minNode = it->second;
		}
		if(dist <= sep) {
			Qj.push_back(*it);
		}
		std::vector<CoverTreeNode*>::const_iterator it2;
		for(it2=children.begin();it2!=children.end();++it2) {
			dist = p.distance((*it2)->getObservation());
			if(dist<minDist) {
				minDist = dist;
				minNode = *it2;
				if(dist == 0.0) parent = it->second;
			}
			if(dist <= sep) {
				Qj.push_back(std::make_pair(dist,*it2));
			}
		}
	}
	if(level>_minLevel) remove_rec(p,coverSets,level-1,multi);
	if(minNode->hasObservation(p)) {
		//the multi flag indicates the point we removed is from a
		//node containing multiple points, and we have removed it,
		//so we don't need to do anything else.
		if(multi) return;
		if(!minNode->isSingle()) {
			minNode->removeObservation(p);
			multi=true;
			return;
		}
		if(parent!=NULL) parent->removeChild(level, minNode);
		std::vector<CoverTreeNode*> children = minNode->getChildren(level-1);
		std::vector<distNodePair>& Q = coverSets[level-1];
		if(Q.size()==1 && Q[0].second==minNode) {
			Q.pop_back();
		} else {
			for(unsigned int i=0;i<Q.size();i++) {
				if(Q[i].second==minNode) {
					Q[i]=Q.back();
					Q.pop_back();
					break;
				}
			}
		}
		std::vector<CoverTreeNode*>::const_iterator it;
		for(it=children.begin();it!=children.end();++it) {
			int i = level-1;
			Observation q = (*it)->getObservation();
			double minDQ = 1.e+50;
			CoverTreeNode* minDQNode = NULL;
			double sep = pow(base,i);
			bool br=false;
			while(true) {
				std::vector<distNodePair>&
				  Q = coverSets[i];
				std::vector<distNodePair>::const_iterator it2;
				minDQ = 1.e+50;
				for(it2=Q.begin();it2!=Q.end();++it2) {
					double d = q.distance(it2->second->getObservation());
					if(d<minDQ) {
						minDQ = d;
						minDQNode = it2->second;
						if(d <=sep) {
							br=true;
							break;
						}
					}
				}
				minDQ=1.e+50;
				if(br) break;
				Q.push_back(std::make_pair((*it)->distance(p),*it));
				i++;
				sep = pow(base,i);
			}
			//minDQNode->getObservation().print();
			//std::cout << " is level " << i << " parent of ";
			//(*it)->getObservation().print();
			if (minDQNode != NULL)
			 minDQNode->addChild(i,*it);
		}
		if(parent!=NULL) {
			delete minNode;
			_numNodes--;
		}
	}
}/*}}}*/

void   Covertree::CoverTreeNode::addChild(int level, CoverTreeNode* p){/*{{{*/
	_childMap[level].push_back(p);
}/*}}}*/
void   Covertree::CoverTreeNode::addObservation(const Observation& p){/*{{{*/
	if(find(_observations.begin(), _observations.end(), p) == _observations.end())
	 _observations.push_back(p);
}/*}}}*/
Covertree::CoverTreeNode::CoverTreeNode(const Observation& p) {/*{{{*/
	_observations.push_back(p);
}/*}}}*/
double Covertree::CoverTreeNode::distance(const CoverTreeNode& p) const{/*{{{*/
	return _observations[0].distance(p.getObservation());
}/*}}}*/
std::vector<Covertree::CoverTreeNode*> Covertree::CoverTreeNode::getAllChildren() const{/*{{{*/
	std::vector<CoverTreeNode*> children;
	std::map<int,std::vector<CoverTreeNode*> >::const_iterator it;
	for(it=_childMap.begin();it!=_childMap.end();++it) {
		children.insert(children.end(), it->second.begin(), it->second.end());
	}
	return children;
}/*}}}*/
std::vector<Covertree::CoverTreeNode*> Covertree::CoverTreeNode::getChildren(int level) const{/*{{{*/
	std::map<int,std::vector<CoverTreeNode*> >::const_iterator
	  it = _childMap.find(level);
	if(it!=_childMap.end()) {
		return it->second;
	}
	return std::vector<CoverTreeNode*>();
}/*}}}*/
const Observation& Covertree::CoverTreeNode::getObservation() const{/*{{{*/
	return _observations[0]; 
}/*}}}*/
Covertree::CoverTreeNode* Covertree::getRoot() const{/*{{{*/
	return _root;
}/*}}}*/
bool   Covertree::CoverTreeNode::hasObservation(const Observation& p) const{/*{{{*/
	return find(_observations.begin(), _observations.end(), p) != _observations.end();
}/*}}}*/
bool   Covertree::CoverTreeNode::isSingle() const{/*{{{*/
	return _observations.size() == 1;
}/*}}}*/
bool   Covertree::isValidTree() const {/*{{{*/
	if(_numNodes==0)
	 return _root==NULL;

	std::vector<CoverTreeNode*> nodes;
	nodes.push_back(_root);
	for(int i=_maxLevel;i>_minLevel;i--) {
		double sep = pow(base,i);
		std::vector<CoverTreeNode*>::const_iterator it, it2;
		//verify separation invariant of cover tree: for each level,
		//every point is farther than base^level away
		for(it=nodes.begin(); it!=nodes.end(); ++it) {
			for(it2=nodes.begin(); it2!=nodes.end(); ++it2) {
				double dist=(*it)->distance((*it2)->getObservation());
				if(dist<=sep && dist!=0.0) {
					std::cout << "Level " << i << " Separation invariant failed.\n";
					return false;
				}
			}
		}
		std::vector<CoverTreeNode*> allChildren;
		for(it=nodes.begin(); it!=nodes.end(); ++it) {        
			std::vector<CoverTreeNode*> children = (*it)->getChildren(i);
			//verify covering tree invariant: the children of node n at level
			//i are no further than base^i away
			for(it2=children.begin(); it2!=children.end(); ++it2) {
				double dist = (*it2)->distance((*it)->getObservation());
				if(dist>sep) {
					std::cout << "Level" << i << " covering tree invariant failed.n";
					return false;
				}
			}
			allChildren.insert
			  (allChildren.end(),children.begin(),children.end());
		}
		nodes.insert(nodes.begin(),allChildren.begin(),allChildren.end());
	}
	return true;
}/*}}}*/
void   Covertree::CoverTreeNode::removeChild(int level, CoverTreeNode* p){/*{{{*/
	std::vector<CoverTreeNode*>& v = _childMap[level];
	for(unsigned int i=0;i<v.size();i++) {
		if(v[i]==p) {
			v[i]=v.back();
			v.pop_back();
			break;
		}
	}
}/*}}}*/
void   Covertree::CoverTreeNode::removeObservation(const Observation& p){/*{{{*/
	std::vector<Observation>::iterator it =
	  find(_observations.begin(), _observations.end(), p);
	if(it != _observations.end())
	 _observations.erase(it);
}/*}}}*/
