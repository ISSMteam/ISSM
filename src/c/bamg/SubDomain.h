#ifndef _SUBDOMAIN_H_
#define _SUBDOMAIN_H_

#include "./include.h"
#include "./Edge.h"

namespace bamg {

	class Triangle;
	class Mesh;

	class SubDomain {

		public:

			Triangle *head;
			long      ReferenceNumber;
			int       direction;   // -1 or 1
			Edge     *edge;        // to geometrical

			//Methods
			void Set(const Mesh &,long,Mesh &);
	};

}
#endif
