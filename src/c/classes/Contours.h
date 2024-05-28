/*!\brief Declaration of Contours class.
 */ 

#ifndef _CONTAINER_CONTOURS_H_
#define  _CONTAINER_CONTOURS_H_

#include "../datastructures/datastructures.h"
#include "./Contour.h"

class Contours: public DataSet{

	public:

		/*constructors, destructors*/
		Contours();
		~Contours();
};

/*Methods that relate to datasets: */
int ExpWrite(Contours* contours,char* domainname);
template <class doubletype> Contours* ExpRead(char* domainname){ /*{{{*/

	/*intermediary: */
	int          nprof;
	int         *profnvertices = NULL;
	doubletype **pprofx        = NULL;
	doubletype **pprofy        = NULL;

	/*If domainname is an empty string, return empty dataset*/
	if (strcmp(domainname,"")==0){
		nprof=0;
	}
	else{
		ExpRead<doubletype>(&nprof,&profnvertices,&pprofx, &pprofy, NULL,domainname);
	}

	/*now create dataset of contours: */
	Contours *domain=new Contours();

	for(int i=0;i<nprof;i++){
		domain->AddObject(new Contour<doubletype>(i,profnvertices[i],pprofx[i],pprofy[i],1));
	}
	return domain;
} /*}}}*/

#endif //ifndef _CONTOURS_H_
