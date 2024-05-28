/*!\file BamgTriangulatex
 */

#include "./BamgTriangulatex.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../bamg/bamgobjects.h"

using namespace bamg;
using namespace std;

int BamgTriangulatex(int** pindex,int* pnels,double* x,double* y,int nods){

	BamgOpts* bamgopts=new BamgOpts();//use bamgopts->verbose>5 to debug bamg::Mesh()
	Mesh Th(x,y,nods,bamgopts);
	Th.WriteIndex(pindex,pnels);
	delete bamgopts;
	//delete &Th;
	return 0;
}
