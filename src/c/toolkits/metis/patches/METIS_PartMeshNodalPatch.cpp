/*!\file METIS_PartMeshNodalPatch
 * \brief common interface to Metis 4.0 and 5.0
 */

#include <config.h>
#include "../metisincludes.h"
#include "../../../shared/shared.h"

/*METIS prototypes*/
extern "C" {
#if _METIS_VERSION_ == 4
	void METIS_PartMeshNodal(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
#endif
#if _METIS_VERSION_ == 5
	int METIS_PartMeshNodal(idx_t*, idx_t*, idx_t*, idx_t*, idx_t*, idx_t*, idx_t*, real_t*, idx_t*, idx_t*, idx_t*, idx_t*);
	int METIS_SetDefaultOptions(idx_t *options);
#endif
}

void METIS_PartMeshNodalPatch(int numberofelements,int numberofnodes,int* index,int* vweights,int num_procs,int* epart,int* npart){

	#if _METIS_VERSION_ == 4
	/*Our interface originates in the Metis 4.0 version, hence identical calls*/
	int  edgecut=1;
	int  etype  =1; //tria mesh see metis/Programs/Io.c
	int  numflag=0;
	METIS_PartMeshNodal(&numberofelements,&numberofnodes, index,&etype,&numflag,&num_procs,&edgecut, epart, npart); 

	#elif _METIS_VERSION_ == 5

	/*Create options*/
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);

	options[METIS_OPTION_PTYPE]   = METIS_PTYPE_KWAY;     /* partitioning method  */
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;    /* type of objective */
	options[METIS_OPTION_CTYPE]   = METIS_CTYPE_SHEM;     /* matching scheme to be used during coarsening.*/
	options[METIS_OPTION_IPTYPE]  = METIS_IPTYPE_METISRB; /* algorithm used during initial partitioning*/
	options[METIS_OPTION_RTYPE]   = METIS_RTYPE_GREEDY;   /* algorithm used for refinement*/
	options[METIS_OPTION_DBGLVL]  = 0;                    /* amount of progress/debugging information will be printed */
	options[METIS_OPTION_UFACTOR] = 30;                   /* maximum allowed load imbalance among the partitions*/
	options[METIS_OPTION_MINCONN] = 0;                    /* explicitly minimize the maximum connectivity ?*/
	options[METIS_OPTION_CONTIG]  = 0;                    /* force contiguous partitions?*/
	options[METIS_OPTION_SEED]    = -1;                   /* seed for the random number generator*/
	options[METIS_OPTION_NITER]   = 10;                   /* number of iterations for the refinement algorithms*/
	options[METIS_OPTION_NCUTS]   = 1;                    /* number of different partitionings that it will compute*/

	/*create eptr*/
	idx_t  k=0;
	idx_t* eptr=xNew<idx_t>(numberofelements+1);
	eptr[0]=0;
	for(int i=0;i<numberofelements;i++){
		k+=3;
		eptr[i+1]=k;
	}

	/*create tpwgts (Weight per processor)*/
	real_t* tpwgts=xNew<real_t>(num_procs);
	for(int i=0;i<num_procs;i++) tpwgts[i]=1.0/(num_procs);

	/*create vwgt (Weight per node)*/
	idx_t* vwgts=NULL;
	if(vweights){
		vwgts=xNew<idx_t>(numberofnodes);
		for(int i=0;i<numberofnodes;i++) vwgts[i]=reCast<idx_t>(vweights[i]);
	}

	/*Call METIS*/
	idx_t objval;
	int output = METIS_PartMeshNodal(&numberofelements,&numberofnodes,eptr,index,vwgts,NULL,&num_procs,tpwgts,options,&objval,epart,npart);
	if(output!=METIS_OK) _error_("Could not partition mesh");

	/*clean-up*/
	xDelete<idx_t>(vwgts);
	xDelete<idx_t>(eptr);
	xDelete<real_t>(tpwgts);

	#else
	_error_("METIS version not supported yet");
	#endif
}
