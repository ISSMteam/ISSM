/*!\file:  QmuStatisticsx.h
 */ 

#ifndef _QMU_STATISTCSX_H_
#define _QMU_STATISTCSX_H_

#include "../../classes/classes.h"

int ComputeMeanVariance(Parameters* parameters,Results* results,int color, ISSM_MPI_Comm statcomm);
int ComputeSampleSeries(Parameters* parameters,Results* results,int color, ISSM_MPI_Comm statcomm);
int OutputStatistics(Parameters* parameters,Results* results,int color,ISSM_MPI_Comm statcomm);
int ComputeHistogram(Parameters* parameters,Results* results,int color, ISSM_MPI_Comm statcomm);
int readdata(IssmDouble** pdoublemat, int* pdoublematsize, IssmDouble* pdouble, FILE* fid,char* field,int step);
bool DakotaDirStructure(int argc,char** argv);
int DakotaStatistics(int argc,char** argv);

/* local prototypes: */

#endif  /* _QMU_STATISTCSX_H_ */
