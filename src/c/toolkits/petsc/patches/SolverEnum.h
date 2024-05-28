/*!\file: SolverEnum.h
 * \brief prototypes for SolverEnum.h
 */ 

#ifndef _SOLVERENUM_H_
#define  _SOLVERENUM_H_

typedef enum{
	PETSCPACKAGE,
	MUMPSPACKAGE_LU,
	MUMPSPACKAGE_CHOL,
	SPOOLESPACKAGE_LU,
	SPOOLESPACKAGE_CHOL,
	SUPERLUDISTPACKAGE,
} EXTERNALPACKAGES; 

#endif //ifndef _SOLVERENUM_H_
