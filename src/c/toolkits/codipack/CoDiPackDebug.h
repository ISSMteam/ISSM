#ifndef SRC_C_TOOLKITS_CODIPACK_CODIPACKDEBUG_H_
#define SRC_C_TOOLKITS_CODIPACK_CODIPACKDEBUG_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#if defined(_HAVE_CODIPACK_)

#include "./CoDiPackTypes.h"

int CoDiGetUniqueID();

void CoDiStartDumpTape();
void CoDiStopDumpTape();
void CoDiStartDumpEval();
void CoDiStopDumpEval();
void CoDiPauseDumpTape();
void CoDiResumeDumpTape();

bool          CoDiDebugGetOutputPrimal();
bool          CoDiDebugGetOutputReverse();
bool          CoDiDebugGetOutputIdentifiers();
int           CoDiDebugGetOutputPrecission();
std::ostream& CoDiDebugGetOutputStream();

void CoDiDebugSetOutputPrimal(bool value);
void CoDiDebugSetOutputReverse(bool value);
void CoDiDebugSetOutputIdentifiers(bool value);
void CoDiDebugSetOutputPrecission(int value);
void CoDiDebugSetOutputStream(std::ostream& value);

bool CoDiIsDebugOutput();
void CoDiEnableDebugOutput(bool b);
bool CoDiDisableDebugOutput();

void* MatDebugOutputStart(std::string message, int M, int N);
void MatDebugOutputAddRow(void* h, int row, int size, int* cols, CoDiReal* values);
void MatDebugOutputAddRow(void*, int, int, int*, double*);
void MatDebugOutputFinish(void* h);

void VecDebugOutput(std::string message, int M, int m, CoDiReal* values);
void VecDebugOutput(std::string message, int M, int m, double* values);

void ArrayDebugOutput(std::string message, int m, CoDiReal* values);
void ArrayDebugOutput(std::string message, int m, double* values);

#endif /* _HAVE_CODIPACK_ */
#endif /* SRC_C_TOOLKITS_CODIPACK_CODIPACKDEBUG_H_ */
