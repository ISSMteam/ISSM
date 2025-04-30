/*
 * CoDiPackCommon.hpp
 *
 *  Created on: Mai 30, 2016
 *      Author: sagebaum
 */

#ifndef _CODIPACK_COMMON_HPP_
#define _CODIPACK_COMMON_HPP_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile without HAVE_CONFIG_H symbol! run configure first!"
#endif

#if defined(_HAVE_CODIPACK_)
#include "codi.hpp"

template<typename Real, typename Passive>
inline void getVectorPrimal(const Real* vec, Passive* pasVec, int n) {
  for(int i = 0; i < n; ++i) {
    pasVec[i]=vec[i].getValue();
  }
}

template<typename Real, typename Passive>
inline void setVectorPrimal(Real* vec, const Passive* pasVec, int n) {
  for(int i = 0; i < n; ++i) {
    vec[i].value() = pasVec[i];
  }
}

template<typename Real, typename Data>
inline void getVectorGradData(const Real* vec, Data* dataVec, int n) {
  for(int i = 0; i < n; ++i) {
	 #if _CODIPACK_MAJOR_==2
    dataVec[i]=vec[i].getIdentifier();
	 #elif _CODIPACK_MAJOR_==1
	 dataVec[i]=vec[i].getGradientData();
	 #else
	 #error "_CODIPACK_MAJOR_ not supported"
	 #endif
  }
}

template<typename Real, typename Passive, typename Data>
inline void getVectorPrimalAndGradData(const Real* vec, Passive* pasVec, Data* dataVec, int n) {
  for(int i = 0; i < n; ++i) {
    pasVec[i]=vec[i].getValue();
	 #if _CODIPACK_MAJOR_==2
    dataVec[i]=vec[i].getIdentifier();
	 #elif _CODIPACK_MAJOR_==1
	 dataVec[i]=vec[i].getGradientData();
	 #else
	 #error "_CODIPACK_MAJOR_ not supported"
	 #endif
  }
}

template<typename Real, typename Passive, typename Data>
inline void getPrimalAndGradData(const Real& value, Passive& pas, Data& data) {
  pas=value.getValue();
  #if _CODIPACK_MAJOR_==2
  data=value.getIdentifier();
  #elif _CODIPACK_MAJOR_==1
  data=value.getGradientData();
  #else
  #error "_CODIPACK_MAJOR_ not supported"
  #endif
}

template<typename Real, typename Data>
inline void registerVector(Real* vec, Data* dataVec, int n) {
	#if _CODIPACK_MAJOR_==2
	typename Real::Tape& tape = Real::getTape();
	#elif _CODIPACK_MAJOR_==1
	typename Real::TapeType& tape = Real::getGlobalTape();
	#else
	#error "_CODIPACK_MAJOR_ not supported"
	#endif

  for(int i = 0; i < n; ++i) {
    tape.registerInput(vec[i]);
	 #if _CODIPACK_MAJOR_==2
    dataVec[i]=vec[i].getIdentifier();
	 #elif _CODIPACK_MAJOR_==1
	 dataVec[i]=vec[i].getGradientData();
	 #else
	 #error "_CODIPACK_MAJOR_ not supported"
	 #endif
  }
}

template<typename Tape, typename Data, typename Adjoint>
inline void getVectorAdjoint(Tape& tape, const Data* dataVec, Adjoint* adjVec, int n) {
  for(int i = 0; i < n; ++i) {
    Data index = dataVec[i];
    adjVec[i] = tape.getGradient(index);
    tape.gradient(index) = 0.0;
  }
}

template<typename Tape, typename Data, typename Adjoint>
inline void updateVectorAdjoint(Tape& tape, const Data* dataVec, const Adjoint* adjVec, int n) {
  for(int i = 0; i < n; ++i) {
    Data index = dataVec[i];
    if(0 != index) {
      tape.gradient(index) += adjVec[i];
    }
  }
}

template<typename Tape, typename Data, typename Adjoint>
inline void updateAdjoint(Tape& tape, const Data& data, const Adjoint& adj) {
  Data index = data;
  if(0 != index) {
    tape.gradient(index) += adj;
  }
}

#endif

#endif
