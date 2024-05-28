/*
 * CoDiPackGlobal.h
 *
 *  Created on: Jul 20, 2016
 *      Author: a_h_ck
 */

#ifndef SRC_C_TOOLKITS_CODIPACK_CODIPACKGLOBAL_H_
#define SRC_C_TOOLKITS_CODIPACK_CODIPACKGLOBAL_H_

#ifdef HAVE_CONFIG_H
    #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#if defined(_HAVE_CODIPACK_)

#include <iostream>
#include <vector>

struct CoDi_global {
    std::vector<int> input_indices;
    std::vector<int> output_indices;

    void print(std::ostream& out) const {
        out << "-------------------------------------\nCoDi_global:\n  in = [ ";
        for(auto& in_index : input_indices) {
            out << in_index << " ";
        }
        out << "]\n  out = [ ";
        for(auto& out_index : output_indices) {
            out << out_index << " ";
        }
        out << "]\n-------------------------------------\n";
    }
};

#endif /* _HAVE_CODIPACK_ */
#endif /* SRC_C_TOOLKITS_CODIPACK_CODIPACKGLOBAL_H_ */
