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
#include <sstream>

#include "CoDiPackTypes.h"
#include "../../shared/Exceptions/exceptions.h"

struct CoDi_global {

		using Tape = typename CoDiReal::Tape;
		using Identifier = typename CoDiReal::Identifier;

    std::vector<int> input_indices;
    std::vector<int> output_indices;

    void print(std::ostream& out) const {
			CoDiReal::getTape().printStatistics(std::cout);
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

		void registerInput(CoDiReal& value) {
			CoDiReal::getTape().registerInput(value);
			input_indices.push_back(value.getIdentifier());
		}

		void registerOutput(CoDiReal& value) {
			CoDiReal::getTape().registerOutput(value);
			output_indices.push_back(value.getIdentifier());
		}

		void start() {
			clear();
			// TODO: Maybe add preallocation of tape.
			CoDiReal::getTape().setActive();
		}

		void stop() {
			CoDiReal::getTape().setPassive();
		}

		void clear() {
			input_indices.clear();
			output_indices.clear();
			CoDiReal::getTape().reset();
		}

		void getFullGradient(double* vec, size_t size) {
			if (size < input_indices.size()){
				_error_("number of requested input values is larger than the number of registered ones. requested_size: " << size << " input_size:" << input_indices.size());
			}

			Tape& tape = CoDiReal::getTape();
			for(size_t i = 0; i < input_indices.size(); i += 1) {
				vec[i] = tape.getGradient(input_indices[i]);
			}
		}

		void updateFullGradient(double* vec, size_t size) {
			if (size < input_indices.size()){
				_error_("number of requested input values is larger than the number of registered ones. requested_size: " << size << " input_size:" << input_indices.size());
			}

			Tape& tape = CoDiReal::getTape();
			for(size_t i = 0; i < input_indices.size(); i += 1) {
				vec[i] += tape.getGradient(input_indices[i]);
			}
		}

		void setGradient(int index, double const& seed) {
			if (output_indices.size() <= index){
				_error_("index value for output is outside bounds of stored output indices. index: " << index << " output_size:" << output_indices.size());
			}
			CoDiReal::getTape().gradient(output_indices[index]) += seed;
		}


		void setFullGradient(double const * vec, size_t size) {
			if (size < output_indices.size()){
				_error_("number of given output values is larger than the number of registered ones. provided_size: " << size << " output_size:" << input_indices.size());
			}

			Tape& tape = CoDiReal::getTape();
			for(size_t i = 0; i < input_indices.size(); i += 1) {
				CoDiReal::getTape().gradient(output_indices[i]) += vec[i];
			}
		}

		void evaluate() {
			CoDiReal::getTape().evaluate();
		}


};

extern CoDi_global codi_global;

#endif /* _HAVE_CODIPACK_ */
#endif /* SRC_C_TOOLKITS_CODIPACK_CODIPACKGLOBAL_H_ */
