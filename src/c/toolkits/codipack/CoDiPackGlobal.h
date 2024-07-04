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
#include <fstream>
#include <iostream>
#include <chrono>

#include "CoDiPackTypes.h"

class Parameters;

struct CoDi_global {

		using Tape = typename CoDiReal::Tape;
		using Identifier = typename CoDiReal::Identifier;

		using time_point = typename std::chrono::system_clock::time_point;

		std::vector<int> input_indices;
		std::vector<int> output_indices;

		bool has_time_output;
		bool has_memory_output;

		std::ofstream time_output;
		std::ofstream memory_output;

		time_point record_start;
		time_point record_end;

		time_point evaluate_start;
		time_point evaluate_end;

		int run_count;

		// Misc functions.

		void init(Parameters* parameters);
		void print(std::ostream& out) const;

		// Tape management and input/output registration.
		void registerInput(CoDiReal& value);
		void registerOutput(CoDiReal& value);
		void start();
		void stop();
		void clear();

		// Gradient setters, getters and tape evaluation.

		void getFullGradient(double* vec, size_t size);
		void updateFullGradient(double* vec, size_t size);
		void setGradient(int index, double const& seed);
		void setFullGradient(double const * vec, size_t size);
		void evaluate();

	private:
		// Time and memory measurement.
		void recordStart();
		void recordEnd();

		void evaluateStart();
		void evaluateEnd();

		void outputTimeAndMem();
};

extern CoDi_global codi_global;

#endif /* _HAVE_CODIPACK_ */
#endif /* SRC_C_TOOLKITS_CODIPACK_CODIPACKGLOBAL_H_ */
