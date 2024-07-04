/*!\file CoDiPackGlobal.cpp
 * \brief: implementation specific details for the CoDiPack AD tool.
 */

#include "CoDiPackGlobal.h"

#if defined(_HAVE_CODIPACK_)

#include "../../classes/Params/Parameters.h"
#include "../../shared/Exceptions/exceptions.h"

void CoDi_global::init(Parameters* parameters) {
	run_count = 0;

	parameters->FindParam(&has_time_output, AutodiffOutputTimeEnum);
	parameters->FindParam(&has_memory_output, AutodiffOutputTapeMemoryEnum);
	parameters->FindParam(&has_preaccumulation, AutodiffEnablePreaccumulationEnum);

	if (0 == IssmComm::GetRank()) {
		if (has_time_output) {
			time_output.open("ad_time.dat");
			time_output << "Run; Record; Reverse; Total; Tape Memory\n";
			time_output << std::fixed << std::setprecision(3);
		}

		if (has_memory_output) {
			memory_output.open("ad_mem.dat");
			memory_output << "Run; ";
			CoDiReal::getTape().getTapeValues().formatHeader(memory_output);
		}
	}
}

void CoDi_global::print(std::ostream& out) const {
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

void CoDi_global::registerInput(CoDiReal& value) {
	CoDiReal::getTape().registerInput(value);
	input_indices.push_back(value.getIdentifier());
}

void CoDi_global::registerOutput(CoDiReal& value) {
	CoDiReal::getTape().registerOutput(value);
	output_indices.push_back(value.getIdentifier());
}

void CoDi_global::start() {
	clear();
	// TODO: Maybe add preallocation of tape.
	CoDiReal::getTape().setActive();

	recordStart();
}

void CoDi_global::stop() {
	recordEnd();
	CoDiReal::getTape().setPassive();
}

void CoDi_global::clear() {
	input_indices.clear();
	output_indices.clear();
	CoDiReal::getTape().reset();
}

void CoDi_global::getFullGradient(double* vec, size_t size) {
	if (size < input_indices.size()){
		_error_("number of requested input values is larger than the number of registered ones. requested_size: " << size << " input_size:" << input_indices.size());
	}

	Tape& tape = CoDiReal::getTape();
	for(size_t i = 0; i < input_indices.size(); i += 1) {
		vec[i] = tape.getGradient(input_indices[i]);
	}
}

void CoDi_global::updateFullGradient(double* vec, size_t size) {
	if (size < input_indices.size()){
		_error_("number of requested input values is larger than the number of registered ones. requested_size: " << size << " input_size:" << input_indices.size());
	}

	Tape& tape = CoDiReal::getTape();
	for(size_t i = 0; i < input_indices.size(); i += 1) {
		vec[i] += tape.getGradient(input_indices[i]);
	}
}

void CoDi_global::setGradient(int index, double const& seed) {
	if (output_indices.size() <= index){
		_error_("index value for output is outside bounds of stored output indices. index: " << index << " output_size:" << output_indices.size());
	}
	CoDiReal::getTape().gradient(output_indices[index]) += seed;
}


void CoDi_global::setFullGradient(double const * vec, size_t size) {
	if (size < output_indices.size()){
		_error_("number of given output values is larger than the number of registered ones. provided_size: " << size << " output_size:" << input_indices.size());
	}

	Tape& tape = CoDiReal::getTape();
	for(size_t i = 0; i < input_indices.size(); i += 1) {
		CoDiReal::getTape().gradient(output_indices[i]) += vec[i];
	}
}

void CoDi_global::evaluate() {
	evaluateStart();
	CoDiReal::getTape().evaluate();
	evaluateEnd();

	outputTimeAndMem();
}

void CoDi_global::recordStart() {
	if (has_time_output)  {
		record_start = std::chrono::system_clock::now();
	}

	if (has_time_output || has_memory_output) {
		run_count += 1;
	}
}

void CoDi_global::recordEnd() {
	if (has_time_output)  {
		record_end = std::chrono::system_clock::now();
	}

}

void CoDi_global::evaluateStart() {
	if (has_time_output)  {
		evaluate_start = std::chrono::system_clock::now();
	}
}

void CoDi_global::evaluateEnd() {
	if (has_time_output)  {
		evaluate_end = std::chrono::system_clock::now();
	}

}

void CoDi_global::outputTimeAndMem() {

	if (has_time_output || has_memory_output) {
		codi::TapeValues tapeValues = CoDiReal::getTape().getTapeValues();
		tapeValues.combineData();


		if (has_time_output) {
			IssmPDouble times[2];
			times[0] = ((std::chrono::duration<double>) (record_end - record_start)).count();
			times[1] = ((std::chrono::duration<double>) (evaluate_end - evaluate_start)).count();
			ISSM_MPI_Allreduce(ISSM_MPI_IN_PLACE, times,2,ISSM_MPI_PDOUBLE,ISSM_MPI_MAX,IssmComm::GetComm());

			if (0 == IssmComm::GetRank()) {
				time_output << run_count << "; " << times[0] << "; " << times[1] << "; "
							<< (times[0] + times[1]) << "; " << (tapeValues.getUsedMemorySize() / 1024.0/1024.0) << std::endl;
			}
		}

		if (has_memory_output) {
			if (0 == IssmComm::GetRank()) {
				memory_output << run_count << "; ";
				tapeValues.formatRow(memory_output);
			}
		}
	}
}



CoDi_global codi_global = {};

#endif /* _HAVE_CODIPACK_ */
