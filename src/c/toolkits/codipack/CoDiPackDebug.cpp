/*!\file CoDiPackDebug.cpp
 * \brief: implementation specific details for the CoDiPack AD tool.
 */

#include <string>
#include <vector>

#include "CoDiPackDebug.h"
#include "CoDiPackTypes.h"

#include "../mpi/issmmpi.h"

#if defined(_HAVE_ADJOINTPETSC_)
#include <adjoint_petsc/options.h>
#endif


#if defined(_HAVE_CODIPACK_)

using Real = typename CoDiReal::Real;
using Identifier = typename CoDiReal::Identifier;
using Tape = typename CoDiReal::Tape;
using VectorInterface = codi::VectorAccessInterface<Real, Identifier>;
using EventHandle = typename codi::EventSystem<Tape>::Handle;

struct DebugSettings {
		bool outputPrimal           = false;
		bool outputReverse          = true;
		bool outputId               = false;
    bool idFormatLong           = false;

    bool debugEnabled           = false;
    bool dumpTape               = false;
    bool dumpTapePaused         = false;
		int precission              = 12;

    int globalId                = 0;

		std::ostream* stream        = &std::cerr;

    EventHandle dumpEventHandle = {};
};

DebugSettings debugSettings = {};

int CoDiGetUniqueID() {
	debugSettings.globalId += 1;

	return debugSettings.globalId;
}

void writeId(Identifier id) {
  if(debugSettings.outputId) {
    if(debugSettings.idFormatLong) {
      (*debugSettings.stream) << "(" << id << ")";
    }
    else {
      char id_str = 'a';
      if(0 == id) {
        id_str = 'p';
      }
      (*debugSettings.stream) << "(" << id_str << ")";
    }
  }
}

void handleStatementRecord(Tape& tape, Identifier const& lhsIdentifier, Real const& newValue,
                           size_t numActiveVariables, Identifier const* rhsIdentifiers, Real const* jacobians,
                           void* userData) {
  std::ostream& stream = CoDiDebugGetOutputStream();

  if(debugSettings.dumpTape && !debugSettings.dumpTapePaused) {
    stream.setf(std::ios::scientific);
    stream.setf(std::ios::showpos);
    stream.precision(CoDiDebugGetOutputPrecission());

    stream << "primal: " << newValue << "\n";
    // write jacobie values
    for(size_t i = 0; i < numActiveVariables; ++i) {
      stream << "jac: " << i << " " << jacobians[i] << "\n";
    }
  }
}

static void handleStatementEvaluate(Tape& tape, Identifier const& lhsIdentifier, size_t sizeLhsAdjoint,
                                    Real const* lhsAdjoint, void* userData) {
  std::ostream& stream = CoDiDebugGetOutputStream();

  if(debugSettings.dumpTape && !debugSettings.dumpTapePaused) {
    stream.setf(std::ios::scientific);
    stream.setf(std::ios::showpos);
    stream.precision(CoDiDebugGetOutputPrecission());

    for(size_t i = 0; i < sizeLhsAdjoint; ++i) {
      stream << "lhs: " << i << " " << lhsAdjoint[i] << "\n";
    }
  }
}

void CoDiStartDumpTape() {
	if(!CoDiIsDebugOutput()) { return; }

  debugSettings.dumpTape = true;
  debugSettings.dumpEventHandle = codi::EventSystem<Tape>::registerStatementStoreOnTapeListener(handleStatementRecord, nullptr);
}
void CoDiStopDumpTape() {
	if(!CoDiIsDebugOutput()) { return; }

  codi::EventSystemBase<Tape>::deregisterListener(debugSettings.dumpEventHandle);
  debugSettings.dumpTape = false;
}

void CoDiStartDumpEval() {
	if(!CoDiIsDebugOutput()) { return; }

  debugSettings.dumpTape = true;
  debugSettings.dumpEventHandle = codi::EventSystem<Tape>::registerStatementEvaluateListener(handleStatementEvaluate, nullptr);
}

void CoDiStopDumpEval() {
	if(!CoDiIsDebugOutput()) { return; }

  codi::EventSystemBase<Tape>::deregisterListener(debugSettings.dumpEventHandle);
  debugSettings.dumpTape = false;
}

void reverse_dump(Tape* tape, void* d, VectorInterface* vi) {
  debugSettings.dumpTape = !debugSettings.dumpTape;
}

void CoDiPauseDumpTape() {

	if(!CoDiIsDebugOutput()) { return; }
  debugSettings.dumpTapePaused = true;
  CoDiReal::getTape().pushExternalFunction(codi::ExternalFunction<Tape>::create( reverse_dump, nullptr, nullptr));
}
void CoDiResumeDumpTape() {
	if(!CoDiIsDebugOutput()) { return; }

  debugSettings.dumpTapePaused = false;
  CoDiReal::getTape().pushExternalFunction(codi::ExternalFunction<Tape>::create( reverse_dump, nullptr, nullptr));
}

bool          CoDiDebugGetOutputPrimal() {return debugSettings.outputPrimal; }
bool          CoDiDebugGetOutputReverse() {return debugSettings.outputReverse; }
bool          CoDiDebugGetOutputIdentifiers() {return debugSettings.outputId; }
int           CoDiDebugGetOutputPrecission() {return debugSettings.precission; }
std::ostream& CoDiDebugGetOutputStream() {return *debugSettings.stream; }

void CoDiDebugSetOutputPrimal(bool value) {
	debugSettings.outputPrimal = value;
#if defined(_HAVE_ADJOINTPETSC_)
	adjoint_petsc::ADPetscOptionsSetDebugOutputPrimal(value);
#endif
}

void CoDiDebugSetOutputReverse(bool value) {
	debugSettings.outputReverse = value;
#if defined(_HAVE_ADJOINTPETSC_)
	adjoint_petsc::ADPetscOptionsSetDebugOutputReverse(value);
#endif
}

void CoDiDebugSetOutputIdentifiers(bool value) {
	debugSettings.outputId = value;
#if defined(_HAVE_ADJOINTPETSC_)
	adjoint_petsc::ADPetscOptionsSetDebugOutputIdentifiers(value);
#endif
}

void CoDiDebugSetOutputPrecission(int value) {
	debugSettings.precission = value;
#if defined(_HAVE_ADJOINTPETSC_)
	adjoint_petsc::ADPetscOptionsSetDebugOutputPrecission(value);
#endif
}

void CoDiDebugSetOutputStream(std::ostream& value) {
	debugSettings.stream = &value;
#if defined(_HAVE_ADJOINTPETSC_)
	adjoint_petsc::ADPetscOptionsSetDebugOutputStream(value);
#endif
}

bool CoDiIsDebugOutput() {
	return debugSettings.debugEnabled;
}

void CoDiEnableDebugOutput(bool b) {
	if(b) {
		debugSettings.debugEnabled = true;
	}
}
bool CoDiDisableDebugOutput() {
	bool cur = debugSettings.debugEnabled;
	debugSettings.debugEnabled = false;

	return cur;
}

struct MatMatrixEntry {
	int row;
	int col;
	double value;
	int id;
};

inline bool operator<(MatMatrixEntry const& a, MatMatrixEntry const& b) {
	return a.row < b.row || (a.row == b.row && a.col < b.col);
}

struct MatDataEntry {
		int row;
		int col;
		Real value;
		Identifier id;

		MatDataEntry() = default;
		MatDataEntry(int row, int col, Real value, Identifier id) : row(row), col(col), value(value), id(id) {}
		MatDataEntry(int row, int col, Identifier id) : row(row), col(col), id(id) {}
};

inline bool operator<(MatDataEntry const& a, MatDataEntry const& b) {
	return a.row < b.row || (a.row == b.row && a.col < b.col);
}


struct Data_MatDebugOutputReverse {
		int M, N;
		std::vector<MatDataEntry> entries;
		std::string message;
		int id;

		Data_MatDebugOutputReverse(int M, int N, std::string message, int id) : M(M), N(N), entries(), message(message), id(id) {}

		void addEntries(Identifier row, int size, Identifier* cols, CoDiReal* values) {
			size_t start = entries.size();
			for(int i = 0; i < size; i += 1) {
				entries.push_back(MatDataEntry(row, cols[i], values[i].getIdentifier()));
			}

			std::sort(entries.begin() + start, entries.end());
		}

		static void reverse(Tape* tape, void* d, VectorInterface* vi) {
			Data_MatDebugOutputReverse* data = (Data_MatDebugOutputReverse*)d;

			std::ostream& out = *debugSettings.stream;

			int dim = vi->getVectorSize();
			out.setf(std::ios::scientific);
			out.setf(std::ios::showpos);
			out.precision(debugSettings.precission);

			out << data->message << " reverse matrix id: " << data->id << std::endl;
			for(int cur_dim = 0; cur_dim < dim; cur_dim += 1) {
				int my_rank=IssmComm::GetRank();
				for(int i=0;i<IssmComm::GetSize();i++){
					if(my_rank==i){
						if(i==0) {
							out << "Matrix of global size " << data->M << "x" << data->N << std::endl;
						}
						out << "Rank: " << my_rank  << " dim: " << cur_dim <<"\n";
						for(size_t j = 0; j < data->entries.size(); j += 1) {
							MatDataEntry& entry = data->entries[j];
							out << entry.row << " " << entry.col << " " << vi->getAdjoint(entry.id, cur_dim);
              writeId(entry.id);
							out << "\n";
						}
					}
					out.flush();
					ISSM_MPI_Barrier(IssmComm::GetComm());
				}
			}
		}

		static void free(Tape* tape, void* d) {
			Data_MatDebugOutputReverse* data = (Data_MatDebugOutputReverse*)d;

			delete data;
		}

		void push() {
			CoDiReal::getTape().pushExternalFunction(codi::ExternalFunction<Tape>::create(&reverse, this, &free));
		}
};

struct Data_MatDebugOutput {
		int id;

		std::vector<MatDataEntry> entries;

		Data_MatDebugOutputReverse* reverse;

		Data_MatDebugOutput() = default;
};

void* MatDebugOutputStart(std::string message, int M, int N) {
	if(!CoDiIsDebugOutput()) { return nullptr; }
	if(!(debugSettings.outputPrimal || debugSettings.outputReverse)) {
		return nullptr;
	}

	int id = CoDiGetUniqueID();

	Data_MatDebugOutput* handle = new Data_MatDebugOutput();

	if(debugSettings.outputReverse) {
		handle->reverse = new Data_MatDebugOutputReverse(M, N, message, id);
	}

	std::ostream& out = *debugSettings.stream;

	if(debugSettings.outputPrimal) {
	out.setf(std::ios::scientific);
	out.setf(std::ios::showpos);
	out.precision(debugSettings.precission);
	if(0 == IssmComm::GetRank()) {
		out << message << " forward matrix id: " << id << std::endl;
		out << "Matrix of global size " << M << "x" << N << std::endl;
	}
	ISSM_MPI_Barrier(IssmComm::GetComm());
  }

  return handle;
}

void MatDebugOutputAddRow(void* h, int row, int size, int* cols, CoDiReal* values) {
	if(nullptr == h) { return; }

	Data_MatDebugOutput* handle = (Data_MatDebugOutput*)h;

	if(debugSettings.outputPrimal) {
		for(int i = 0; i < size; i += 1) {
			handle->entries.push_back(MatDataEntry(row, cols[i], values[i].getValue(), values[i].getIdentifier()));
		}
	}
	if(debugSettings.outputReverse) {
		handle->reverse->addEntries(row, size, cols, values);
	}
}
void MatDebugOutputAddRow(void*, int, int, int*, double*) {}

void MatDebugOutputFinish(void* h) {
	if(nullptr == h) { return; }

	Data_MatDebugOutput* handle = (Data_MatDebugOutput*)h;

	if(debugSettings.outputPrimal) {
		std::sort(handle->entries.begin(), handle->entries.end());
		std::ostream& out = *debugSettings.stream;


		for(int cur_rank = 0; cur_rank < IssmComm::GetSize(); cur_rank += 1) {
			if(cur_rank == IssmComm::GetRank()) {
				out << "Rank: " << cur_rank << "\n";
				for(size_t i = 0; i < handle->entries.size(); i += 1) {
					MatDataEntry& entry = handle->entries[i];
					out << entry.row << " " << entry.col << " " << entry.value;
          writeId(entry.id);
					out << "\n";
				}
				out.flush();
			}
			ISSM_MPI_Barrier(IssmComm::GetComm());
		}
	}
	if(debugSettings.outputReverse) {
		handle->reverse->push();
	}

	delete handle;
}

struct Data_VecDebugOutputReverse {
		int M;
		std::vector<Identifier> vec_i;
		std::string message;
		int id;

		bool is_array;

		Data_VecDebugOutputReverse(int M, int size, std::string message, int id, bool is_array) : M(M), vec_i(size), message(message), id(id), is_array(is_array) {}

		void getIdentifiers(CoDiReal* values) {
			for(size_t i = 0; i < vec_i.size(); i += 1) {
				vec_i[i] = values[i].getIdentifier();
			}
		}

		static void reverse(Tape* tape, void* d, VectorInterface* vi) {
			Data_VecDebugOutputReverse* data = (Data_VecDebugOutputReverse*)d;

			int dim = vi->getVectorSize();

			std::ostream& out = *debugSettings.stream;

			out.setf(std::ios::scientific);
			out.setf(std::ios::showpos);
			out.precision(debugSettings.precission);

			int my_rank=IssmComm::GetRank();
			if(my_rank == 0) {
				if(data->is_array) {
					out << data->message << " reverse array id: " << data->id << std::endl;
					out << "Array of size " << data->M << std::endl;
				} else {
					out << data->message << " reverse vector id: " << data->id << std::endl;
					out << "Vector of global size M=" << data->M << std::endl;
				}
			}
			for(int cur_dim = 0; cur_dim < dim; cur_dim += 1) {
				for(int i=0;i<IssmComm::GetSize();i++){
					if(my_rank==i){
						if(i==0) {
							out << "Rank: " << my_rank << " dim: " << cur_dim <<"\n";
						}
						for(size_t j=0;j<data->vec_i.size();j++){
							out << vi->getAdjoint(data->vec_i[j], cur_dim);
              writeId(data->vec_i[j]);
							out << "\n";
						}
					}
					out.flush();
					ISSM_MPI_Barrier(IssmComm::GetComm());
				}
			}
		}

		static void free(Tape* tape, void* d) {
			Data_VecDebugOutputReverse* data = (Data_VecDebugOutputReverse*)d;

			delete data;
		}

		void push() {
			CoDiReal::getTape().pushExternalFunction(codi::ExternalFunction<Tape>::create(&reverse, this, &free));
		}
};

void VecDebugOutputImpl(std::string message, int M, int m, CoDiReal* values, bool is_array) {
	if(!CoDiIsDebugOutput()) { return; }
	if(!(debugSettings.outputPrimal || debugSettings.outputReverse)) {
		return;
	}

	int id = CoDiGetUniqueID();

	if(debugSettings.outputPrimal) {
		std::ostream& out = *debugSettings.stream;

		out.setf(std::ios::scientific);
		out.setf(std::ios::showpos);
		out.precision(debugSettings.precission);
		if(0 == IssmComm::GetRank()) {
			if(is_array) {
				out << message << " forward array id: " << id << std::endl;
				out << "Array of size " << M << std::endl;
			} else {
				out << message << " forward vector id: " << id << std::endl;
				out << "Vector of global size M=" << M << std::endl;
			}

			for(int cur_rank = 0; cur_rank < IssmComm::GetSize(); cur_rank += 1) {
				if(cur_rank == IssmComm::GetRank()) {
					out << "Rank: " << cur_rank << "\n";
					for(size_t i = 0; i < m; i += 1) {
						out << values[i].getValue();
            writeId(values[i].getIdentifier());
						out << "\n";
					}
					out.flush();
				}
				ISSM_MPI_Barrier(IssmComm::GetComm());
			}
		}
		ISSM_MPI_Barrier(IssmComm::GetComm());
	}

	if(debugSettings.outputReverse) {
		Data_VecDebugOutputReverse* data = new Data_VecDebugOutputReverse(M, m, message, id, is_array);
		data->getIdentifiers(values);
		data->push();
	}
}

void VecDebugOutput(std::string message, int M, int m, CoDiReal* values) {
	VecDebugOutputImpl(message, m, m, values, false);
}
void VecDebugOutput(std::string, int, int, double*) {}

void ArrayDebugOutput(std::string message, int m, CoDiReal* values) {
	VecDebugOutputImpl(message, m, m, values, true);
}
void ArrayDebugOutput(std::string, int, double*) {}

#endif /* _HAVE_CODIPACK_ */
