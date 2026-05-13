/*!\file StressBalanceEmulatorx
 * \brief: evaluate stressbalance emulator and update velocity inputs
 */

#include "./StressBalanceEmulatorx.h"
#include "../../toolkits/toolkits.h"
#include "../modules.h"

#include <cmath>
#include <sstream>
#include <sys/types.h>
#include <vector>

#include <pybind11/numpy.h>
namespace py = pybind11;

void StressBalanceEmulatorx(FemModel* femmodel){/*{{{*/

	if(!femmodel->parameters->Exist(StressbalanceEmulatorEnum)) _error_("StressbalanceEmulatorEnum not found; stressbalance emulator was not initialized");
	Param* emulator_param = femmodel->parameters->FindParamObject(StressbalanceEmulatorEnum);
	if(emulator_param->ObjectEnum()!=EmulatorParamEnum) _error_("Parameter should be EmulatorParam");
	EmulatorParam* emulator = xDynamicCast<EmulatorParam*>(emulator_param);

	int rank             = IssmComm::GetRank();
	int numberofvertices = femmodel->vertices->NumberOfVertices();
	IssmDouble yts;
	femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);

	/*Recover geometry and model inputs*/
	IssmDouble* thickness = NULL;
	IssmDouble* base      = NULL;
	IssmDouble* frictionc = NULL;
	IssmDouble* mask      = NULL;
	IssmDouble* vx        = NULL;
	IssmDouble* vy        = NULL;
	int thickness_size = 0;
	int base_size      = 0;
	int frictionc_size = 0;
	int mask_size      = 0;
	int vx_size        = 0;
	int vy_size        = 0;

	GetVectorFromInputsx(&thickness,&thickness_size,femmodel,ThicknessEnum);
	GetVectorFromInputsx(&base,&base_size,femmodel,BaseEnum);
	GetVectorFromInputsx(&frictionc,&frictionc_size,femmodel,FrictionCEnum);
	GetVectorFromInputsx(&mask,&mask_size,femmodel,MaskIceLevelsetEnum);
	GetVectorFromInputsx(&vx,&vx_size,femmodel,VxEnum);
	GetVectorFromInputsx(&vy,&vy_size,femmodel,VyEnum);

	if(thickness_size!=numberofvertices || base_size!=numberofvertices || frictionc_size!=numberofvertices || mask_size!=numberofvertices || vx_size!=numberofvertices || vy_size!=numberofvertices){
		_error_("StressBalanceEmulatorx expects vertex-based Thickness, Base, FrictionC, MaskIceLevelset, Vx, and Vy inputs");
	}

	std::vector<IssmDouble> vx_all(numberofvertices,0.);
	std::vector<IssmDouble> vy_all(numberofvertices,0.);
	std::vector<IssmDouble> vel_all(numberofvertices,0.);

	/*Evaluate emulator on rank 0*/
	if(rank==0){
		const size_t nvertices = static_cast<size_t>(numberofvertices);
		const size_t nfeatures = 6;
		std::vector<float> features(nvertices*nfeatures);

		for(size_t i=0;i<nvertices;i++){
			features[i*nfeatures + 0] = static_cast<float>(mask[i]/400000.0);
			features[i*nfeatures + 1] = static_cast<float>(thickness[i]/5000.0);
			features[i*nfeatures + 2] = static_cast<float>((base[i]-1200.0)/5000.0);
			features[i*nfeatures + 3] = static_cast<float>(frictionc[i]/12000.0);
			features[i*nfeatures + 4] = static_cast<float>(vx[i]*yts/1000.0);
			features[i*nfeatures + 5] = static_cast<float>(vy[i]*yts/1000.0);
		}

		try{
			std::vector<ssize_t> features_shape = {static_cast<ssize_t>(nvertices),static_cast<ssize_t>(nfeatures)};
			py::array_t<float> features_np(features_shape,features.data());
			py::object output_obj = emulator->mod.attr("predict_logits_np")(features_np);
			py::array_t<float,py::array::c_style | py::array::forcecast> output_array(output_obj);
			py::buffer_info output_info = output_array.request();

			if(output_info.ndim!=2){
				_error_("StressBalanceEmulatorx output is not a 2D array");
			}
			size_t rows = static_cast<size_t>(output_info.shape[0]);
			size_t cols = static_cast<size_t>(output_info.shape[1]);
			if(rows!=nvertices || cols!=2){
				std::ostringstream message;
				message << "StressBalanceEmulatorx unexpected output shape: (" << rows << ", " << cols << "), expected (" << nvertices << ", 2)";
				_error_(message.str());
			}

			const float* output = static_cast<const float*>(output_info.ptr);
			for(size_t i=0;i<nvertices;i++){
				vx_all[i] = static_cast<IssmDouble>(output[i*2 + 0]*1000.0/yts);
				vy_all[i] = static_cast<IssmDouble>(output[i*2 + 1]*1000.0/yts);
			}
		}
		catch(const py::error_already_set& e){
			_error_("StressBalanceEmulatorx Python exception: " << e.what());
		}
	}

	xDelete<IssmDouble>(thickness);
	xDelete<IssmDouble>(base);
	xDelete<IssmDouble>(frictionc);
	xDelete<IssmDouble>(mask);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(vy);

	/*Broadcast emulator output and update model inputs*/
	if(numberofvertices>0){
		ISSM_MPI_Comm comm = IssmComm::GetComm();
		ISSM_MPI_Bcast(vx_all.data(),numberofvertices,ISSM_MPI_DOUBLE,0,comm);
		ISSM_MPI_Bcast(vy_all.data(),numberofvertices,ISSM_MPI_DOUBLE,0,comm);

		for(int i=0;i<numberofvertices;i++){
			vel_all[i] = sqrt(vx_all[i]*vx_all[i] + vy_all[i]*vy_all[i]);
		}

		InputUpdateFromVectorx(femmodel,vx_all.data(),VxEnum,VertexSIdEnum);
		InputUpdateFromVectorx(femmodel,vy_all.data(),VyEnum,VertexSIdEnum);
		InputUpdateFromVectorx(femmodel,vel_all.data(),VelEnum,VertexSIdEnum);

		std::vector<IssmDouble> pressure(numberofvertices,0.);
		InputUpdateFromVectorx(femmodel,pressure.data(),PressureEnum,VertexSIdEnum);
	}

	return;
}/*}}}*/
