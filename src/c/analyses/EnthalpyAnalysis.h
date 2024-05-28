/*! \file EnthalpyAnalysis.h 
 *  \brief: header file for generic external result object
 */

#ifndef _EnthalpyAnalysis_
#define _EnthalpyAnalysis_

/*Headers*/
#include "./Analysis.h"
#include "../classes/classes.h"

class EnthalpyAnalysis: public Analysis{

	public:
		/*Model processing*/
		void CreateConstraints(Constraints* constraints,IoModel* iomodel);
		void CreateLoads(Loads* loads, IoModel* iomodel);
		void CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr=false);
		int  DofsPerNode(int** doflist,int domaintype,int approximation);
		void UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type);
		void UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum);

		/*Finite element Analysis*/
		static void       ApplyBasalConstraints(IssmDouble* serial_spc,Element* element);
		static void       ComputeBasalMeltingrate(FemModel* femmodel);
		static void       ComputeBasalMeltingrate(Element* element);
		void              Core(FemModel* femmodel);
		void              PreCore(FemModel* femmodel);
		ElementVector*    CreateDVector(Element* element);
		ElementMatrix*    CreateJacobianMatrix(Element* element);
		ElementMatrix*    CreateKMatrix(Element* element);
		ElementMatrix*    CreateKMatrixVolume(Element* element);
		ElementMatrix*    CreateKMatrixShelf(Element* element);
		ElementVector*    CreatePVector(Element* element);
		ElementVector*    CreatePVectorVolume(Element* element);
		ElementVector*    CreatePVectorSheet(Element* element);
		ElementVector*    CreatePVectorShelf(Element* element);
		static void       DrainWaterfraction(FemModel* femmodel);
 		static void       ComputeWaterfractionDrainage(FemModel* femmodel);
 		static void       DrainageUpdateWatercolumn(FemModel* femmodel);
 		static void       DrainageUpdateEnthalpy(FemModel* femmodel);
		static IssmDouble EnthalpyDiffusionParameter(Element* element,IssmDouble enthalpy,IssmDouble pressure);
		static IssmDouble EnthalpyDiffusionParameterVolume(Element* element,int enthalpy_enum);
		static void       GetBasalConstraints(Vector<IssmDouble>* vec_spc,Element* element);
		static void       GetBasalConstraintsSteadystate(Vector<IssmDouble>* vec_spc,Element* element);
		static void       GetBasalConstraintsTransient(Vector<IssmDouble>* vec_spc,Element* element);
		void              GetBConduct(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss);
		void              GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element);
		static int        GetThermalBasalCondition(Element* element, IssmDouble enthalpy, IssmDouble enthalpy_up, IssmDouble pressure, IssmDouble pressure_up, IssmDouble watercolumn, IssmDouble meltingrate);
		static IssmDouble GetWetIceConductivity(Element* element, IssmDouble enthalpy, IssmDouble pressure);
		void              GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index);
		void              InputUpdateFromSolution(IssmDouble* solution,Element* element);
		static void       PostProcessing(FemModel* femmodel);
		static IssmDouble PureIceEnthalpy(Element* element,IssmDouble pressure);
		static IssmDouble TMeltingPoint(Element* element,IssmDouble pressure);
		static void       UpdateBasalConstraints(FemModel* femmodel);
		void              UpdateConstraints(FemModel* femmodel);
};
#endif
