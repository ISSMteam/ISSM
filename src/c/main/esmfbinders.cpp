/*!\file:  esmfbinder.cpp
 * \brief: ESMF binders for ISSM. Binders developed initially for the GEOS-5 framework.
 */ 

#include "./issm.h"

/*GEOS 5 specific declarations:*/
const int GCMForcingNumTerms = 1;
const int GCMForcingTerms[GCMForcingNumTerms]= { SMBgcmEnum}; 
const int ISSMOutputNumTerms = 1;
const int ISSMOutputTerms[ISSMOutputNumTerms]= { SurfaceEnum };

extern "C" {

	FemModel *femmodel;

	/*GEOS 5*/
      void InitializeISSM(int argc, char** argv, int* pnumberofelements, int* pnumberofnodes, MPI_Fint* Fcomm){ /*{{{*/
		int numberofelements;
        int numberofnodes;
          
        /* convert Fortran MPI comm to C MPI comm */
        MPI_Comm Ccomm = MPI_Comm_f2c(*Fcomm);             
                
        /*Initialize femmodel from arguments provided command line: */
		femmodel = new FemModel(argc,argv,Ccomm);

		/*Get number of nodes and elements from the mesh: */
		numberofelements=femmodel->elements->NumberOfElements();
        numberofnodes=femmodel->vertices->Size();

		/*Bypass SMB model, will be provided by GCM! */
		femmodel->parameters->SetParam(SMBgcmEnum,SmbEnum); 

        /*Restart file: */
		femmodel->Restart();

		/*Assign output pointers: */
		*pnumberofelements=numberofelements;
        *pnumberofnodes=numberofnodes;
	} /*}}}*/

	void RunISSM(IssmDouble dt, IssmDouble* gcmforcings, IssmDouble* issmoutputs){ /*{{{*/

		int numberofelements;
		IssmDouble yts;
		IssmDouble rho_ice;
		IssmDouble area;
		IssmDouble start_time,final_time;

		/*Figure out number of elements: */
		numberofelements=femmodel->elements->Size();

		/*Fetch some necessary constants: */
		femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);

		/*Setup gcm forcings as element-wise input: {{{ */
		for (int f=0;f<GCMForcingNumTerms;f++){

			int forcing_type=GCMForcingTerms[f];

			for (int i=0;i<femmodel->elements->Size();i++){
				Element* element=dynamic_cast<Element*>(femmodel->elements->GetObjectByOffset(i));

				switch(forcing_type){
					case SMBgcmEnum:
						/*{{{*/
						{
						/*Recover rho_ice: */
						rho_ice=element->FindParam(MaterialsRhoIceEnum);

						/*Recover area of element: */
						area=element->SurfaceArea();

						/*Recover smb forcing from the gcm forcings: */
						IssmDouble smbforcing=*(gcmforcings+f*numberofelements+i); 

						/*Convert to SI. The smbforcing from GEOS-5 in kg/s, and we transform it into m/s: */
						smbforcing=smbforcing/(rho_ice*area);

						/*Add into the element as new forcing :*/
						element->AddInput(SmbMassBalanceEnum,&smbforcing,P0Enum);
						}
						/*}}}*/
						break; 
					default: 
						{ _error_("Unknown forcing type " << forcing_type << "\n"); }
						break;
				}
			}
		}

		/*}}}*/

		/*Retrieve ISSM outputs and pass them back to the Gcm : {{{*/
		for (int f=0;f<ISSMOutputNumTerms;f++){

			int output_type=ISSMOutputTerms[f];

			for (int i=0;i<femmodel->elements->Size();i++){
				Element* element=dynamic_cast<Element*>(femmodel->elements->GetObjectByOffset(i));

				switch(output_type){
					case SurfaceEnum:
						/*{{{*/
						{

						IssmDouble surface;

						/*Recover surface from the ISSM element: */
						Input* surface_input = element->GetInput(SurfaceEnum); _assert_(surface_input);
						surface_input->GetInputAverage(&surface);

						*(issmoutputs+f*numberofelements+i) = surface;

						}
						/*}}}*/
						break; 
					default: 
						{ _error_("Unknown output type " << output_type << "\n"); }
						break;
				}
			}
		}

		/*}}}*/

		/*Before running, setup the time interval: */
		femmodel->parameters->FindParam(&start_time,TimeEnum);
		final_time=start_time+dt;
		femmodel->parameters->SetParam(final_time,TimesteppingFinalTimeEnum); //we are bypassing ISSM's initial final time!

		/*Now, run: */
		femmodel->Solve();

		/*For the next time around, save the final time as start time */
		femmodel->parameters->SetParam(final_time,TimesteppingStartTimeEnum);
	} /*}}}*/

	void FinalizeISSM(){ /*{{{*/

		/*Output results: */
		OutputResultsx(femmodel);

		/*Check point: */
		femmodel->CheckPoint();

		/*Wrap up: */
		delete femmodel; femmodel=NULL;
	} /*}}}*/

    void GetNodesISSM(int* nodeIds, IssmDouble* nodeCoords){ 
        /*obtain nodes of mesh for creating ESMF version in Fortran interface*/
        /*nodeIds are the global Id's of the nodes and nodeCoords are the    */
        /*(x,y) coordinates, as described in the ESMF reference document     */
        int sdim = 2; // spatial dimension (2D map-plane)
        for (int i=0;i<femmodel->vertices->Size();i++){
            Vertex* vertex = xDynamicCast<Vertex*>(femmodel->vertices->GetObjectByOffset(i));
            *(nodeIds+i)     = vertex->Sid()+1;
            *(nodeCoords+sdim*i+0) = vertex->x;
            *(nodeCoords+sdim*i+1) = vertex->y;
        }
    }

    void GetElementsISSM(int* elementIds, int* elementConn){
        /*obtain elements of mesh for creating ESMF version in Fortran interface*/
        /*Element connectivity (elementConn) contains the indices of the nodes  */
        /*that form the element as described in the ESMF reference document     */ 
        for(int i=0;i<femmodel->elements->Size();i++){
            Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
            *(elementIds + i)    = element->Sid()+1;
            *(elementConn + i*3+0) = element->vertices[0]->Lid()+1;
            *(elementConn + i*3+1) = element->vertices[1]->Lid()+1;
            *(elementConn + i*3+2) = element->vertices[2]->Lid()+1;
        }
    }

}
