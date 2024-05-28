/*!\file:  IssmDirectApplicInterface.* This code is only valid for Dakota versions lower than 5!
 *
 * \brief: derived DirectApplicInterface class declaration and implementation, taylored to the ISSM cores. 
 * This class is registered into the interface database of Dakota, and is used to directly call ISSM cores 
 * from Dakota. 
 *
 * This routine helps running ISSM and Dakota in library mode, for Dakota versions that are <=5, and which 
 * do not really support parallelism. In library mode, Dakota does not 
 * run as an execuatble. Its capabilities are linked into the ISSM software, and ISSM calls dakota routines 
 * directly from the dakota library. dakota_core.cpp is the code that is in charge of calling those routines. 
 *
 * Prior to versions 6 and more, Dakota had its own way of running in parallel (for embarassingly parallel jobs). 
 * We do not want that, as ISSM knows exactly how to run "really parallel" jobs that use all CPUS. To bypass Dakota's parallelism, 
 * we overloaded the constructor for the parallel library (see the Dakota patch in the externalpackages/dakota
 * directory). This overloaded constructor fires up Dakota serially on CPU 0 only! We take care of broadcasting 
 * to the other CPUS, hence ISSM is running in parallel, and Dakota serially on CPU0. 
 *
 * Now, how does CPU 0 drive all other CPUS to carry out sensitivity analysese? By synchronizing its call to 
 * our ISSM cores (stressbalance_core, thermal_core, transient_core, etc ...) on CPU 0 with all other CPUS. 
 * This explains the structure of dakota_core.cpp, where cpu 0 runs Dakota, the Dakota pluggin fires up DakotaSpawnCore.cpp, 
 * while the other CPUS are waiting for a broadcast from CPU0, once they get it, they also fire up 
 * DakotaSpawnCore. In the end, DakotaSpawnCore is fired up on all CPUS, with CPU0 having Dakota inputs, that it will 
 * broacast to other CPUS. 
 *
 * Now, how does dakota call the DakotaSpawnCore routine? The DakotaSpawnCore is embedded into the IssmDirectApplicInterface object 
 * which is derived from the Direct Interface Dakota objct. This is the only way to run Dakota in library 
 * mode (see their developper guide for more info). Dakota registers the IssmDirectApplicInterface object into its own 
 * database, and calls on the embedded DakotaSpawnCore from CPU0. 
 *
 */ 

/*Issm Configuration: {{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif 
/*}}}*/

#if !defined(_WRAPPERS_) && defined(_HAVE_DAKOTA_) && (_DAKOTA_MAJOR_ <= 5) //this only works for Dakota <=5, which had no effective parallel capabilities yet.

/*Dakota include files:{{{*/
#if (_DAKOTA_MAJOR_ < 5 || (_DAKOTA_MAJOR_ == 5 && _DAKOTA_MINOR_ < 3))
#include <DirectApplicInterface.H>
#include <DakotaResponse.H>
#include <ParamResponsePair.H>
#include <system_defs.h>
#include <ProblemDescDB.H>
#include <ParallelLibrary.H>
#else
#include <DirectApplicInterface.hpp>
#include <DakotaResponse.hpp>
#include <ParamResponsePair.hpp>
#include <ProblemDescDB.hpp>
#include <ParallelLibrary.hpp>
#endif
/*}}}*/

int  DakotaSpawnCore(double* d_responses, int d_numresponses, double* d_variables, char** d_variables_descriptors,int d_numvariables, void* void_femmodel,int counter);

/*IssmDirectApplicInterface class */
namespace SIM {
	class IssmDirectApplicInterface: public Dakota::DirectApplicInterface{
		public:
			/*these fields are used by core solutions: */
			void *femmodel;
			int   counter;
			/*Constructors/Destructors*/
			IssmDirectApplicInterface(const Dakota::ProblemDescDB& problem_db,void* in_femmodel):Dakota::DirectApplicInterface(problem_db){/*{{{*/
				femmodel = in_femmodel;
				counter  = 0;
			}/*}}}*/
			~IssmDirectApplicInterface(){/*{{{*/
				/* Virtual destructor handles referenceCount at Interface level. */ 
			}/*}}}*/
		protected:
			/*execute the input filter portion of a direct evaluation invocation*/
			//int derived_map_if(const Dakota::String& if_name);
			/*execute an analysis code portion of a direct evaluation invocation*/
			int derived_map_ac(const Dakota::String& driver){/*{{{*/

				int i;
				IssmDouble* variables=NULL;
				char** variable_descriptors=NULL;
				char*  variable_descriptor=NULL;
				IssmDouble* responses=NULL;

				/*increae counter: */
				counter++;

				/*Before launching analysis, we need to transfer the dakota inputs into Issm 
				 *readable variables: */

				/*First, the variables: */
				variables=xNew<IssmDouble>(numACV);
				for(i=0;i<numACV;i++){
					variables[i]=xC[i];
				}
				/*The descriptors: */
				variable_descriptors=xNew<char*>(numACV);
				for(i=0;i<numACV;i++){
					std::string label=xCLabels[i];
					variable_descriptor=xNew<char>(strlen(label.c_str())+1);
					memcpy(variable_descriptor,label.c_str(),(strlen(label.c_str())+1)*sizeof(char));

					variable_descriptors[i]=variable_descriptor;
				}

				/*Initialize responses: */
				responses=xNewZeroInit<IssmDouble>(numFns);

				/*run core solution: */
				DakotaSpawnCore(responses,numFns, variables,variable_descriptors,numACV,femmodel,counter);

				/*populate responses: */
				for(i=0;i<numFns;i++){
					fnVals[i]=responses[i];
				}

				/*Free resources:*/
				xDelete<IssmDouble>(variables);
				for(i=0;i<numACV;i++){
					variable_descriptor=variable_descriptors[i];
					xDelete<char>(variable_descriptor);
				}
				xDelete<char*>(variable_descriptors);
				xDelete<IssmDouble>(responses);

				return 0;
			}/*}}}*/
			/*execute the output filter portion of a direct evaluation invocation*/
			//int derived_map_of(const Dakota::String& of_name);
			/*add for issm: */
			int GetCounter(){/*{{{*/
				return counter;
			}/*}}}*/
		private:
	};
} 
#endif
