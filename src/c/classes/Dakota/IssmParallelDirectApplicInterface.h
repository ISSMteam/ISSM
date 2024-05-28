/*!\file:  IssmParallelDirectApplicInterface.h. This code is only valid for Dakota versions higher than 6!
 *
 * \brief: derived ParallelDirectApplicInterface class declaration and implementation, taylored to ISSM.
 * This class is registered into the interface database of Dakota, and is used to directly call ISSM cores 
 * from Dakota. 
 *
 * This routine helps running ISSM and Dakota in library mode, for Dakota versions that are >=6, and that fully 
 * support parallelism.  The setup is radically different than from version <6! Now, dakota runs the show more. 
 * The reason is that dakota now controls the parallelism in a master/slave setup, and hands over to ISSM a  bunch 
 * of slave communicators, which are then used to run our simulations. Because ISSM is now ESMF compliant, we can 
 * use these communicators to create separate identical FemModel instances on each slave communicator! This allows 
 * us to scale to large jobs (think 1000's of cpus), which we split into multiple sub-slave communicators, which 
 * run the sampling (or forward different, local reliability, optimization you name it) simulations on each slave. 
 * 
 * This is all bootstraped from the main issm_dakota main, (see c/main directory), which is heavily inspired on the
 * main found in the dakota/src/library_mode.cpp code. We also have to create an ISSM code that registers into the 
 * dakota database, which is capable of running ISSM. This is derived from the Dakota class called 
 * ParallelDirectApplicInterface. 
 */ 
#ifndef _ISSMPARALLELDIRECTAPPLICINTERFACE_
#define _ISSMPARALLELDIRECTAPPLICINTERFACE_

/*Issm Configuration: {{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
/*}}}*/

#if !defined(_WRAPPERS_) && defined(_HAVE_DAKOTA_) && _DAKOTA_MAJOR_ >= 6

#include <DirectApplicInterface.hpp>
class FemModel;

namespace SIM {
	class IssmParallelDirectApplicInterface: public Dakota::DirectApplicInterface{

		private: 
			FemModel* femmodel_init;
		public:
			IssmParallelDirectApplicInterface(const Dakota::ProblemDescDB& problem_db, const MPI_Comm& evaluation_comm, int argc, char** argv);
			~IssmParallelDirectApplicInterface();
		protected:
			/// execute an analysis code portion of a direct evaluation invocation
			int derived_map_ac(const Dakota::String& ac_name);
	};
}
/*}}}*/
#endif
#endif
