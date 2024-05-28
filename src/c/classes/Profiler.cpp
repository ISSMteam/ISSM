/*!\file Profiler.c
 * \brief: implementation of the Profiler object
 */

/*Include files*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Profiler.h"
#include "../toolkits/toolkits.h"

/*Profiler constructors and destructors:*/
Profiler::Profiler(){/*{{{*/
	for(int i=0;i<MAXPROFSIZE;i++){
		this->time[i]          = 0.;
		this->time_start[i]    = 0.;
		this->flops[i]         = 0.;
		this->flops_start[i]   = 0.;
		this->memory[i]        = 0.;
		this->memory_start[i]  = 0.;
		this->running[i]       = false;
		this->used[i]          = false;
	}
} /*}}}*/
Profiler::~Profiler(){/*{{{*/
	/*Nothing to delete, everything is statically allocated*/
} /*}}}*/
Object* Profiler::copy(){/*{{{*/
	/*First do simple copy: */
	Profiler* output=new Profiler();

	for(int i=0;i<MAXPROFSIZE;i++){
		output->time[i]  =this->time[i];
		output->flops[i] =this->flops[i];
		output->memory[i]=this->memory[i];
	}

	return (Object*)output;
}
/*}}}*/

/*Object virtual functions definitions:*/
void Profiler::DeepEcho(void){/*{{{*/
	this->Echo();
}/*}}}*/
void Profiler::Echo(void){/*{{{*/

	_printf_("Profiler:\n");
	for(int i=0;i<MAXPROFSIZE;i++){
		_printf_("    Tag "<<i<<":\n");
		_printf_("       flops:   "<<this->flops[i]<<"\n");
		_printf_("       memory:  "<<this->memory[i]<<"\n");
		_printf_("       time:    "<<this->time[i]<<"\n");
		_printf_("       running: "<<this->time[i]<<"\n");
	}

}
/*}}}*/
int  Profiler::Id(void){ /*{{{*/
	return -1; 
}
/*}}}*/
void Profiler::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	IssmPDouble* pointer = NULL;
	bool*       bpointer = NULL;

	int object_enum = ProfilerEnum;
	marshallhandle->call(object_enum);
	pointer = &this->time[0];
	marshallhandle->call(pointer,MAXPROFSIZE);
	pointer = &this->flops[0];
	marshallhandle->call(pointer,MAXPROFSIZE);
	pointer = &this->memory[0];
	marshallhandle->call(pointer,MAXPROFSIZE);
	bpointer = &this->running[0];
	marshallhandle->call(bpointer,MAXPROFSIZE);

} /*}}}*/
int  Profiler::ObjectEnum(void){/*{{{*/
	return ProfilerEnum;
}/*}}}*/

/*Profiler routines:*/
IssmPDouble  Profiler::TotalFlops(int tag){/*{{{*/

	/*Get tag*/
	_assert_(tag>=0); 
	_assert_(tag<MAXPROFSIZE); 
	if(this->running[tag]) _error_("Tag "<<tag<<" has not been stopped");

	return this->flops[tag];
}/*}}}*/
IssmPDouble  Profiler::TotalTime(int tag){/*{{{*/

	/*Get tag*/
	_assert_(tag>=0); 
	_assert_(tag<MAXPROFSIZE); 
	if(this->running[tag]) _error_("Tag "<<tag<<" has not been stopped");

	#ifdef _HAVE_MPI_
	return this->time[tag];
	#else
	return this->time[tag]/CLOCKS_PER_SEC;
	#endif
}
/*}}}*/
int Profiler::TotalTimeModHour(int tag){/*{{{*/

	IssmPDouble delta = this->TotalTime(tag);
	return int((reCast<int,IssmPDouble>(delta))/3600);

}
/*}}}*/
int Profiler::TotalTimeModMin(int tag){/*{{{*/

	IssmPDouble delta = this->TotalTime(tag);
	return int(int(reCast<int,IssmPDouble>(delta))%3600/60);
}
/*}}}*/
int Profiler::TotalTimeModSec(int tag){/*{{{*/

	IssmPDouble delta = this->TotalTime(tag);
	return int(reCast<int,IssmPDouble>(delta)%60);
}
/*}}}*/
IssmPDouble  Profiler::Memory(int tag){/*{{{*/

	/*Get initial flops*/
	_assert_(tag>=0); 
	_assert_(tag<MAXPROFSIZE); 
	if(this->running[tag]) _error_("Tag "<<tag<<" has not been stopped");

	return this->memory[tag];
}
/*}}}*/
void  Profiler::Start(int tag,bool dontmpisync){/*{{{*/

	/*Check tag*/
	_assert_(tag>=0); 
	_assert_(tag<MAXPROFSIZE); 
	if(this->running[tag]) _error_("Tag "<<tag<<" is already running");

	/*If mpisync requested, make sure all the cpus are at the same point in the execution: */
	if(!dontmpisync){
		ISSM_MPI_Barrier(IssmComm::GetComm()); 
	}

	/*Capture time: */
	#ifdef _HAVE_MPI_
	IssmPDouble t=ISSM_MPI_Wtime();
	#else
	IssmPDouble t=(IssmPDouble)clock();
	#endif

	/*Capture flops: */
	IssmPDouble f = 0.;
	IssmPDouble m = 0.;
	#ifdef _HAVE_PETSC_
		PetscGetFlops(&f);
		PetscMemoryGetCurrentUsage(&m);
	#else
		/*do nothing for now:*/
	#endif

	/*Plug into this->time: */
	_assert_(tag>=0); 
	_assert_(tag<MAXPROFSIZE); 
	this->time_start[tag]   = t;
	this->flops_start[tag]  = f;
	this->memory_start[tag] = m;

	/*turn on running*/
	this->running[tag] = true;
	this->used[tag]    = true;
}/*}}}*/
void  Profiler::Stop(int tag,bool dontmpisync){/*{{{*/

	/*Check tag*/
	_assert_(tag>=0); 
	_assert_(tag<MAXPROFSIZE); 
	if(!this->running[tag]) _error_("Tag "<<tag<<" is not running");

	/*If mpisync requested, make sure all the cpus are at the same point in the execution: */
	if(!dontmpisync){
		ISSM_MPI_Barrier(IssmComm::GetComm()); 
	}

	/*Capture time: */
	#ifdef _HAVE_MPI_
	IssmPDouble t=ISSM_MPI_Wtime();
	#else
	IssmPDouble t=(IssmPDouble)clock();
	#endif

	/*Capture flops: */
	IssmPDouble f = 0.;
	IssmPDouble m = 0.;
	#ifdef _HAVE_PETSC_
	PetscGetFlops(&f);
	PetscMemoryGetCurrentUsage(&m);
	#else
	/*do nothing for now:*/
	#endif

	/*Plug into this->time: */
	_assert_(tag>=0); 
	_assert_(tag<MAXPROFSIZE); 
	this->time[tag]   += t - this->time_start[tag];
	this->flops[tag]  += f - this->flops_start[tag];
	this->memory[tag] += m - this->memory_start[tag];

	/*turn off running*/
	this->running[tag] = false;
}/*}}}*/
bool  Profiler::Used(int tag){/*{{{*/

	/*Check tag*/
	_assert_(tag>=0); 
	_assert_(tag<MAXPROFSIZE); 
	return this->used[tag];
}/*}}}*/
