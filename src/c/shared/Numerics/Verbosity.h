/*!\file:Verbosity.h
 * \brief: Deal with verbosity levels
 */ 

#ifndef _VERBOSITY_H_
#define _VERBOSITY_H_

/*List of Verbosity levels (Add your own and Synchronize: must begin with "Verbose")*/
bool VerboseMProcessor(void);
bool VerboseModule(void);
bool VerboseSolution(void);
bool VerboseSolver(void);
bool VerboseConvergence(void);
bool VerboseControl(void);
bool VerboseQmu(void);
bool VerboseAutodiff(void);
bool VerboseSmb(void);

/*Setup Verbosity level*/
void SetVerbosityLevel(int level);
int  GetVerbosityLevel(void);

#endif
