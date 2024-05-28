/*\file diskio.h
 *\brief: I/O for ISSM from disk
 */

#ifndef _DISK_IO_H_
#define _DISK_IO_H_

#include <stdio.h>

FILE* pfopen(char* filename,const char* format,bool errorout=true);
FILE* pfopen0(char* filename,const char* format);
void  pfclose(FILE* fid,char* filename);
void WriteLockFile(char* filename);

#endif	/* _IO_H_ */
