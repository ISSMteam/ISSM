/*\file io.h
 *\brief: I/O for ISSM
 */

#ifndef _ISSM_IO_H_
#define _ISSM_IO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif 

#include "./Comm/IssmComm.h"
#include "./Disk/diskio.h"
#include "./Print/Print.h"
#include "./Marshalling/Marshalling.h"
#include "./Marshalling/IoCodeConversions.h"

#endif	/* _IO_H_ */
