/* \file Exceptions.cpp
 * \brief: implementation of the exceptions.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <cstring>
#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>
#include "./exceptions.h"

ErrorException::ErrorException(const string & what_arg){/*{{{*/

	int len;
	len           = strlen(what_arg.c_str())+1;
	what_str      = new char[len];
	memcpy(what_str,what_arg.c_str(),len);

	file_name     = NULL;
	function_name = NULL;
	file_line     = 0;

}/*}}}*/
ErrorException::ErrorException(int what_rank,const string& what_file, const string& what_function,int what_line, const string& what_arg){/*{{{*/

	/*Intermediaries*/
	int len;

	this->rank     = what_rank;
	this->file_line= what_line;

	len = strlen(what_arg.c_str())+1;
	this->what_str = new char[len];
	memcpy(this->what_str,what_arg.c_str(),len);

	len = strlen(what_file.c_str())+1;
	this->file_name = new char[len];
	memcpy(this->file_name,what_file.c_str(),len);

	len = strlen(what_function.c_str())+1;
	this->function_name = new char[len];
	memcpy(this->function_name,what_function.c_str(),len);

	/*Uncomment if messages do not print properly*/
	//this->Report();
}/*}}}*/
ErrorException::~ErrorException() throw(){/*{{{*/
	delete [] what_str;
	delete [] file_name;
	delete [] function_name;
}/*}}}*/
const char* ErrorException::what() const throw(){/*{{{*/
	return what_str;
}/*}}}*/
void ErrorException::Report() const{/*{{{*/

	cerr <<"\n[" << this->rank<< "] ??? Error using ==> " << this->file_name << ":" << this->file_line << 
	       "\n[" << this->rank<< "] " << this->function_name << " error message: " << what() << "\n" << endl;

	return;
}/*}}}*/
const char* ErrorException::WrapperReport() const{/*{{{*/

	/*Output*/
	std::ostringstream buffer;

	buffer << "\nError in ==> " << this->file_name << ":" << file_line << "\n";
	buffer << this->function_name << " error message: " << this->what_str;

	/*Convert std::ostringstream to std::string and then create char* */
	std::string buffer2 = buffer.str();
	int   message_len = strlen(buffer2.c_str())+1;
	char* message = new char[message_len];
	snprintf(message, message_len,"%s",buffer2.c_str());

	return message;
}/*}}}*/
