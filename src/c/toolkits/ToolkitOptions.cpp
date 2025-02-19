/*! \file ToolkitOptions.cpp
 * \brief  file containing the methods for ToolkitOptions.h
 */

#include <string.h>
#include "./ToolkitOptions.h"
#include "../shared/Numerics/types.h"
#include "../shared/Exceptions/exceptions.h"
#include "../shared/MemOps/MemOps.h"

#ifdef _DO_NOT_LOAD_GLOBALS_
char* ToolkitOptions::toolkittype;
char* ToolkitOptions::toolkitoptions;
#endif

void  ToolkitOptions::Init(const char* toolkit_in,const char* options){ /*{{{*/

	/*First, avoid a leak: */
	xDelete<char>(toolkitoptions);
	xDelete<char>(toolkittype);

	/*copy options into toolkitoptions:*/
	_assert_(toolkit_in);
	_assert_(options);
	int toolkittype_len    = strlen(toolkit_in)+1;
	int toolkitoptions_len =  strlen(options)+1;
	toolkittype    = xNew<char>(toolkittype_len); 
	toolkitoptions = xNew<char>(toolkitoptions_len);
	snprintf(toolkittype, toolkittype_len,"%s",toolkit_in);
	snprintf(toolkitoptions, toolkitoptions_len,"%s",options);
}/*}}}*/
void  ToolkitOptions::Init(){ /*{{{*/
	toolkittype    = NULL;
	toolkitoptions = NULL;
}/*}}}*/
void  ToolkitOptions::Delete(){ /*{{{*/

	xDelete<char>(toolkitoptions);
	xDelete<char>(toolkittype);

}/*}}}*/
char* ToolkitOptions::GetToolkitType(){  /*{{{*/

	if(toolkittype==NULL) _error_("toolkittype not set (may be a mex?)");
	int   toolkittype_out_len = strlen(toolkittype)+1;
	char* toolkittype_out = xNew<char>(toolkittype_out_len); 
	snprintf(toolkittype_out, toolkittype_out_len, "%s", toolkittype);
	return toolkittype_out;
}/*}}}*/
char* ToolkitOptions::GetToolkitOptionValue(const char* option){  /*{{{*/

	return TokenValue(toolkitoptions,option);

}/*}}}*/
char* TokenValue(char* tokenlist,const char* target){ /*{{{*/

	/*output:*/
	char* value=NULL;

	/*intermediary: */
	char *token         = NULL;
	char *tokenlistcopy = NULL;

	/*First, because tokenizing destroys a string, copy what we have: */
	if(tokenlist==NULL) _error_("tokenlist not set (may be a mex?)");
	int tokenlist_len = strlen(tokenlist)+1;
	tokenlistcopy= xNew<char>(tokenlist_len);
	snprintf(tokenlistcopy, tokenlist_len,"%s",tokenlist);

	/*Now go through list of tokens, and look for  target, return value: */
	token=strtok(tokenlistcopy, " ");
	while(token != NULL) {

		/*Is this token starting with "-", if so, compare to our target: */
		if (strncmp(token,"-",1)==0){
			if (strcmp(token+1,target)==0){
				/*Ok, we found our target. Get next token: */
				token = strtok(NULL, " ");
				/*This token could actually be another option start with "-", just be sure: */
				if (strncmp(token,"-",1)==0){
					/*ok, we hit another option, which means our target value is "":*/
					int value_len = strlen("")+1;
					value = xNew<char>(value_len);
					snprintf(value, value_len,"%s","");
					continue;
				}
				else{
					/*this token is the value we are looking for, copy: */
					int value_len = strlen(token)+1;
					value= xNew<char>(value_len);
					snprintf(value, value_len,"%s",token);
				}
			}
			else{
				/*we found the wrong target. Go to the next option: */
				token = strtok(NULL, " ");
				if (strncmp(token,"-",1)==0){
					/*this is indeed an option, continue: */
					continue;
				}
				else{
					/*this is the value of the option, discard it: */
				}
			}
		}
		else _error_("token list should start with an option, not a value");

		/*Get new token and continue*/
		token = strtok(NULL, " ");
	}

	/*Clean up and return*/
	xDelete<char>(tokenlistcopy);
	return value;
}
/*}}}*/
