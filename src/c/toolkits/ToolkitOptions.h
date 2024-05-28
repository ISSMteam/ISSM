/* \file ToolkitOptions.h
 * \brief  create a class with a static string of options, and static methods to access it
 * This is a way of protecting access to the toolkit options, and to make it accessible everywhere
 * in the code.
 */

#ifndef _TOOLKIT_OPTIONS_H
#define _TOOLKIT_OPTIONS_H

class ToolkitOptions {

	private:
		static char* toolkittype;
		static char* toolkitoptions;

	public:
		static void  Init(const char* type_in,const char* options);
		static void  Init(void);
		static void  Delete(void);
		static char* GetToolkitType(void);
		static char* GetToolkitOptionValue(const char* option);
};

char* TokenValue(char* tokenlist,const char* target);

#endif  /* _TOOLKIT_OPTIONS_H */
