/*!\file KMLFileUtils.cpp
 * \brief: utilities for kml file reading.
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KMLFileReadUtils.h"
#include "./KML_Object.h"
#include "../shared/shared.h"
/*}}}*/

char* KMLFileToken(FILE* fid,/*{{{*/
				   int* pncom=NULL,char*** ppcom=NULL){

/*  get the next token (tag or field) in the file  */

	bool    inew=1,itag=0,ifield=0;
	int     c;
	int     ibuf=0,buflen=1024,bufblk=1024;
	char    *buffer=NULL,*bufferc=NULL,**pcom2=NULL;

	buffer=xNew<char>(buflen);
	buffer[0]='\0';

/*  read kml file character-by-character  */

//  note that fgets includes newline
//	fgets(buffer,buflen,fid);

	while ((c=getc(fid)) != EOF) {
		/*  ignore leading blanks  */
		if (inew && isspace(c))
			continue;

		/*  distinguish between tag or field  */
		if (!itag && !ifield) {

			/*  distinguish between tag or comment  */
			if (c == '<') {
				ungetc(c,fid);
				if (!(bufferc=KMLFileTokenComment(fid))) {
					c=getc(fid);
					itag=1;
				}
				else {
					if (pncom && ppcom) {
						(*pncom)++;
						pcom2=xNew<char*>(*pncom);
						memcpy(pcom2,*ppcom,(*pncom-1)*sizeof(char*));
						xDelete<char*>(*ppcom);
						*ppcom=pcom2;
						pcom2=NULL;
//						*ppcom=(char **) xrealloc(*ppcom,*pncom*sizeof(char*));
						(*ppcom)[*pncom-1]=bufferc;
					}
					else
						xDelete<char>(bufferc);
					inew=1;
					continue;
				}
			}
			else
				ifield=1;
			inew=0;
			KMLFileTokenBuffer(&buffer,&ibuf,&buflen,
							   c,
							   bufblk);
		}

		/*  accumulate tag, not including newlines  */
		else if (itag) {
			if (c != '\n') {
				inew=0;
				KMLFileTokenBuffer(&buffer,&ibuf,&buflen,
								   c,
								   bufblk);
				if (c == '>')
					break;
			}
			else
				inew=1;
		}

		/*  accumulate field, including newlines  */
		else if (ifield) {
			/*  distinguish between another tag or comment  */
			if (c == '<') {
				ungetc(c,fid);
				if (!(bufferc=KMLFileTokenComment(fid)))
					break;
				else
					if (pncom && ppcom) {
						(*pncom)++;
						pcom2=xNew<char*>(*pncom);
						memcpy(pcom2,*ppcom,(*pncom-1)*sizeof(char*));
						xDelete<char*>(*ppcom);
						*ppcom=pcom2;
						pcom2=NULL;
//						*ppcom=(char **) xrealloc(*ppcom,*pncom*sizeof(char*));
						(*ppcom)[*pncom-1]=bufferc;
					}
					else
						xDelete<char>(bufferc);
			}
			else {
				inew=0;
				KMLFileTokenBuffer(&buffer,&ibuf,&buflen,
								   c,
								   bufblk);
				if (c == '\n')
					inew=1;
			}
		}

	}

/*  remove trailing blanks or newline  */

	while (ibuf > 0)
		if (isspace(buffer[ibuf-1]))
			ibuf--;
		else {
			buffer[ibuf]='\0';
			break;
		}

//	if      (itag)
//		_printf0_("tag buffer (length=" << ibuf << "):\n");
//	else if (ifield)
//		_printf0_("field buffer (length=" << ibuf << "):\n");
//	_printf0_(buffer << "\n");

	if (!ibuf)
		xDelete<char>(buffer);

	return(buffer);
}
/*}}}*/
char* KMLFileTokenComment(FILE* fid){/*{{{*/

/*  check for comment in the file and read it  */

	bool    inew=1;
	int     i;
	int     c;
	int     ibuf=0,buflen=1024,bufblk=1024;
	char*   buffer=NULL;

	buffer=xNew<char>(buflen);
	buffer[0]='\0';

/*  read kml file character-by-character  */

	while ((c=getc(fid)) != EOF) {
		/*  ignore leading blanks  */
		if (inew && isspace(c))
			continue;

		inew=0;
		KMLFileTokenBuffer(&buffer,&ibuf,&buflen,
						   c,
						   bufblk);

		/*  check for comment  */
		if (ibuf <= 4) {
			if ((ibuf == 1 && buffer[0] != '<') ||
				(ibuf == 2 && buffer[1] != '!') ||
				(ibuf == 3 && buffer[2] != '-') ||
				(ibuf == 4 && buffer[3] != '-')) {
				for (i=ibuf-1; i>=0; i--)
					ungetc(buffer[i],fid);
				xDelete<char>(buffer);
				return(buffer);
			}
		}

		/*  accumulate comment, including newlines  */
		else
			if (buffer[ibuf-3]=='-' && buffer[ibuf-2]=='-' && buffer[ibuf-1]=='>')
				break;
	}

/*  remove trailing blanks or newline  */

	while (ibuf > 0)
		if (isspace(buffer[ibuf-1]))
			ibuf--;
		else {
			buffer[ibuf]='\0';
			break;
		}

//	_printf0_("comment buffer (length=" << ibuf << "):\n");
//	_printf0_(buffer << "\n");

	if (!ibuf)
		xDelete<char>(buffer);

	return(buffer);
}
/*}}}*/
void KMLFileTokenBuffer(char** pbuffer,int* pibuf,int* pbuflen,/*{{{*/
						int c,
						int bufblk){

/*  add the specified character to the token buffer  */

	char*   buffer2=NULL;

/*  check buffer length and realloc if necessary  */

	if (*pibuf+2 > *pbuflen) {
		*pbuflen+=bufblk;
		buffer2=xNew<char>(*pbuflen);
		memcpy(buffer2,*pbuffer,(*pbuflen-bufblk)*sizeof(char));
		xDelete<char>(*pbuffer);
		*pbuffer=buffer2;
		buffer2=NULL;
//		*pbuffer=(char *) xrealloc(*pbuffer,*pbuflen*sizeof(char));
	}

/*  add character and terminator  */

	(*pbuffer)[(*pibuf)++]=c;
	(*pbuffer)[ *pibuf   ]='\0';

	return;
}
/*}}}*/
char* KMLFileTagName(char* pname,/*{{{*/
					 char* ktag){

	return(KMLFileTagName(pname,NULL,0,
						  ktag));
}
/*}}}*/
char* KMLFileTagName(char* pname,int *m,int maxlen,/*{{{*/
					 char* ktag){

/*  for the given tag buffer, read and store the name  */

	char*   ktagi;
	char*   ktokn;

	if (strncmp(&ktag[0],"<"        ,1) || strncmp(&ktag[strlen(ktag)-1],">",1))
		_error_("KMLFileTagName -- Missing tag delimiters in " << ktag << ".\n");

/*  strtok modifies ktag, so work on copy  */

	ktagi=xNew<char>(strlen(ktag)+1);
	memcpy(ktagi,ktag,(strlen(ktag)+1)*sizeof(char));

/*  skip opening delimeter and find subsequent blank or closing delimiter  */

	ktokn=strtok(ktagi,"< >");
//	_printf0_("KMLFileTagName -- initial token=\"" << ktokn << "\".\n");

	if (!pname) {
		if (maxlen)
			pname=xNew<char>(maxlen       +1);
		else
			pname=xNew<char>(strlen(ktokn)+1);
	}

	if (maxlen && (maxlen < strlen(ktokn))) {
		_printf0_("KMLFileTagName -- string field too short for " << ktag << ".\n");
		_printf0_("KMLFileTagName -- \"" << ktokn << "\" truncated to " << maxlen << " characters.\n");
		strncpy(pname,ktokn,maxlen);
	}
	else
		memcpy(pname,ktokn,(strlen(ktokn)+1)*sizeof(char));

	xDelete<char>(ktagi);

	if (m)
		*m=strlen(pname);

	return(pname);
}
/*}}}*/
int KMLFileTagAttrib(KML_Object* kobj,/*{{{*/
					 char* ktag){

/*  for the given tag buffer, read and store the attributes  */

	char*   ktagi;
	char*   ktokn;
	char*   ktokv;
	char    quote[]={'\"','\0'};
	int     isolo=0;

/*  strtok modifies ktag, so work on copy  */

	ktagi=xNew<char>(strlen(ktag)+1);
	memcpy(ktagi,ktag,(strlen(ktag)+1)*sizeof(char));

/*  loop through tag to find all attributes  */

	/*  return first non blank and move past subsequent blank  */
	ktokn=strtok(ktagi," ");
//	_printf0_("KMLFileTagAttrib -- initial token=\"" << ktokn << "\".\n");

	/*  return next non " =?/>" and move past subsequent " =?/>"  */
	while((ktokn=strtok(NULL," =?/>"))){

		/*  return next non quote and move past subsequent quote  */
		ktokv=strtok(NULL,quote);
//		_printf0_("KMLFileTagAttrib -- attribute " << ktokn << "=\"" << ktokv << "\".\n");

/*  add the attribute to the dataset  */

		if (kobj)
			kobj->AddAttrib(ktokn,ktokv);
	}

	xDelete<char>(ktagi);

/*  check for xml declaration, dtd declaration, or solo tag  */

	if ((!strncmp(&ktag[0],"<?"       ,2) && !strncmp(&ktag[strlen(ktag)-2],"?>",2)) ||
		(!strncmp(&ktag[0],"<!DOCTYPE",9) && !strncmp(&ktag[strlen(ktag)-1], ">",1)) ||
		(!strncmp(&ktag[0],"<"        ,1) && !strncmp(&ktag[strlen(ktag)-2],"/>",2)))
		isolo=1;
//	_printf0_("KMLFileTagAttrib -- isolo=" << isolo << ".\n");

	return(isolo);
}
/*}}}*/
int KMLFileTokenParse(int* pival,/*{{{*/
					  char* ktag,
					  FILE* fid){

	char*   kstr;

/*  get next token and convert to appropriate format  */

	if (!(kstr=KMLFileToken(fid,
							NULL,NULL)) ||
		(kstr[0] == '<'))
		_error_("KMLFileTokenParse -- Missing integer field for " << ktag << ".\n");

	sscanf(kstr,"%d",pival);
	xDelete<char>(kstr);

/*  get additional token and compare to closing tag  */

	if(ktag){
		if (!(kstr=KMLFileToken(fid,
								NULL,NULL)) ||
			(kstr[0] != '<') ||
			(kstr[1] != '/') ||
			(strncmp(&(kstr[2]),&(ktag[1]),strlen(ktag)-1)))
		  {_error_("KMLFileTokenParse -- Missing closing tag for " << ktag << ".\n");}
		else{
			xDelete<char>(kstr);
		}
	}

//	_printf0_("KMLFileTokenParse -- " << ktag << "=" << *pival << ".\n");

	return(0);
}
/*}}}*/
int KMLFileTokenParse(bool* pbval, char* ktag, FILE* fid){/*{{{*/

	int     ival;
	char*   kstr;

/*  get next token and convert to appropriate format  */

	if (!(kstr=KMLFileToken(fid,
							NULL,NULL)) ||
		(kstr[0] == '<'))
	  {_error_("KMLFileTokenParse -- Missing bool field for " << ktag << ".\n");}

	sscanf(kstr,"%d",&ival);
	*pbval=(bool)ival;
	xDelete<char>(kstr);

/*  get additional token and compare to closing tag  */

	if(ktag){
		if (!(kstr=KMLFileToken(fid,
								NULL,NULL)) ||
			(kstr[0] != '<') ||
			(kstr[1] != '/') ||
			(strncmp(&(kstr[2]),&(ktag[1]),strlen(ktag)-1)))
		  {_error_("KMLFileTokenParse -- Missing closing tag for " << ktag << ".\n");}
		else
			xDelete<char>(kstr);
	}

//	_printf0_("KMLFileTokenParse -- " << ktag << "=" << (*pbval ? "true" : "false") << ".\n");

	return(0);
}
/*}}}*/
char* KMLFileTokenParse(char* pstr,/*{{{*/
						char* ktag,
						FILE* fid){

	return(KMLFileTokenParse(pstr,NULL,0,
							 ktag,
							 fid));
}
/*}}}*/
char* KMLFileTokenParse(char* pstr,int *m,int maxlen,/*{{{*/
						char* ktag,
						FILE* fid){

	char*   kstr;

/*  get next token and allocate if necessary  */

	if (!(kstr=KMLFileToken(fid,
							NULL,NULL)) ||
		(kstr[0] == '<'))
		_error_("KMLFileTokenParse -- Missing string field for " << ktag << ".\n");

	if (!pstr) {
		if (maxlen)
			pstr=xNew<char>(maxlen      +1);
		else
			pstr=xNew<char>(strlen(kstr)+1);
	}

	if (maxlen && (maxlen < strlen(kstr))) {
		_printf0_("KMLFileTokenParse -- string field too short for " << ktag << ".\n");
		_printf0_("KMLFileTokenParse -- \"" << kstr << "\" truncated to " << maxlen << " characters.\n");
		strncpy(pstr,kstr,maxlen);
	}
	else
		memcpy(pstr,kstr,(strlen(kstr)+1)*sizeof(char));

	xDelete<char>(kstr);

	if (m)
		*m=strlen(pstr);

/*  get additional token and compare to closing tag  */

	if(ktag){
		if (!(kstr=KMLFileToken(fid,
								NULL,NULL)) ||
			(kstr[0] != '<') ||
			(kstr[1] != '/') ||
			(strncmp(&(kstr[2]),&(ktag[1]),strlen(ktag)-1)))
		  {_error_("KMLFileTokenParse -- Missing closing tag for " << ktag << ".\n");}
		else
			xDelete<char>(kstr);
	}

//	_printf0_("KMLFileTokenParse -- " << ktag << "=\"" << pstr << "\".\n");

	return(pstr);
}
/*}}}*/
int KMLFileTokenParse(float* pfval,/*{{{*/
					  char* ktag,
					  FILE* fid){

	char*   kstr;

/*  get next token and convert to appropriate format  */

	if (!(kstr=KMLFileToken(fid,
							NULL,NULL)) ||
		(kstr[0] == '<'))
	  {_error_("KMLFileTokenParse -- Missing integer field for " << ktag << ".\n");}

	sscanf(kstr,"%g",pfval);
	xDelete<char>(kstr);

/*  get additional token and compare to closing tag  */

	if(ktag){
		if (!(kstr=KMLFileToken(fid,
								NULL,NULL)) ||
			(kstr[0] != '<') ||
			(kstr[1] != '/') ||
			(strncmp(&(kstr[2]),&(ktag[1]),strlen(ktag)-1)))
		  {_error_("KMLFileTokenParse -- Missing closing tag for " << ktag << ".\n");}
		else
			xDelete<char>(kstr);
	}

//	_printf0_("KMLFileTokenParse -- " << ktag << "=" << *pfval << ".\n");

	return(0);
}
/*}}}*/
int KMLFileTokenParse(double* pdval,/*{{{*/
					  char* ktag,
					  FILE* fid){

	char*   kstr;

/*  get next token and convert to appropriate format  */

	if (!(kstr=KMLFileToken(fid,
							NULL,NULL)) ||
		(kstr[0] == '<'))
		_error_("KMLFileTokenParse -- Missing integer field for " << ktag << ".\n");

	sscanf(kstr,"%lg",pdval);
	xDelete<char>(kstr);

/*  get additional token and compare to closing tag  */

	if(ktag){
		if (!(kstr=KMLFileToken(fid,
								NULL,NULL)) ||
			(kstr[0] != '<') ||
			(kstr[1] != '/') ||
			(strncmp(&(kstr[2]),&(ktag[1]),strlen(ktag)-1)))
		  {_error_("KMLFileTokenParse -- Missing closing tag for " << ktag << ".\n");}
		else
			xDelete<char>(kstr);
	}

//	_printf0_("KMLFileTokenParse -- " << ktag << "=" << *pdval << ".\n");

	return(0);
}
/*}}}*/
int KMLFileTokenParse(double **pdval,int* m,int maxlen,/*{{{*/
					  char* ktag,
					  FILE* fid){

	int     i=-1;
	char*   kstr;
	char*   ktok;
	double* dval2=NULL;
	char    delim[]={' ',',','\f','\n','\r','\t','\v','\0'};

/*  get next token and allocate if necessary  */

	if (!(kstr=KMLFileToken(fid,
							NULL,NULL)) ||
		(kstr[0] == '<'))
		_error_("KMLFileTokenParse -- Missing double [m] field for " << ktag << ".\n");

	if(!*pdval){
		if (maxlen)
			*pdval=xNew<IssmPDouble>(maxlen            );
		else
			*pdval=xNew<IssmPDouble>((strlen(kstr)+1)/2);
	}

/*  loop through string to get all values  */

	ktok=strtok(kstr,delim);
	while (ktok) {
		i++;
		if (maxlen && (maxlen < i+1))
			_error_("KMLFileTokenParse -- Double [m] field too short for " << ktag << ".\n");
		sscanf(ktok,"%lg",&((*pdval)[i]));
		ktok=strtok(NULL,delim);
	}
	xDelete<char>(kstr);

	if (!maxlen)
		dval2=xNew<double>(i+1);
		memcpy(dval2,*pdval,(i+1)*sizeof(double));
		xDelete<double>(*pdval);
		*pdval=dval2;
		dval2=NULL;
//		*pdval=(double *) xrealloc(*pdval,(i+1)*sizeof(double));

	if (m)
		*m=i+1;

/*  get additional token and compare to closing tag  */

	if(ktag){
		if (!(kstr=KMLFileToken(fid,
								NULL,NULL)) ||
			(kstr[0] != '<') ||
			(kstr[1] != '/') ||
			(strncmp(&(kstr[2]),&(ktag[1]),strlen(ktag)-1)))
		  {_error_("KMLFileTokenParse -- Missing closing tag for " << ktag << ".\n");}
		else
			xDelete<char>(kstr);
	}

//	_printf0_("KMLFileTokenParse -- " << ktag << "=...\n");
//	for (j=0; j<=i; j++)
//		_printf0_("   [" << j << "]: " << (*pdval)[j] << "g\n");

	return(0);
}
/*}}}*/
int KMLFileTokenParse(double **pdval,int* m,int n,int maxlen,/*{{{*/
					  char* ktag,
					  FILE* fid){

	int     i=-1,j=-1;
	char*   kstr;
	char*   ktok;
	double* dval2=NULL;
	char    delim[]={' ',',','\f','\n','\r','\t','\v','\0'};

/*  get next token and allocate if necessary  */

	if (!(kstr=KMLFileToken(fid,
							NULL,NULL)) ||
		(kstr[0] == '<'))
		_error_("KMLFileTokenParse -- Missing double [m x n] field for " << ktag << ".\n");

	if(!*pdval){
		if (maxlen)
			*pdval=xNew<IssmPDouble>(maxlen*n          );
		else
			*pdval=xNew<IssmPDouble>((strlen(kstr)+1)/2);
	}

/*  loop through string to get all values  */

	ktok=strtok(kstr,delim);
	while (ktok) {
		i++;
		if (maxlen && (maxlen*n < i+1))
			_error_("KMLFileTokenParse -- Double [m x n] field too short for " << ktag << ".\n");
		j=(j+1) % n;
		sscanf(ktok,"%lg",&((*pdval)[i]));
		ktok=strtok(NULL,delim);
	}
	xDelete<char>(kstr);

	if (!maxlen)
		dval2=xNew<double>((i+1)*n);
		memcpy(dval2,*pdval,((i+1)*n)*sizeof(double));
		xDelete<double>(*pdval);
		*pdval=dval2;
		dval2=NULL;
//		*pdval=(double *) xrealloc(*pdval,((i+1)*n)*sizeof(double));

	if (m)
		*m=((i+1)+(n-1))/n;

	if ((i+1) % n)
		_printf0_("KMLFileTokenParse -- Double [m x n] field for " << ktag << " does not have multiple of n values.\n");

/*  get additional token and compare to closing tag  */

	if(ktag){
		if (!(kstr=KMLFileToken(fid,
								NULL,NULL)) ||
			(kstr[0] != '<') ||
			(kstr[1] != '/') ||
			(strncmp(&(kstr[2]),&(ktag[1]),strlen(ktag)-1)))
		  {_error_("KMLFileTokenParse -- Missing closing tag for " << ktag << ".\n");}
		else
			xDelete<char>(kstr);
	}

//	_printf0_("KMLFileTokenParse -- " << ktag << "=...\n");
//	for (j=0; j<=i; j++)
//		_printf0_("   [" << j << "]: " << (*pdval)[j] << "g\n");

	return(0);
}
/*}}}*/
int KMLFileTagSkip(char* ktag, FILE* fid){/*{{{*/

	char*   kstr;

/*  note that tags of the same type can be nested inside each other, so for each
	opening tag, must find corresponding closing tag  */

	_printf0_("KMLFileTagSkip -- input tag " << ktag << ".\n");

/*  if next token is a closing tag, compare to input  */

	while((kstr=KMLFileToken(fid,NULL,NULL))){
		if((kstr[0] == '<') && (kstr[1] == '/') && (!strncmp(&(kstr[2]),&(ktag[1]),(strcspn(ktag," >")-1)/sizeof(char)))){
			_printf0_("KMLFileTagSkip -- closing tag " << kstr << ".\n");
			xDelete<char>(kstr);
			return(0);
		}

/*  if next token is an opening tag, call recursively  */

		else if ((kstr[0] == '<') &&
				 (kstr[1] != '/')) {
			_printf0_("KMLFileTagSkip -- opening tag " << kstr << ".\n");
			KMLFileTagSkip(kstr,
						   fid);
		}

/*  if next token is a closing tag, error out  */

		else if ((kstr[0] == '<') &&
				 (kstr[1] == '/')) {
			_error_("KMLFileTagSkip -- Unexpected closing tag " << kstr << ".\n");
		}

		xDelete<char>(kstr);
	}

	_error_("KMLFileTokenParse -- Corresponding closing tag for " << ktag << " not found.\n");

	return(0);
}
/*}}}*/
