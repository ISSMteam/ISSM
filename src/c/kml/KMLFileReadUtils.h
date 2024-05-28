/*!\file:  KMLFileReadUtils.h
 * \brief: header file for kml file reading utilities.
 */ 

#ifndef _KMLFILEREADUTILS_H
#define _KMLFILEREADUTILS_H

/*Headers:{{{*/
#include "../shared/shared.h"
#include "../datastructures/datastructures.h"
class KML_Object;
/*}}}*/

/* local prototypes: */
char* KMLFileToken(FILE* fid, int* pncom,char*** ppcom);
char* KMLFileTokenComment(FILE* fid);
void KMLFileTokenBuffer(char** pbuffer,int* pibuf,int* pbuflen, int c, int bufblk);
char* KMLFileTagName(char* pname, char* ktag);
char* KMLFileTagName(char* pname,int *m,int maxlen, char* ktag);
int KMLFileTagAttrib(KML_Object* kobj, char* ktag);
int KMLFileTokenParse(int* pival, char* ktag, FILE* fid);
int KMLFileTokenParse(bool* pbval, char* ktag, FILE* fid);
char* KMLFileTokenParse(char* pstr, char* ktag, FILE* fid);
char* KMLFileTokenParse(char* pstr,int *m,int maxlen, char* ktag, FILE* fid);
int KMLFileTokenParse(float* pfval, char* ktag, FILE* fid);
int KMLFileTokenParse(double* pdval, char* ktag, FILE* fid);
int KMLFileTokenParse(double **pdval,int* m,int maxlen, char* ktag, FILE* fid);
int KMLFileTokenParse(double **pdval,int* m,int n,int maxlen, char* ktag, FILE* fid);
int KMLFileTagSkip(char* ktag, FILE* fid);

#endif  /* _KMLFILEREADUTILS_H */
