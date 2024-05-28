/*********************************************************************
Copyright 2008, 2010 Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
retains certain rights in this software.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of Sandia Corporation nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
***********************************************************************/

/* nidr.c */

#ifndef NIDR_H	/* for $DAKOTA/src/nidr.c */
#include "nidr.h"
#endif
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "avltree.h"

#ifndef NIDR_SQUAWKMAX
#define NIDR_SQUAWKMAX 10
#endif

#ifndef NO_NIDR_DYNLIB /*{*/
typedef KeyWord *(*KW_ADD)(void);
static KeyWord *kwfind(const char *, KeyWord *, int, int *);
static KeyWord *toomany(const char *, KeyWord *, int);
#ifdef _WIN32 /*{{*/
#include <windows.h>
#define dlopen(x,y) LoadLibrary(x)
#define find_dlsym(a,b,c) (a = (KW_ADD)GetProcAddress((HINSTANCE)(b),c))
#define dlclose(x) FreeLibrary((HMODULE)x)
#define NO_DLERROR
#else /*}{*/
#include <dlfcn.h>
#include <unistd.h>
#define find_dlsym(a,b,c) (a = (KW_ADD)dlsym(b,c))
#undef NO_DLERROR
#endif /*}}*/
#endif /*}*/

 extern KeyWord Dakota_Keyword_Top;
 extern int nidrLineNumber;
 static KeyWord* Keyword_Top = &Dakota_Keyword_Top;
 static void *KW_g;
 void (*nidr_comment)(const char*);
 static void nidr_keyword_finish(void);
 static Comment *OutsideComment;
 static void kw_finish2(void), kw_finish3(void);
 static void kw_setup1(KeyWord *);
 static FILE *dumpfile;
 static KeyWord **ToClear, **ToClear0, **ToClearEnd;
 static int dumplev, nsquawk, nparse_errors, primary, strict;

 int NIDR_disallow_missing_start = 1;

 enum {n_KWStack0 = 64};

 static KWinfo KWStack0[n_KWStack0];

 static Uint n_KWStack = n_KWStack0;

 static KeyWord *curid, *curkw;
 static KWinfo	*KWStack = KWStack0,
		*KWStackBot = KWStack0,
		*KWStackEnd = KWStack0 + n_KWStack0;

 static Values KWval, KWvalmax;
 static Real *KWvalbuf;
 static Uint nKWvalbuf;

 typedef struct Sbuf Sbuf;
 enum { n_KWsbuf = 8192 };
 struct Sbuf {
	char buf[n_KWsbuf];
	Sbuf *next;
	};

 typedef struct KWseen KWseen;

 struct
KWseen {
	const char *name;
	KeyWord *kw;
	KWseen *mnext, *mprev;	/* for identifiers so far unrequited when kw == 0 */
				/* kw != 0 ==> mprev = first child, and mnext = next sibling */
	KWseen *parent;
	KWseen **lcn;		/* &mprev field of last child; lcn == 0 when this */
				/* keyword was entered into the AVL tree because */
				/* its parent was seen. */
	Comment *comment;
	char **svals;
	Real *rvals;
	size_t nvals;
	};

 static KWseen *KW_cur;
 NIDR_KWlib *NIDR_Libs;

 void
nidr_lib_cleanup(void)
{
	KeyWord *kw;
	NIDR_KWlib *L, *L1;

	L1 = NIDR_Libs;
	NIDR_Libs = 0;
	while((L = L1)) {
		if (L->oldtop)
			Keyword_Top = L->oldtop;
		if ((kw = L->kw0)) {
			kw->f.vs = 0;
			kw->kind &= ~KWKind_Loaded;
			}
#ifndef NO_NIDR_DYNLIB /*{{*/
		dlclose(L->h);
#else  /*}{*/
		/* botch("dlclose is NOT SUPPORTED for current configuration"); */
		fprintf(stderr, "\ndlclose is NOT SUPPORTED for current configuration");
#endif /*}}*/
		L1 = L->next;
		free(L);
		}
	}

 static Sbuf KWsbuf0, *KWsbuf = &KWsbuf0;
 static char *KWsbuf1 = KWsbuf0.buf, *KWsbufe = KWsbuf0.buf + n_KWsbuf;
 static KWseen *curkws;
 static const char *valkind[3] = {"integer","numeric","string"};

 int
nidr_parse_error(void)
{
	int n;
	if ((n = nsquawk - NIDR_SQUAWKMAX) > 0)
		fprintf(stderr, "\n%d error message%s suppressed.\n",
			n, "s" + (n == 1));
	return nsquawk + nparse_errors;
	}

 void
nidr_signal_parse_error(void)
{ ++nparse_errors; }

 void
nidr_tolower(char *s)
{
	for(; *s; ++s)
		*s = tolower(*s);
	}

 static void
botch(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	fprintf(stderr, "\nBotch:  ");
	vfprintf(stderr, fmt, ap);
	fputs(".\n", stderr);
	va_end(ap);
	exit(1);
	}

 static void
squawk(const char *fmt, ...)
{
	va_list ap;
	if (++nsquawk <= NIDR_SQUAWKMAX) {
		fprintf(stderr, "Input line %d: ", nidrLineNumber);
		va_start(ap, fmt);
		vfprintf(stderr, fmt, ap);
		fputs(".\n", stderr);
		va_end(ap);
		}
	}

#ifdef NIDR_MALLOC_DEBUG
 typedef struct MallocDebug MallocDebug;
 struct MallocDebug
{
	MallocDebug *next, *prev;
	char *where;
	int nalloc;
	};

 static MallocDebug MDtop = {&MDtop, &MDtop, 0, 0};
 int MallocDebugCount, MallocDebugCount1;

 static void*
Alloc(const char *where, size_t len)
{
	MallocDebug *md = malloc(len + sizeof(MallocDebug));
	if (!md) {
		fprintf(stderr, "malloc(%lu) failure in %s\n", (unsigned long)len, where);
		exit(1);
		}
	(md->next = MDtop.next)->prev = md;
	(md->prev = &MDtop)->next = md;
	md->where = where;
	if ((md->nalloc = ++MallocDebugCount) == MallocDebugCount1)
		printf("Hit %d\n", md->nalloc);
	return (void*)(md + 1);
	}

 static void
MallocDebugFree(void *v)
{
	MallocDebug *md = (MallocDebug *)v - 1;
	md->next->prev = md->prev;
	md->prev->next = md->next;
	free(md);
	}
#define free MallocDebugFree

#else //!NIDR_MALLOC_DEBUG
 static void*
Alloc(const char *where, size_t len)
{
	void *rv = malloc(len);
	if (!rv) {
		fprintf(stderr, "malloc(%lu) failure in %s\n", (unsigned long)len, where);
		exit(1);
		}
	return rv;
	}
#endif //NIDR_MALLOC_DEBUG

 const char *
nidr_basename(const char *p)
{
	const char *b;

#ifdef _WIN32
	if (p[0] && p[1] == ':')
		p += 2;
	else if (p[0] == '\\' && p[1] == '\\')
		for(p += 2; *p; )
			switch(*p++) {
			 case '/':
			 case '\\':
				goto break2;
			}
 break2:
#endif
	b = p;
	while(*p)
		switch(*p++) {
		 case '/':
#ifdef _WIN32
		 case '\\':
#endif
			b = p;
		}
	return b;
	}

 const char *nidr_exedir;
 int nidr_exedir_len = -1; /* allow resetting to -1 for debugging */

#ifndef _WIN32
 static int
Is_executable(uid_t myuid, gid_t mygid, struct stat *sb)
{
	if (sb->st_uid == myuid) {
		if (sb->st_mode & S_IXUSR)
			return 1;
		}
	else if (sb->st_gid == mygid) {
		if (sb->st_mode & S_IXGRP)
			return 1;
		}
	else if (sb->st_mode & S_IXOTH)
		return 1;
	return 0;
	}
#endif

 int
nidr_save_exedir(const char *argv0, int pathadj)
{
	/* pathadj & 1 ==> add exedir to $PATH */
	/* pathadj & 2 ==> add . to $PATH */
	/* (in both cases if not already there) */

    /* These conditionals don't seem to work (perhaps expected) for Cygwin 
	   binaries run from Windows command prompt as the compile-time is 
	   unix-style, but runtime the path is windows-like.  For now comment
	   out warning when on Cygwin build */
#ifdef _WIN32
#define Pathname "Path"
#define Sep ';'
#define Slash '\\'
#define Executable(x) !stat(x,&sb)
#else
#define Pathname "PATH"
#define Sep ':'
#define Slash '/'
#define Executable(x) !stat(x,&sb) && Is_executable(myuid, mygid, &sb)
#endif
	char *buf, buf0[4096], *s;
	const char *av0, *p, *p0, *p1, *p2;
	int dotseen, finaldot, rc;
	size_t buflen, L, L1, L2, L3;
	struct stat sb;
	static const char dotslash[3] = { '.', Slash, 0 };
#ifdef _WIN32
	int c;
	pathadj &= 1;	/* . is implicitly in $PATH under _WIN32 */
#else
	gid_t mygid;
	uid_t myuid;
#endif
	if (nidr_exedir_len != -1) {
		fprintf(stderr, "\nIgnoring extra call on nidr_save_argv0()\n");
		return 1;
		}
	nidr_exedir_len = 0;
	if (!(av0 = argv0))
		return 2;
	if (!(p = getenv(Pathname))) {
		fprintf(stderr, "\nnidr_save_exedir: no $%s\n", Pathname);
		return 3;
		}
	dotseen = finaldot = rc = 0;
	buf = buf0;
	buflen = sizeof(buf0);
	p0 = p2 = p;
	p1 = nidr_basename(av0);
	if ((L = p1 - av0) > 0) {
		memcpy(s = (char*)Alloc("nidr_save_argv0", L+1), av0, L);
		s[L] = 0;
		nidr_exedir = s;
		nidr_exedir_len = (int)L;
#ifdef _WIN32
		for(L1 = 0; L1 < L; ++L1)
			if (s[L1] == '/')
				s[L1] = '\\';
#endif
		if (!pathadj)
			return 0;
		if (*p == Sep)
			dotseen = 1;
		while(*p) {
			if (*p == Sep && (p[1] == Sep || p[1] == 0)) {
				dotseen = 1;
				break;
				}
			++p;
			}
		if (s[0] == '.' && s[1] == Slash && L == 2) {
#ifdef _WIN32
			return 0;
#else
			if (!dotseen)
				goto dot_add;
#endif
			}
		L1 = L - 1;
		for(p = p0;;) {
			while(*p == Sep)
				++p;
			if (!*p)
				break;
			if (!strncmp(p, s, L1)) {
#ifdef _WIN32
				return 0;
#else
				if (!(pathadj &= ~1) || *p0 == Sep)
					return 0;
				for(p = p0;;) {
					if (*p == Sep || (*p == '.' && p[1] == Sep))
						return 0;
					while(*p != Sep) {
						if (!*p)
							goto dot_add;
						++p;
						}
					if (!*++p)
						return 0;
					}
#endif
				}
			while(*++p != Sep)
				if (!*p)
					goto break2;
			if (!*++p) {
#ifdef _WIN32
				finaldot = 1; /* leave at 0 for !_WIN32 */
#endif
				break;
				}
			}
 break2:
		L1 = strlen(Pathname);
		L2 = strlen(p0);
		if (!pathadj & 2)
			dotseen = 1;
		L3 = L1 + L2 + L + 3;
		s = (char*)Alloc("nidr_save_argv0", L3);
		memcpy(s, Pathname, L1);
		s[L1++] = '=';
		memcpy(s+L1, p0, L2);
		L1 += L2;
		if (!finaldot)
			s[L1++] = Sep;
		memcpy(s+L1, nidr_exedir, --L);
		L1 += L;
		if (!dotseen)
			s[L1++] = Sep;
		s[L1] = 0;	/* omit final slash */
		putenv(s);
		return 0;
		}
	L = strlen(av0);
#ifdef _WIN32
	if (L < 5 || av0[L-4] != '.'
	|| ((c = av0[L-3]) != 'e' && c != 'E')
	|| ((c = av0[L-2]) != 'x' && c != 'X')
	|| ((c = av0[L-1]) != 'e' && c != 'E')) {
		memcpy(s = (char*)Alloc("nidr_save_argv0", L + 5), av0, L);
		strcpy(s+L, ".exe");
		L += 4;
		av0 = s;
		}
	if (Executable(av0)) {
		/* handle implicit . */
		dotseen = 1;
		nidr_exedir = dotslash;
		}
	else /* do for loop */
#else

/* Fix for C99 */
#ifndef __cplusplus
extern uid_t geteuid(void);
extern uid_t getegid(void);
#endif /* __cplusplus */

	myuid = geteuid();
	mygid = getegid();
#endif
	for(; *p; p = p2) {
		for(p1 = p;; ++p1) {
			if (*p1 == Sep) {
				p2 = p1 + 1;
				if (!*p2)
					finaldot = 1;
				break;
				}
			if (!*p1) {
				p2 = p1;
				break;
				}
			}
		if (p1 == p || (*p == '.' && p1 == p + 1)) {
			if (!dotseen) {
				dotseen = 1;
				if (Executable(av0)) {
					nidr_exedir = dotslash;
					break;
					}
				}
			}
		else {
			L1 = p1 - p;
			L2 = L + L1 + 2;
			if (L2 > buflen) {
				if (buf != buf0)
					free(buf);
				buf = (char*)Alloc("nidr_save_argv0", L2);
				buflen = L2;
				}
			memcpy(buf, p, L1);
			buf[L1++] = Slash;
			strcpy(buf+L1, av0);
			if (Executable(buf)) {
				s = (char*)Alloc("nidr_save_argv0", L1+1);
				memcpy(s, buf, L1);
				s[L1] = 0;
				nidr_exedir = s;
				nidr_exedir_len = (int)L1;
				pathadj &= ~1;
				break;
				}
			}
		}
	if (dotseen)
		pathadj &= ~2;
	if (!finaldot && *p2) {
		while(p2[1])
			++p2;
		if (*p2 == Sep)
			finaldot = 1;
		}
	if (finaldot && !dotseen && !nidr_exedir) {
		dotseen = 1;
		if (Executable(av0))
			nidr_exedir = dotslash;
		}
	if (nidr_exedir == dotslash)
		nidr_exedir_len = 2;
	else {
		if (pathadj & 2 && !finaldot) {
#ifndef _WIN32
 dot_add:
#endif
			L = strlen(p0);
			L1 = strlen(Pathname);
			L2 = L + L1 + 3;
			s = (char*)Alloc("nidr_save_argv0", L2);
			memcpy(s, Pathname, L1);
			s[L1++] = '=';
			memcpy(s+L1, p0, L);
			s[L += L1] = Sep;
			s[L+1] = 0;
			putenv(s);
			}
		if (!nidr_exedir) {
/* Suppress warning for Cygwin and Windows, where path isn't resolved correctly above */
#if !defined(__CYGWIN__) && !defined(_MSC_VER)
			fprintf(stderr, "\nnidr_save_exedir: could not find \"%s\" in $%s\n",
				av0, Pathname);
#endif
			rc = 4;
			}
		}
	if (buf != buf0)
		free(buf);
	return rc;
	}

 void *
nidr_dlopen(const char *libname)
{
#ifdef NO_NIDR_DYNLIB /*{{*/
	botch("dlopen for \"%s\" is NOT SUPPORTED", libname);
	return (void*)libname;
#else /*}{*/
	char buf0[4096], *buf;
	const char *b;
	size_t buflen, L, L1;
	void *h;

	b = nidr_basename(libname);
	if (b > libname)
		return dlopen(libname, RTLD_NOW);
	buf = buf0;
	buflen = sizeof(buf0);
	L = strlen(libname);
	if (nidr_exedir) {
		L1 = L + nidr_exedir_len + 1;
		if (L1 > buflen) {
			buf = (char*)Alloc("nidr_dlopen", L1);
			buflen = L1;
			}
		memcpy(buf, nidr_exedir, nidr_exedir_len);
		strcpy(buf + nidr_exedir_len, libname);
		if ((h = dlopen(buf, RTLD_NOW)))
			goto ret;
		}
	if (!(h = dlopen(libname, RTLD_NOW))) {
		L1 = L + 3;
		if (L1 > buflen) {
			buf = (char*)Alloc("nidr_dlopen", L1);
			buflen = L1;
			}
		buf[0] = '.';
		buf[1] = Slash;
		strcpy(buf+2, libname);
		if (!(h = dlopen(buf, RTLD_NOW)))
			h = dlopen(libname, RTLD_NOW); 	/* for dlerror */
		}
 ret:
	if (buf != buf0)
		free(buf);
	return h;
#endif  /*}}*/
	}

#undef Executable
#undef Slash
#undef Sep
#undef Pathname


 struct
Comment {
	int k;		/* subscript for comfree */
	size_t avail;	/* bytes left (from tnext) */
	char *text;	/* text of comment */
	char *tnext;	/* null byte at end of comment */
	Comment *fnext;	/* next free Comment */
	};

 enum { Comment_kmax = 7 };

 static Comment *comfree[Comment_kmax+1];
 static size_t Comment_maxlen[Comment_kmax+1];

 static void
comment_free(Comment *c)
{
	int k = c->k;

	if (k > Comment_kmax)
		free(c);
	else {
		c->fnext = comfree[k];
		comfree[k] = c;
		}
	}

 static Comment*
alloc_comment(int k, size_t L)
{
	Comment *c;

	for(; k <= Comment_kmax; ++k) {
		if (L <= Comment_maxlen[k]) {
			L = Comment_maxlen[k];
			if ((c = comfree[k])) {
				comfree[k] = c->fnext;
				goto have_c;
				}
			break;
			}
		}
	c = (Comment*)Alloc("save_comment", L + sizeof(Comment) + 1);
	c->k = k;
	c->text = (char*)(c+1);
 have_c:
	c->avail = L;
	c->tnext = c->text;
	return c;
	}

 static void
save_comment(const char *s)
{
	Comment *c, *c1, **cp;
	size_t L, L1;

	L = strlen(s);
	cp = curid ? &curid->comment : curkws ? &curkws->comment : &OutsideComment;
	if ((c = *cp)) {
		if (c->avail >= L)
			goto cupdate;
		L1 = c->tnext - c->text;
		c1 = alloc_comment(c->k + 1, L + L1);
		memcpy(c1->text, c->text, L1);
		c1->tnext = c1->text + L1;
		c1->avail -= L1;
		comment_free(c);
		c = c1;
		}
 	else
		c = alloc_comment(0, L);
 cupdate:
	memcpy(c->tnext, s, L+1);
	c->tnext += L;
	c->avail -= L;
	*cp = c;
	}

 static void
comment_setup(void)
{
	int i;
	size_t L;
	nidr_comment = save_comment;
	/* "- 1" to allow for terminating '\0' */
	for(L = 64; L <= sizeof(Comment) - 1; L <<= 1);
	for(i = 0; i <= Comment_kmax; ++i, L <<= 1)
		Comment_maxlen[i] = L - sizeof(Comment) - 1;
	}

 static void
comment_reset(void)
{
	Comment *c, *c1;
	int i;

	for(i = 0; i <= Comment_kmax; ++i) {
		c1 = comfree[i];
		comfree[i] = 0;
		while((c = c1)) {
			c1 = c->fnext;
			free(c);
			}
		}
	nidr_comment = 0;
	}

 static void
dumpcomment(Comment **cp)
{
	Comment *c = *cp;
	*cp = 0;
	fprintf(dumpfile, "%s", c->text);
	comment_free(c);
	}

 static void
dumpname(int hasval, KeyWord *kw)
{
	const char *fmt[2] = { "%s", "%s =" };
	int i;
	if (OutsideComment)
		dumpcomment(&OutsideComment);
	if (primary)
		kw += kw->paoff;
	for(i = 0; i < dumplev; ++i)
		putc(' ', dumpfile);
	fprintf(dumpfile,fmt[hasval],kw->name);
	if (!hasval) {
		if (kw->comment)
			dumpcomment(&kw->comment);
		else if (kw != curkw)
			putc('\n', dumpfile);
		}
	}

 static void
dumpstring(const char *s0)
{
	const char *s;
	int c, n1, n2, q;

	n1 = n2 = 0;
	for(s = s0;;)
		switch(*s++) {
		  case 0: goto break2;
		  case '\'':
			++n1;
			break;
		  case '"':
			++n2;
		  }
 break2:
	q = '\'';
	if (n1 > n2)
		q = '"';
	putc(' ', dumpfile);
	putc(q, dumpfile);
	s = s0;
	while((c = *s++)) {
		if (c == q)
			putc(q, dumpfile);
		putc(c, dumpfile);
		}
	putc(q, dumpfile);
	}

 static void
dumpvals0(KeyWord *kw)
{
	Real *r;
	const char **sp;
	int i, *ip, indent, j, n;

	ip = 0; /* shut up warning of possible use without initialization */
	sp = 0; /* ditto */
	if (!(r = KWval.r) && !(ip = KWval.i) && !(sp = KWval.s))
		return;
	n = KWval.n;
	putc((indent = n > 1) ? '\n' : ' ', dumpfile);
	for(i = 0;;) {
		if (indent) {
			putc('\t', dumpfile);
			for(j = 0; j < dumplev; ++j)
				putc(' ', dumpfile);
			}
		if (r)
			fprintf(dumpfile, "%.15g", r[i]);
		else if (ip)
			fprintf(dumpfile, "%d", ip[i]);
		else
			dumpstring(sp[i]);
		if (++i >= n)
			break;
		indent = 1;
		putc('\n', dumpfile);
		}
	if (kw->comment)
		dumpcomment(&kw->comment);
	else
		putc('\n', dumpfile);
	}

 static void (*dumpvals)(KeyWord *kw) = dumpvals0;

 static void
dumpvals1(KeyWord *kw)
{
	Real *r;
	const char **sp;
	int i, *ip, n;

	ip = 0; /* shut up warning of possible use without initialization */
	sp = 0; /* ditto */
	if ((r = KWval.r) || (ip = KWval.i) || (sp = KWval.s)) {
		n = KWval.n;
		for(i = 0; i < n; ++i) {
			if (r)
				fprintf(dumpfile, " %.15g", r[i]);
			else if (ip)
				fprintf(dumpfile, " %d", ip[i]);
			else
				dumpstring(sp[i]);
			}
		}
	if (kw->comment)
		dumpcomment(&kw->comment);
	else
		putc('\n', dumpfile);
	}

 char *
nidr_KWscopy(const char *s)
{
	Sbuf *sb;
	char *rv;

	size_t L = strlen(s) + 1;
	if (L >= n_KWsbuf)
		botch("String too long in KWscopy");
	if (KWsbufe - KWsbuf1 < L) {
		if (!KWsbuf->next) {
			KWsbuf->next = sb = (Sbuf*)Alloc("KWscopy", sizeof(Sbuf));
			sb->next = 0;
			}
		KWsbuf = KWsbuf->next;
		KWsbuf1 = KWsbuf->buf;
		KWsbufe = KWsbuf1 + n_KWsbuf;
		}
	strcpy(KWsbuf1, s);
	rv = KWsbuf1;
	KWsbuf1 += L;
	return rv;
	}

 static void
KWvalbuf_inc(void)
{
	Real *r;
	Uint n;

	n = nKWvalbuf << 1;
	r = (Real*)Alloc("KWvalbuf", n*sizeof(Real));
	memcpy(r, KWvalbuf, nKWvalbuf*sizeof(Real));
	free(KWvalbuf);
	KWvalbuf = r;
	nKWvalbuf = n;
	KWvalmax.n <<= 1;
	if (KWval.r) {
		KWval.r = r;
		KWvalmax.r = r + n;
		}
	else if (KWval.i) {
		KWval.i = (int*) r;
		KWvalmax.i = (int*)(r + n);
		}
	else if (KWval.s) {
		KWval.s = (const char**)r;
		KWvalmax.s = (const char**)(r + n);
		}
	else
		botch("Unexpected case in KWvalbuf_inc");
	}

/* KWval.rstate values...
 *	value	form seen
 *	0	v
 *	1	L:u
 *	2	L:s:u
 *	3	n*v
 *	4	n*L:u
 *	5	n*L:s:u
 */

 static void
finish_rexpand(void)
{
	int i, k, n, os;
	Real sgn, st, u, v, x;

	os = KWval.rstate;
	KWval.rstate = 0;
	n = KWval.n;
	k = 1;
	if (os >= 3) {
		KWval.n = n -= os-1;
		k = KWval.r[n];
		if (k != KWval.r[n]) {
			squawk("Noninteger replication factor %.17g", KWval.r[n]);
			return;
			}
		else if (k < 1) {
			squawk("Nonpositive replication factor %d", k);
			return;
			}
		++n;
		os -= 3;
		}
	else
		KWval.n = n -= os + 1;
	v = KWval.r[n++];
	u = st = 0.;	/* Shut up warning of not being initialized. */
			/* Both will be assigned before being used. */
	switch(os) {
	  case 0:
		n = KWval.n;
		for(i = 0; i < k; ++i) {
			if (n >= KWvalmax.n)
				KWvalbuf_inc();
			KWval.r[n++] = v;
			}
		KWval.n = n;
		return;
	  case 1:
		st = 1;
		u = KWval.r[n];
		break;
	  case 2:
		st = KWval.r[n];
		if (st == 0.) {
			squawk("Invalid stride == zero.");
			return;
			}
		u = KWval.r[n+1];
	  }
	sgn = 1.;
	if (st < 0.)
		sgn = -1.;
	if (sgn*(u - v) < 0.) {
		squawk("Empty sequence.");
		return;
		}
	n = KWval.n;
	do {
		for(i = 0; sgn*(u - (x = v + i*st)) >= 0.; ++i) {
			if (n >= KWvalmax.n)
				KWvalbuf_inc();
			KWval.r[n++] = x;
			}
		}
		while(--k > 0);
	KWval.n = n;
	}

 static void
rexpand(int state)
{
	int os;

	os = KWval.rstate;
	KWval.rstate = 0;
	switch(state) {
	  case 1: /* just saw *v */
		if (os == 0)
			KWval.rstate = 3;
		else
			squawk("Unexpected '*'");
		break;
	  case 2: /* just saw :v */
		if (os == 2 || os == 5)
			squawk("Unexpected ':'");
		else
			KWval.rstate = os + 1;
		break;
	  }
	}

 static void
nidr_bufr_strict(Real r, int state)
{
	int k, n;

	if (KWval.s) {
		squawk("expected a quoted string, but found a number");
		return;
		}
	if (KWval.rstate && !state)
		finish_rexpand();
	if (!KWval.r && !KWval.i) {
		squawk("No values may be specified for %s", KWStack->kw->name);
		return;
		}
	if ((n = KWval.n) >= KWvalmax.n)
		KWvalbuf_inc();
	if (KWval.r)
		KWval.r[n] = r;
	else {
		k = (int)r;
		if (k != r)
			squawk("truncating %.17g to %d", r, k);
		KWval.i[n] = k;
		}
	++KWval.n;
	if (state | KWval.rstate)
		rexpand(state);
	}

 static void
nidr_bufs_strict(const char *s)
{
	if (!KWval.s) {
		if (KWval.r)
			squawk("Expected a number, but found a quoted string");
		else
			squawk("Misplaced quoted string");
		return;
		}
	if (KWval.n >= KWvalmax.n)
		KWvalbuf_inc();
	KWval.s[KWval.n++] = s;
	}

 void
nidr_reset(void)
{
	/* Originally did this in case KWKind_Str of kw_setup(), */
	/* but this leads to confusion with erroneous input. */
	if (curkw)
		nidr_keyword_finish();
	KWsbuf = &KWsbuf0;
	KWsbuf1 = KWsbuf0.buf;
	KWsbufe = KWsbuf0.buf + n_KWsbuf;
	}

 NIDR_KWlib *
nidr_lib_record(void *h, const char *libname)
{
	NIDR_KWlib *Lib;
	size_t L;

	L = strlen(libname) + 1;
	Lib = (NIDR_KWlib*)Alloc("NIDR_lib_record", sizeof(NIDR_KWlib) + L);
	memset(Lib, 0, sizeof(NIDR_KWlib));
	memcpy(Lib->libname = (char*)(Lib+1), libname, L);
	if (!(Lib->next = NIDR_Libs))
		atexit(nidr_lib_cleanup);
	NIDR_Libs = Lib;
	Lib->h = h;
	return Lib;
	}

 static KeyWord*
kw_insert(KeyWord *kw, int *tryagain)
{
#ifdef NO_NIDR_DYNLIB /*{{*/
	botch("Loading library \"%s\" for %s is disallowed",
		kw->f.vf, kw->name);
#else /*}{*/
	KW_ADD kwa;
	KeyWord *kw0, *kw1, *kw2;
	NIDR_KWlib *Lib;
	Uint u1, ui;
	const char *lname, *s;
	int newtop, nmatch;
	void *h;

	if (tryagain)
		*tryagain = 0;
	if (kw->kind & KWKind_Loaded)
		return (KeyWord*)kw->f.vs;
	h = nidr_dlopen(lname = (const char*)kw->f.vf);
	if (!h) {
#ifndef NO_DLERROR
		if ((s = dlerror()))
			botch("Cannot open library \"%s\" for %s:\n\t%s",
				lname, kw->name, s);
		else
#endif
			botch("Cannot open library \"%s\" for %s",
				lname, kw->name);
		}
	if (!find_dlsym(kwa, h, "keyword_add"))
		botch("dlsym(\"keyword_add\") failed for %s", lname);
	kw1 = (*kwa)();
	if (!(s = kw1->name)) {
		s = "<NULL>";
		goto namebotch;
		}
	newtop = 0;
	if (strcmp(s, kw->name)) {
		if (!KW_cur && !strcmp(s,"KeywordTop") && kw1->kind & KWKind_Dynmult)
			newtop = 1;
		else if (tryagain && kw1->nkw > 0
		 && (kw2 = kwfind(kw->name, kw1->kw, kw1->nkw, &nmatch))) {
			if (nmatch > 1) {
				toomany(kw->name, kw1, nmatch);
				botch("Too many matches in library %s", lname);
				}
			*tryagain = 1;
			}
		else
 namebotch:
			botch("Library %s: expected top keyword to be %s but got %s",
				lname, kw->name, s);
		}
	ui = kw->kind  & (KWKind_Mask|KWKind_List);
	u1 = kw1->kind & (KWKind_Mask|KWKind_List);
	if (ui != u1)
		botch("Library %s: expected kind %u for %s, but got %u",
			lname, ui, s, u1);
	Lib = nidr_lib_record(h, lname);
	Lib->kw0 = kw0 = kw;
	memcpy(&Lib->kw, kw, sizeof(KeyWord));
	kw = &Lib->kw;
	kw->kw = kw1->kw;
	kw->nkw = kw1->nkw;
	kw->f = kw1->f;
	kw0->f.vs = (void*)kw;
	kw0->kind |= KWKind_Loaded;
	if (newtop) {
		Lib->oldtop = Keyword_Top;
		Keyword_Top = kw;
		kw->kind |= KWKind_Dynmult;
		}
#endif	/*}}*/
	return kw;
	}

 static void
kwnext_setup(KeyWord *kw, Uint n)
{
	KeyWord *kwe;

	if (kw->kwnext || (n <= 1 && kw->name))
		return;
#ifndef NO_NIDR_DYNLIB /*{*/
	if (kw->kind & KWKind_Extended) {
		KeyWordx *kx1, *kxe;
		for(kx1 = (KeyWordx*)kw; !kx1->kw.name; ++kx1)
			kx1->kw.kwnext = (KeyWord*)(kx1 + 1);
		for(kxe = kx1 + n - 1; kx1 < kxe; ++kx1)
			kx1->kw.kwnext = (KeyWord*)(kx1 + 1);
		return;
		}
#endif
	for(; !kw->name; ++kw)
		kw->kwnext = kw + 1;
	for(kwe = kw + n - 1; kw < kwe; ++kw)
		kw->kwnext = kw + 1;
	}

 static void
KWStack_inc(void)
{
	KWinfo *kwi;
	Uint nn;
	size_t len;

	nn = n_KWStack << 1;
	kwi = (KWinfo*)Alloc("kw_setup", len = nn*sizeof(KWinfo));
	memcpy(kwi, KWStackBot, len >> 1);
	if (KWStackBot != KWStack0)
		free(KWStackBot);
	KWStackBot = kwi;
	KWStackEnd = kwi + nn;
	KWStack = kwi + n_KWStack;
	n_KWStack = nn;
	}

 static KeyWord*
kw_setup(KeyWord *kw, void *g, const char *name)
{
	KWinfo *kwi;
	KeyWord **alt, *kw1, **req;

	Uint k, nalt, nn, nreq;
	int *altct, deferred;
	size_t len;

	deferred = 0;
	if (kw->kind & KWKind_Dynlib) {
		if (kw->kw)
			deferred = 1;
		else
			kw = kw_insert(kw, 0);
		}
 top:
	if ((kw1 = kw->kw)) {
		kwnext_setup(kw1, kw->nkw);
		if (kw->kind & KWKind_Dynmult)
			return kw;
		while(!kw1->name) {
			if (!(kw1->kind & KWKind_Stacked)) {
				kw1->kind |= KWKind_Stacked;
				kw_setup(kw1, g, name);
				}
			kw1 = kw1->kwnext;
			}
		}
	if (!curkw) {
		KWStack = KWStackBot = KWStack0;
		KWStackEnd = KWStack0 + n_KWStack0;
		curkw = kw;
		}
	else if (++KWStack >= KWStackEnd)
		KWStack_inc();
	kwi = KWStack;
	kwi->name = name;
	kwi->kw = kw;
	kwi->kw1 = kw1;
	nalt = nreq = 0;
	if (kw1)
		for(; kw1; kw1 = kw1->kwnext) {
			if (nalt < kw1->alt)
				nalt = kw1->alt;
			if (nreq < kw1->req)
				nreq = kw1->req;
			}
	kwi->nalt = nalt;
	kwi->nreq = nreq;
	alt = req = 0;
	altct = 0;
	if ((nn = nalt + nreq) > 0) {
		nn += 2;
		alt = (KeyWord**)Alloc("kw_setup(alt)",
				len = nn*sizeof(KeyWord*) + (nalt+1)*sizeof(int));
		memset(alt, 0, len);
		req = alt + nalt + 1;
		altct = (int*)(req + nreq + 1);
		/* altct[0], alt[0] and req[0] = "don't care" slots */
		}
	kwi->alt = alt;
	kwi->req = req;
	kwi->altct = altct;
	if (nreq)
		for(kw1 = kwi->kw1; kw1; kw1 = kw1->kwnext)
			req[kw1->req] = kw1;
	if (nalt)
		for(kw1 = kwi->kw1; kw1; kw1 = kw1->kwnext)
			if (kw1->kind & KWKind_primary)
				++altct[kw1->alt];
	kwi->g = g;
	KWval.n = 0;
	KWval.i = 0;
	KWval.r = 0;
	KWval.s = 0;
	if ((k = kw->kind & KWKind_Mask)) {
		if (!KWvalmax.r)
			KWvalbuf = (Real *)Alloc("kw_setup(KWvalbuf)",
						(nKWvalbuf = 128)*sizeof(Real));
		switch(k) {

		  case KWKind_Int:
			KWval.i = (int*)KWvalbuf;
			KWvalmax.n = (nKWvalbuf*sizeof(Real))/sizeof(int);
			KWvalmax.i = KWval.i + KWvalmax.n;
			break;

		  case KWKind_Real:
			KWval.r = KWvalbuf;
			KWvalmax.r = KWvalbuf + (KWvalmax.n = nKWvalbuf);
			break;

		  case KWKind_Str:
			KWval.s = (const char**)KWvalbuf;
			KWvalmax.n = (nKWvalbuf*sizeof(Real))/sizeof(char*);
			KWvalmax.s = KWval.s + KWvalmax.n;
		  }
		}
	if (deferred) {
		kw = kw_insert(kw, 0);
		deferred = 0;
		goto top;
		}
	if (!(kwi->needstart = kw->kind & KWKind_Mask)) {
		if (kw->name) {
			if (dumpfile)
				dumpname(0, kw);
			++dumplev;
			}
		if (kw->f.start)
			(*kw->f.start)(kw->name, 0, &KWStack->g, kw->f.vs);
		}
	else if (!kw->f.start && NIDR_disallow_missing_start)
		botch("No start routine for %s", kw->name);
	return kw;
	}

 static KeyWord *
kwfind(const char *name, KeyWord *keywds, int n, int *nmatch)
{
	KeyWord *kn, *kn1;
	int k, n0, n1, n2, nn;
	size_t L;

	*nmatch = 0;
	if (n <= 0)
		return 0;
	L = strlen(name);
	n0 = 0;
	nn = n;
#ifndef NO_NIDR_DYNLIB /*{*/
	if (n > 0 && keywds->kind & KWKind_Extended) {
	    while(n > 0) {
		n1 = n >> 1;
		kn = (KeyWord*)((KeyWordx*)keywds + n1);
		k = strncmp(name, kn->name, L);
		if (k < 0)
			n = n1;
		else if (k > 0) {
			n -= ++n1;
			n0 += n1;
			keywds = (KeyWord*)((KeyWordx*)kn + 1);
			}
		else {
			/* Found -- check for range of matches. */
			/* Here we use linear search, as we expect */
			/* the range to be small. */
			n = n1 + n0;
			n2 = n + 1;
			if (kn->name[L]) {
				for(kn1 = kn; n2 < nn; ++n2) {
					kn1 = (KeyWord*)((KeyWordx*)kn1 + 1);
					if (strncmp(name, kn1->name, L))
						break;
					if (!kn1->name[L])
						goto found1;
					}
				kn1 = kn;
				while(n > 0) {
					kn1 = (KeyWord*)((KeyWordx*)kn1 - 1);
					if (strncmp(name, kn1->name, L))
						break;
					if (!kn1->name[L])
						goto found1;
					--n;
					kn = kn1;
					}
				}
			*nmatch = n2 - n;
			return kn;
			}
		}
	    }
	else
#endif	/*}*/
	while(n > 0) {
		n1 = n >> 1;
		kn = keywds + n1;
		k = strncmp(name, kn->name, L);
		if (k < 0)
			n = n1;
		else if (k > 0) {
			n -= ++n1;
			n0 += n1;
			keywds = kn + 1;
			}
		else {
			/* Found -- check for range of matches. */
			/* Here we use linear search, as we expect */
			/* the range to be small. */
			n = n1 + n0;
			n2 = n + 1;
			if (kn->name[L]) {
				for(kn1 = kn; n2 < nn; ++n2) {
					++kn1;
					if (strncmp(name, kn1->name, L))
						break;
					if (!kn1->name[L])
						goto found1;
					}
				kn1 = kn;
				while(n > 0) {
					--kn1;
					if (strncmp(name, kn1->name, L))
						break;
					if (!kn1->name[L]) {
 found1:
						*nmatch = 1;
						return kn1;
						}
					--n;
					kn = kn1;
					}
				}
			*nmatch = n2 - n;
			return kn;
			}
		}
	return 0;	/* not found */
	}

 static KeyWord *
toomany(const char *name, KeyWord *kw, int nmatch)
{
	int i;
	squawk("\"%s\" is ambiguous; possible matches..", name);
	if (nsquawk <=  NIDR_SQUAWKMAX)
		for(i = 0; i < nmatch; i++, kw++)
			fprintf(stderr, "\t%s\n", kw->name);
	return 0;
	}

 KeyWord *
nidr_keyword(const char *name)
{
	int nmatch;
	KeyWord *kw, *kw1;

	kw = kwfind(name, Keyword_Top->kw, Keyword_Top->nkw, &nmatch);
	if (nmatch > 1)
		return toomany(name, kw, nmatch);
	else if (kw) {
		if (!(kw1 = curkw)) {
			kw = kw_setup(kw, KW_g, name);
			if (kw->kind & KWKind_Dynmult)
				return kw;
			}
		if (!strict) {
			if (kw1)
				nidr_keyword_finish();
			kw_setup1(kw);
			}
		}
	return kw;
	}

 static void
valcheck(KeyWord *kw)
{
	Real *r;
	int *z;
	int i, k, n;

	n = KWval.n;
	switch(k = kw->kind & KWKind_Mask) {
	  case KWKind_Int:
		z = KWval.i;
		if (kw->kind & KWKind_strictLb) {
			for(i = 0; i < n; ++i)
				if (z[i] <= kw->Lb) {
					squawk("%s must be > %.0f", kw->name, kw->Lb);
					break;
					}
				}
		else if (kw->kind & KWKind_caneqLb) {
			for(i = 0; i < n; ++i)
				if (z[i] < kw->Lb) {
					squawk("%s must be >= %.0f", kw->name, kw->Lb);
					break;
					}
				}
		if (kw->kind & KWKind_strictUb) {
			for(i = 0; i < n; ++i)
				if (z[i] >= kw->Ub) {
					squawk("%s must be < %.0f", kw->name, kw->Ub);
					break;
					}
				}
		else if (kw->kind & KWKind_caneqUb) {
			for(i = 0; i < n; ++i)
				if (z[i] > kw->Ub) {
					squawk("%s must be >= %.0f", kw->name, kw->Ub);
					break;
					}
				}
		break;
	  case KWKind_Real:
		r = KWval.r;
		if (kw->kind & KWKind_strictLb) {
			for(i = 0; i < n; ++i)
				if (r[i] <= kw->Lb) {
					squawk("%s must be > %g", kw->name, kw->Lb);
					break;
					}
				}
		else if (kw->kind & KWKind_caneqLb) {
			for(i = 0; i < n; ++i)
				if (r[i] < kw->Lb) {
					squawk("%s must be >= %g", kw->name, kw->Lb);
					break;
					}
				}
		if (kw->kind & KWKind_strictUb) {
			for(i = 0; i < n; ++i)
				if (r[i] >= kw->Ub) {
					squawk("%s must be < %g", kw->name, kw->Ub);
					break;
					}
				}
		else if (kw->kind & KWKind_caneqUb) {
			for(i = 0; i < n; ++i)
				if (r[i] > kw->Ub) {
					squawk("%s must be >= %g", kw->name, kw->Ub);
					break;
					}
				}
		break;
	  default:
		botch("Bug: unexpected (kw->kind & KWKind_Mask) = %d in valcheck",n);
	  }
	}

 static void
read_lib(const char *libname, KeyWord *kw)
{
#ifdef NO_NIDR_DYNLIB /*{{*/
	botch("LIBNAME is disallowed: cannot read \"%s\"", libname);
#else /*}{*/
	KeyWord *kw1;
	KW_ADD kwa;
	NIDR_KWlib *Lib;
	void *h;

	h = nidr_dlopen(libname);
	if (!h) {
#ifndef NO_DLERROR
		const char *s;
		if ((s = dlerror()))
			botch("Cannot open library \"%s\":\n\t%s", libname, s);
		else
#endif
			botch("Cannot open library \"%s\"", libname);
		}
	if (!find_dlsym(kwa, h, "keyword_add"))
		botch("dlsym(\"keyword_add\") failed for %s", libname);
	kw1 = (*kwa)();
	Lib = nidr_lib_record(h, libname);
	kw->nkw = kw1->nkw;
	kw->kw = kw1->kw;
	kw->f = kw1->f;
	kw->kind |= KWKind_Loaded;
#endif	/*}}*/
	}

 static void
nidr_id_strict_finish(KWinfo *kwi, KeyWord *kw, const char *name)
{
	KeyWord *kw1;
	int n;

	if (kw->alt) {
		if ((kw1 = kwi->alt[n = kw->alt])) {
			if (strcmp(kw1->name, name))
				squawk("%s and %s are mutually exclusive",
					kw1->name, name);
			else
				squawk("%s was already specified", name);
			}
		else
			kwi->alt[n] = kw;
		}
	if (kw->req) {
		if (kwi->req[n = kw->req])
			kwi->req[n] = 0;
		else if (!kw->alt)
			squawk("%s specified more than once", name);
		}
	}

 static KWinfo *
dispatch_val(KWinfo *kwi)
{
	KeyWord *kw = kwi->kw;

	kwi->needstart = 0;
	if (KWval.n) {
		if (KWval.rstate)
			finish_rexpand();
		if (dumpfile) {
			dumpname(1, kw);
			dumpvals(kw);
			}
		if (kw->kind & (KWKind_Lb|KWKind_Ub))
			valcheck(kw);
		if (kw->f.start)
			(*kw->f.start)(kw->name, &KWval, &kwi->g, kw->f.vs);
		else if ((kw->kind & (KWKind_Libname | KWKind_Loaded)) == KWKind_Libname) {
			read_lib(KWval.s[0], kw);
			kw = kw_setup(kw, kwi->g, kw->name);
			if (kw->f.start)
				(*kw->f.start)(kw->name, &KWval, &kwi->g, kw->f.vs);
			if (kw == kwi->kw) {
				*kwi = *KWStack;
				--KWStack;
				}
			else
				kwi = KWStack;
			}
		KWval.n = 0;
		}
	else if ((kw->kind & (KWKind_Libname|KWKind_Loaded))
			  != (KWKind_Libname|KWKind_Loaded))
		squawk("expected %sone %s value for %s",
			kw->kind & KWKind_List ? "at least " : "",
			valkind[(kw->kind & KWKind_Mask)-1], kw->name);
	++dumplev;
	return kwi;
	}

 static void
oneof(KeyWord *kw, int alt, int n)
{
	KeyWord *kw1;

	squawk("One of the following %d entities\nmust be specified for %s..",
		n, kw->name);
	for(kw1 = kw->kw; !kw1->name; kw1 = kw1->kwnext);
	for(; kw1; kw1 = kw1->kwnext)
		if (kw1->alt == alt && kw1->kind & KWKind_primary)
			fprintf(stderr, "\t%s\n", kw1->name);
	}

 static void
missing_chk(KeyWord *kw1, KWinfo *kwi)
{
	KeyWord *kw0, *kw2, **req;
	Uint a;
	char seen0[1024], *seen;
	const char *kwname;
	int n;
	size_t nreq;

	/* only issue one error message per missing keyword */

	nreq = 0;
	for(kw0 = kw1; kw1; kw1 = kw1->kwnext)
		if (nreq < kw1->req)
			nreq = kw1->req;
	seen = seen0;
	if (++nreq > sizeof(seen0))
		seen = (char*)Alloc("missing_chk", nreq);
	memset(seen, 0, nreq);
	req = kwi->req;
	for(kw1 = kw0; kw1; kw1 = kw1->kwnext) {
		if (kw1->req && req[kw1->req] && !seen[kw1->req] && kw1->kind & KWKind_primary) {
			seen[kw1->req] = 1;
			a = -1;
			if (!kw1->alt || (n = kwi->altct[a = kw1->alt]) <= 1) {
				if (!(kwname = kwi->name))
					kwname = "<NIDRBUG>";
				for(kw2 = kw1;;) {
					if (kw2->alt == a && kw2->kind & KWKind_primary)
						break;
					if (!(kw2 = kw2->kwnext))
						botch("Bug in missing_chk");
					}
				squawk("%s must be specified for %s",
					kw2->name, kwname);
				}
			else
				oneof(kwi->kw, kw1->alt, n);
			}
		}
	if (seen != seen0)
		free(seen);
	}

 static void
finalize(KWinfo *kwi)
{
	KeyWord *kw, *kw1, **req;

	kw = kwi->kw;
	kw->kind &= ~KWKind_Stacked;
	if (kwi->needstart)
		kwi = dispatch_val(kwi);
	if (kw->name)
		--dumplev;
	if (kw->f.final)
		(*kw->f.final)(kw->name, 0, &kwi->g, kw->f.vf);
	if (kwi->alt) {
		if (kwi->nreq) {
			req = kwi->req;
			for(kw1 = kwi->kw1; kw1; kw1 = kw1->kwnext)
				if (kw1->req && req[kw1->req]) {
					missing_chk(kw1, kwi);
					break;
					}
			}
		free(kwi->alt);
		}
	}

 static KeyWord *
nidr_identifier_strict(const char *name)
{
	KWinfo *kwi, *kwi1;
	KeyWord *kw, *kw1;
	int nmatch;
	size_t height;

	if (!curkw)
		botch("curkw = 0 in nidr_identifier");
	kwi = KWStack;
	if (kwi->needstart)
		kwi = dispatch_val(kwi);
	for(kwi1 = kwi;;) {
		kw1 = kwi->kw;
		if ((kw = kwfind(name, kwi->kw1, kw1->nkw, &nmatch)))
			break;
		if (kwi == KWStackBot)
			return 0;
		if ((--kwi)->kw->name && !(kwi->kw->kind & KWKind_Loaded))
			kwi1 = kwi;
		}
	if (nmatch > 1)
		return toomany(name, kw, nmatch);
	while(KWStack > kwi1)
		finalize(KWStack--);
	if ((kw->kind & (KWKind_Libname | KWKind_Loaded)) == KWKind_Libname) {
		nidr_id_strict_finish(kwi, kw, name);
		if (!KWvalmax.r)
			KWvalbuf = (Real *)Alloc("nidr_identifier_strict",
						(nKWvalbuf = 128)*sizeof(Real));
		if (++KWStack >= KWStackEnd)
			KWStack_inc();
		kwi = KWStack;
		kwi->kw = kw;
		kwi->needstart = 1;
		KWval.s = (const char**)KWvalbuf;
		KWvalmax.n = (nKWvalbuf*sizeof(Real))/sizeof(char*);
		KWvalmax.s = KWval.s + KWvalmax.n;
		}
	else {
		height = kwi - KWStackBot;
		kw = kw_setup(kw, kwi->g, name);
		kwi = KWStackBot + height; /* in case kw_setup reallocated KWStack */
		nidr_id_strict_finish(kwi, kw, name);
		}
	return kw;
	}

 static void
nidr_keyword_finish(void)
{
	if (!strict)
		kw_finish2();
	for(;;--KWStack) {
		finalize(KWStack);
		if (KWStack == KWStackBot)
			break;
		}
	if (!strict)
		kw_finish3();
	curid = curkw = 0;
	}

 const char*
nidr_keyword_name(void)
{ return curkw ? curkw->name : "<none>"; }

/* Some of the above assumes strict nesting according to dakota.input.nspec. */
/* Code here is meant to relax this assumption, allowing more flexibility in */
/* the order of identifiers within a DAKOTA "keyword". */

 typedef struct KWpair KWpair;
 typedef struct KWmblk KWmblk;

 struct
KWmblk {
	KWmblk *next;
	size_t len;
	/* memory of length len immediately follows */
	};

 struct
KWpair {
	KeyWord *kw;
	KWseen *kws;
	};

 enum{ KWmblk_gulp = 32000 };

 static AVL_Tree *AVLT, *AVLKWP;
 static KWseen **KW_p, **KW_pe, KWmissing, *KWs0;
 static KWmblk *KWmblk0, *KWmblk1;
 static const char *KWmem0, *KWmem1;

 typedef struct
AVLCmpInfo {
	KWseen **found[2];
	int nfound;
	int inexact;
	} AVLCmpInfo;

 static int
avlcmp(void *v, KWseen **a, KWseen **b)
{
	AVLCmpInfo *AI = (AVLCmpInfo*)v;
	KWseen *ksa, *ksb;
	const char *s, *t;

	s = (ksa = *a)->name;
	t = (ksb = *b)->name;
	for(; *s == *t; ++s, ++t)
		if (!*s)
			return 0;
	if ((!*s && !ksa->kw && ksb->kw)
	  ||(!*t && !ksb->kw && ksa->kw)) {
		/* inexact match */
		if (AI->nfound == 0
		|| (AI->nfound == 1 && AI->found[0] != b))
			AI->found[AI->nfound++] = b;
		return AI->inexact;
		}
	return *s - *t;
	}

 static int
kwpcmp(void *v, KWpair *a, KWpair *b)
{
	if (a->kw == b->kw)
		return 0;
	return a->kw > b->kw ? 1 : -1;
	}

 static void
KWmeminit(void)
{
	KWmblk0 = KWmblk1 = (KWmblk*)Alloc("KWmeminit",
			sizeof(KWmblk) + KWmblk_gulp);
	KWmem0 = (char*)(KWmblk0 + 1);
	KWmem1 = KWmem0 + KWmblk_gulp;
	KWmblk0->len = KWmblk_gulp;
	KWmblk0->next = 0;
	KWmissing.mnext = KWmissing.mprev = &KWmissing;
	KW_cur = 0;
	memset(&KWval, 0, sizeof(KWval));
	KWvalbuf = (Real *)Alloc("kw_setup(KWValbuf)", (nKWvalbuf = 128)*sizeof(Real));
	ToClear = ToClear0 = (KeyWord**)Alloc("kw_setup(ToClear)", 256*sizeof(KeyWord*));
	ToClearEnd = ToClear0 + 256;
	}

 static void
KWmembump(size_t L)
{
	KWmblk *mb, *mb1;
	size_t L1;

	for(L1 = KWmblk_gulp; L1 < L; L1 <<= 1);
	if ((mb = mb1 = KWmblk1->next) && L1 <= mb->len)
		L1 = mb->len;
	else {
		KWmblk1->next = mb = (KWmblk*)Alloc("KWmembump", L1 + sizeof(KWmblk));
		mb->len = L1;
		mb->next = mb1;
		}
	KWmblk1 = mb;
	KWmem0 = (char*)(mb+1);
	KWmem1 = KWmem0 + L1;
	}

 static void *
KWgetmem(size_t L)	/* for aligned memory */
{
	void *rv;

	L = (L + sizeof(Real) - 1) & ~(sizeof(Real) - 1);
	if (KWmem1 - KWmem0 < L)
		KWmembump(L);
	rv = (void*)KWmem0;
	KWmem0 += L;
	return rv;
	}

 static KWseen **
KWhash(const char *s, KeyWord *kw)
{
	AVLCmpInfo AI;
	KWseen KW0, *KW0p, *kws, **kwsp;
	char **ps;
	const char *sa, *sb;

	AI.nfound = 0;
	AI.inexact = -1;
	AVL_setv(AVLT, &AI);
	KW0.name = s;
	KW0.kw = kw;
	KW0p = &KW0;
	curkws = 0;
	if ((kwsp = (KWseen**)AVL_find((const Element*)&KW0p, AVLT)))
		return kwsp;
	if (AI.nfound) {
		if (AI.nfound == 1) {
			AI.inexact = 1;
			AVL_find((const Element*)&KW0p, AVLT);
			if (AI.nfound == 1) {
				if (kw && (kw->kind & (KWKind_Libname | KWKind_Loaded))
						== KWKind_Libname
				 && (ps = (*AI.found[0])->svals))
					read_lib(ps[0], kw);
				return AI.found[0];
				}
			}
		sa = (*AI.found[0])->name;
		sb = (*AI.found[1])->name;
		if (kw)
			squawk("Both '%s' and '%s' match '%s'",
				sa, sb, s);
		else
			squawk("'%s' is ambiguous:\n\tit matches both '%s' and '%s'",
				s, sa, sb);
		return AI.found[0];
		}
	kws = (KWseen*)KWgetmem(sizeof(KWseen));
	memset(kws, 0, sizeof(KWseen));
	if ((kws->kw = kw))
		s = kw->name;
	else {
		curkws = kws;
		kws->mnext = &KWmissing;
		KWmissing.mprev = (kws->mprev = KWmissing.mprev)->mnext = kws;
		s = nidr_KWscopy(s);
		}
	kws->name = s;
	if (KW_p >= KW_pe) {
		KW_p = (KWseen**)KWgetmem(32*sizeof(KWseen*));
		KW_pe = KW_p + 32;
		}
	*(kwsp = KW_p++) = kws;
	AVL_insert((const Element*)kwsp, AVLT);
	return kwsp;
	}

 static void
mixed_squawk(void)
{
	squawk("values for %s cannot be both strings and numbers",
		KW_cur->name);
	}

 static void
nidr_bufr_relaxed(Real r, int state)
{
	int n;

	if (KWval.rstate && !state)
		finish_rexpand();
	if (!(n = KWval.n)) {
		KWval.r = KWvalbuf;
		KWvalmax.r = KWvalbuf + (KWvalmax.n = nKWvalbuf);
		}
	else if (KWval.s) {
		mixed_squawk();
		return;
		}
	if (n >= KWvalmax.n)
		KWvalbuf_inc();
	KWval.r[KWval.n++] = r;
	if (state | KWval.rstate)
		rexpand(state);
	}

 static void
nidr_bufs_relaxed(const char *s)
{
	int n;

	if (!(n = KWval.n)) {
		KWval.s = (const char**)KWvalbuf;
		KWvalmax.n = (nKWvalbuf*sizeof(Real))/sizeof(char*);
		KWvalmax.s = KWval.s + KWvalmax.n;
		}
	else if (KWval.r) {
		mixed_squawk();
		return;
		}
	if (n >= KWvalmax.n)
		KWvalbuf_inc();
	KWval.s[KWval.n++] = s;
	}

 static void kw_setup2(KWseen*);

 static void
kw_finish1(KWseen *kws)
{
	KeyWord *kw;
	int n;
	size_t L;

	if (KWval.rstate)
		finish_rexpand();
	kws->nvals = n = KWval.n;
	KWval.n = 0;
	if (KWval.r) {
		L = n*sizeof(Real);
		memcpy(kws->rvals = (Real*)KWgetmem(L), KWval.r, L);
		KWval.r = 0;
		}
	else if (KWval.s) {
		L = n*sizeof(char*);
		memcpy(kws->svals = (char**)KWgetmem(L), KWval.s, L);
		if ((kw = kws->kw) && kw->kind & KWKind_Libname) {
			read_lib(KWval.s[0], kw);
			if (kw->kw)
				kw_setup2(kws);
			}
		KWval.s = 0;
		}
	}

 static void*
Alloc1(size_t len)
{
	void *rv = malloc(len);
	if (!rv) {
		fprintf(stderr, "malloc(%lu) failure in Alloc1\n", (unsigned long)len);
		exit(1);
		}
	return rv;
	}

 static void
AVL_Clear(void)
{
	while(ToClear > ToClear0)
		(*--ToClear)->kind &= ~KWKind_Hashed;
	AVL_Tree_free(&AVLT);
	if (AVLKWP)
		AVL_Tree_free(&AVLKWP);
	}

 static void
kw_setup1(KeyWord *kw)
{
	KWseen *kws, *kws1;
	KeyWord *kw1;

	if ((kw1 = kw->kw))
		kwnext_setup(kw1, kw->nkw);
	if (!KWmblk0)
		KWmeminit();
	if (AVLT)
		AVL_Clear();
	AVLT = AVL_Tree_alloc(0, (AVL_Elcomp)avlcmp, Alloc1);
	KW_cur = KWs0 = kws = (KWseen*)KWgetmem(sizeof(KWseen));
	memset(kws, 0, sizeof(KWseen));
	kws->name = kw->name;
	kws->kw = kw;
	kws->lcn = &kws->mprev;
	if (kw1) {
		while(!kw1->name)
			kw1 = kw1->kwnext;
		for(; kw1; kw1 = kw1->kwnext) {
			kws1 = *KWhash(kw1->name, kw1);
			kws1->parent = kws;
			}
		}
	}

 static KWseen**
kw_setup3(KWseen **kwtodo1, KWseen *kws, KeyWord *kw)
{
	KWseen *kws1, **kwsp;

	for(; kw; kw = kw->kwnext) {
		kwsp = KWhash(kw->name, kw);
		kws1 = *kwsp;
		if (kws1->comment) {
			kw->comment = kws1->comment;
			kws1->comment = 0;
			}
		if (kws1->parent) {
			kws1 = (KWseen*)KWgetmem(sizeof(KWseen));
			memset(kws1, 0, sizeof(KWseen));
			kws1->kw = kw;
			kws1->name = kw->name;
			*kwsp = kws1;
			}
		kws1->parent = kws;
		if (!kws1->kw) {
			kws1->mprev->mnext = kws1->mnext;
			kws1->mnext->mprev = kws1->mprev;
			*kwtodo1 = kws1;
			kwtodo1 = kws1->lcn = &kws1->mprev;
			*kws->lcn = kws1;
			kws->lcn = &kws1->mnext;
			}
		kws1->kw = kw;
		}
	return kwtodo1;
	}

 static void
bumpToClear(void)
{
	KeyWord **ntc;
	size_t L, L1;

	L = ToClearEnd - ToClear0;
	L1 = L << 1;
	ntc = (KeyWord**)Alloc("bumpToClear", L1*sizeof(KeyWord*));
	memcpy(ntc, ToClear0, L*sizeof(KeyWord*));
	free(ToClear0);
	ToClear0 = ntc;
	ToClear  = ntc + L;
	ToClearEnd = ntc + L1;
	}

 static void
kw_setup2(KWseen *kws)
{
	KWpair kwp, *pkwp;
	KWseen *kws1, *kws2, *kws3, *kwtodo, **kwtodo1, **pkws;
	KeyWord *kw, *kw1;

	kwtodo1 = &kwtodo;
	for(;;) {
		kw = kws->kw;
		if ((kw1 = kw->kw)) {
			kwnext_setup(kw1, kw->nkw);
			kws2 = kws;
			while(!kw1->name) {
				if (!AVLKWP)
					AVLKWP = AVL_Tree_alloc(0, (AVL_Elcomp)kwpcmp, Alloc1);
				if (kw1->kind & KWKind_Hashed) {
					kwp.kw = kw1->kw;
					kwp.kws = 0;
					pkwp = (KWpair*)AVL_find((const Element*)&kwp, AVLKWP);
					kws2 = pkwp->kws;
					}
				else {
					if (ToClear >= ToClearEnd)
						bumpToClear();
					*ToClear++ = kw1;
					kw1->kind |= KWKind_Hashed;
					pkwp = (KWpair*)KWgetmem(sizeof(KWpair) + sizeof(KWseen));
					kws1 = (KWseen*)(pkwp + 1);
					pkwp->kw = kw1->kw;
					pkwp->kws = kws1;
					memset(kws1, 0, sizeof(KWseen));
					kws1->kw = kw1;
					kws1->name = kws->name;
					kws1->lcn = &kws1->mprev;
					kws1->parent = kws2;
					*kws2->lcn = 0;
					for(pkws = &kws2->mprev;
						(kws3 = *pkws) && !kws3->name;
						pkws = &kws3->mnext);
					kws1->mnext = *pkws;
					if (pkws == kws2->lcn)
						kws2->lcn = &kws1->mnext;
					kws2 = *pkws = kws1;
					kwnext_setup(kw1->kw, kw1->nkw);
					kwtodo1 = kw_setup3(kwtodo1, kws1, kw1->kw);
					AVL_insert((const Element*)pkwp, AVLKWP);
					}
				kw1 = kw1->kwnext;
				}
			if (kw->nkw)
				kwtodo1 = kw_setup3(kwtodo1, kws2, kw1);
			}
		*kwtodo1 = 0;
		if (!kwtodo)
			break;
		kws = kwtodo;
		kw = kws->kw;
		if (!(kwtodo = kwtodo->mprev))
			kwtodo1 = &kwtodo;
		}
	}

 static KeyWord *
nidr_identifier_relaxed(const char *name)
{
	KWseen *kws, *kws1;
	KeyWord *kw;
	int tryagain;

	kw_finish1(KW_cur);
 top:
	KW_cur = kws = *KWhash(name, 0);
	if ((kw = kws->kw)) {
		curid = kw;
		if (kws->lcn)
			squawk("'%s' already seen", kw->name);
		else {
			if (kws->comment) {
				kw->comment = kws->comment;
				kws->comment = 0;
				}
			kws->lcn = &kws->mprev;
			kws1 = kws->parent;
			*kws1->lcn = kws;
			kws1->lcn = &kws->mnext;
			if (kw->kw)
				kw_setup2(kws);
			if (kw->kind & KWKind_Dynlib) {
				kw = kw_insert(kw, &tryagain);
				if (kw->kw) {
					kws->kw = kw;
					kw_setup2(kws);
					}
				if (tryagain)
					goto top;
				}
			}
		}
	return (KeyWord*)kws;	/* just needs to be nonzero; won't be dereferenced */
	}

 static void
num_expected(KeyWord *kw, int n)
{
	squawk("expected numerical value%s for %s, not quoted strings",
		"s" + (n == 1), kw->name);
	}

 static void
kw_process(KWseen *kws)
{
	KWseen *kws1;
	KeyWord *kw;
	Real *r;
	Uint k;
	int i, n;

	kw = kws->kw;
	if (kw->name) {
		if (kws != KWs0 && !nidr_identifier_strict(kw->name))
			botch("nidr_identifier_strict did not find \"%s\"", kw->name);
		if ((n = KWval.n = kws->nvals)) {
			KWval.i = 0;
			KWval.r = 0;
			KWval.s = 0;
			KWval.rstate = 0;
			switch(k = kw->kind & KWKind_Mask) {
			  case 0:
				squawk("No values may be specified for %s", kw->name);
				break;

			  case KWKind_Int:
				if (!(r = kws->rvals)) {
					num_expected(kw,n);
					break;
					}
				KWval.i = (int*)KWvalbuf;
				for(i = 0; i < n; i++)
					KWval.i[i] = (int)r[i];
				break;

			  case KWKind_Real:
				if (!(r = kws->rvals)) {
					num_expected(kw,n);
					break;
					}
				KWval.r = r;
				break;

			  case KWKind_Str:
				if (!(KWval.s = (const char **)kws->svals))
					squawk("expected string value%s for %s",
						"s" + (n == 1), kw->name);
			  }
			}
		}
	*kws->lcn = 0;
	for(kws1 = kws->mprev; kws1; kws1 = kws1->mnext)
		kw_process(kws1);
	}

 static void
kw_finish2(void)
{
	KWseen *kws, *kwe;

	kw_finish1(KW_cur);
	kwe = &KWmissing;
	for(kws = KWmissing.mnext; kws != kwe; kws = kws->mnext) {
		squawk("unrecognized identifier '%s'", kws->name);
		}
	KWmissing.mnext = KWmissing.mprev = &KWmissing;
	kw_process(KWs0);
	KWs0 = 0;
	AVL_Clear();
	}

 static void
kw_finish3(void)
{
	KWmblk1 = KWmblk0;
	KWmem0 = (char*)(KWmblk0 + 1);
	KWmem1 = KWmem0 + KWmblk_gulp;
	KW_p = KW_pe = 0;
	}

 void (*nidr_bufr)(Real,int) = nidr_bufr_relaxed;
 void (*nidr_bufs)(const char*) = nidr_bufs_relaxed;
 KeyWord *(*nidr_identifier)(const char*) = nidr_identifier_relaxed;

 void
nidr_set_strict(int n)
{
	if ((strict = n)) {
		nidr_bufr = nidr_bufr_strict;
		nidr_bufs = nidr_bufs_strict;
		nidr_identifier = nidr_identifier_strict;
		}
	else {
		nidr_bufr = nidr_bufr_relaxed;
		nidr_bufs = nidr_bufs_relaxed;
		nidr_identifier = nidr_identifier_relaxed;
		}
	}

 int
nidr_cleanup(void)
{
	KWmblk *mb, *mb1;
	Sbuf *sb, *sb1;

	if (curkw)
		nidr_keyword_finish();
	if (dumpfile) {
		if (OutsideComment)
			dumpcomment(&OutsideComment);
		if (dumpfile != stdout) {
			fclose(dumpfile);
			dumpfile = 0;
			}
		if (nidr_comment)
			comment_reset();
		}
	if (ToClear0) {
		free(ToClear0);
		ToClear = ToClear0 = 0;
		}
	if ((mb1 = KWmblk0)) {
		KWmblk0 = 0;
		do {
			mb = mb1;
			mb1 = mb->next;
			free(mb);
			} while(mb1);
		}
	if (KWvalbuf) {
		free(KWvalbuf);
		KWvalbuf = 0;
		}
	if ((sb1 = KWsbuf0.next)) {
		KWsbuf0.next = 0;
		do {
			sb = sb1;
			sb1 = sb->next;
			free(sb);
			} while(sb1);
		}
	if (AVLT)
		AVL_Clear();
	return nidr_parse_error();
	}

 void
nidr_setup(const char *parser, FILE *df)
{
	const char *s;
	int comkeep, oneline;

	if (!(s = parser))
		return;
	if (!strncmp(s,"nidr",4))
		s += 4;
	if (!strncmp(parser,"strict",6)) {
		nidr_set_strict(1);
		s += 6;
		}
	comkeep = oneline = 0;
	if (*s == '-') for(;;) {
		switch(*++s) {
		  case '1':
			++oneline;
			continue;
		  case 'p':
			++primary;
			continue;
		  case 'c':
			++comkeep;
			continue;
		  }
		break;
		}
	if (df)
		dumpfile = df;
	else if (s[0] == ':' && s[1]) {
		if (s[1] == '-' && !s[2])
			dumpfile = df = stdout;
		else {
			dumpfile = df = fopen(++s,"w");
			if (!dumpfile) {
				fprintf(stderr, "Cannot open \"%s\"\n", s);
				exit(1);
				}
			}
		}
	if (df) {
		if (oneline)
			dumpvals = dumpvals1;
		if (comkeep)
			comment_setup();
		}
	}
