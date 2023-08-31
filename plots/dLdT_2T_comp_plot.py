/*
 * $LynxId: HTParse.c,v 1.78 2016/11/24 15:29:50 tom Exp $
 *
 *		Parse HyperText Document Address		HTParse.c
 *		================================
 */

#include <HTUtils.h>
#include <HTParse.h>

#include <LYUtils.h>
#include <LYLeaks.h>
#include <LYStrings.h>
#include <LYCharUtils.h>
#include <LYGlobalDefs.h>

#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#else
#ifdef __MINGW32__
#include <malloc.h>
#endif /* __MINGW32__ */
#endif

#ifdef USE_IDNA
#include <idna.h>
#include <idn-free.h>
#endif

#define HEX_ESCAPE '%'

struct struct_parts {
    char *access;
    char *host;
    char *absolute;
    char *relative;
    char *search;		/* treated normally as part of path */
    char *anchor;
};

#if 0				/* for debugging */
static void show_parts(const char *name, struct struct_parts *parts, int line)
{
    if (TRACE) {
	CTRACE((tfp, "struct_parts(%s) %s@%d\n", name, __FILE__, line));
	CTRACE((tfp, "   access   '%s'\n", NONNULL(parts->access)));
	CTRACE((tfp, "   host     '%s'\n", NONNULL(parts->host)));
	CTRACE((tfp, "   absolute '%s'\n", NONNULL(parts->absolute)));
	CTRACE((tfp, "   relative '%s'\n", NONNULL(parts->relative)));
	CTRACE((tfp, "   search   '%s'\n", NONNULL(parts->search)));
	CTRACE((tfp, "   anchor   '%s'\n", NONNULL(parts->anchor)));
    }
}
#define SHOW_PARTS(name) show_parts(#name, &name, __LINE__)
#else
#define SHOW_PARTS(name)	/* nothing */
#endif

/*	Strip white space off a string.				HTStrip()
 *	-------------------------------
 *
 * On exit,
 *	Return value points to first non-white character, or to 0 if none.
 *	All trailing white space is OVERWRITTEN with zero.
 */
char *HTStrip(char *s)
{
#define SPACE(c) ((c == ' ') || (c == '\t') || (c == '\n'))
    char *p;

    for (p = s; *p; p++) {	/* Find end of string */
	;
    }
    for (p--; p >= s; p--) {
	if (SPACE(*p))
	    *p = '\0';		/* Zap trailing blanks */
	else
	    break;
    }
    while (SPACE(*s))
	s++;			/* Strip leading blanks */
    return s;
}

/*	Scan a filename for its constituents.			scan()
 *	-------------------------------------
 *
 * On entry,
 *	name	points to a document name which may be incomplete.
 * On exit,
 *	absolute or relative may be nonzero (but not both).
 *	host, anchor and access may be nonzero if they were specified.
 *	Any which are nonzero point to zero terminated strings.
 */
static void scan(char *name,
		 struct struct_parts *parts)
{
    char *after_access;
    char *p;

    parts->access = NULL;
    parts->host = NULL;
    parts->absolute = NULL;
    parts->relative = NULL;
    parts->search = NULL;	/* normally not used - kw */
    parts->anchor = NULL;

    /*
     * Scan left-to-right for a scheme (access).
     */
    after_access = name;
    for (p = name; *p; p++) {
	if (*p == ':') {
	    *p = '\0';
	    parts->access = name;	/* Access name has been specified */
	    after_access = (p + 1);
	    break;
	}
	if (*p == '/' || *p == '#' || *p == ';' || *p == '?')
	    break;
    }

    /*
     * Scan left-to-right for a fragment (anchor).
     */
    for (p = after_access; *p; p++) {
	if (*p == '#') {
	    parts->anchor = (p + 1);
	    *p = '\0';		/* terminate the rest */
	    break;		/* leave things after first # alone - kw */
	}
    }

    /*
     * Scan left-to-right for a host or absolute path.
     */
    p = after_access;
    if (*p == '/') {
	if (p[1] == '/') {
	    parts->host = (p + 2);	/* host has been specified    */
	    *p = '\0';		/* Terminate access           */
	    p = StrChr(parts->host, '/');	/* look for end of host name if any */
	    if (p != NULL) {
		*p = '\0';	/* Terminate host */
		parts->absolute = (p + 1);	/* Root has been found */
	    } else {
		p = StrChr(parts->host, '?');
		if (p != NULL) {
		    *p = '\0';	/* Terminate host */
		    parts->search = (p + 1);
		}
	    }
	} else {
	    parts->absolute = (p + 1);	/* Root found but no host */
	}
    } else {
	parts->relative = (*after_access) ?
	    after_access : NULL;	/* NULL for "" */
    }

    /*
     * Check schemes that commonly have unescaped hashes.
     */
    if (parts->access && parts->anchor &&
    /* optimize */ StrChr("lnsdLNSD", *parts->access) != NULL) {
	if ((!parts->host && strcasecomp(parts->access, "lynxcgi")) ||
	    !strcasecomp(parts->access, "nntp") ||
	    !strcasecomp(parts->access, "snews") ||
	    !strcasecomp(parts->access, "news") ||
	    !strcasecomp(parts->access, "data")) {
	    /*
	     * Access specified but no host and not a lynxcgi URL, so the
	     * anchor may not really be one, e.g., news:j462#36487@foo.bar, or
	     * it's an nntp or snews URL, or news URL with a host.  Restore the
	     * '#' in the address.
	     */
	    /* but only if we have found a path component of which this will
	     * become part. - kw  */
	    if (parts->relative || parts->absolute) {
		*(parts->anchor - 1) = '#';
		parts->anchor = NULL;
	    }
	}
    }
}				/*scan */

#if defined(HAVE_ALLOCA) && !defined(LY_FIND_LEAKS)
#define LYalloca(x)        alloca(x)
#define LYalloca_free(x)   {}
#else
#define LYalloca(x)        malloc(x)
#define LYalloca_free(x)   free(x)
#endif

static char *strchr_or_end(char *string, int ch)
{
    char *result = StrChr(string, ch);

    if (result == 0) {
	result = string + strlen(string);
    }
    return result;
}

/*
 * Given a host specification that may end with a port number, e.g.,
 *	foobar:123
 * point to the ':' which begins the ":port" to make it simple to handle the
 * substring.
 *
 * If no port is found (or a syntax error), return null.
 */
char *HTParsePort(char *host, int *portp)
{
    int brackets = 0;
    char *result = NULL;

    *portp = 0;
    if (host != NULL) {
	while (*host != '\0' && result == 0) {
	    switch (*host++) {
	    case ':':
		if (brackets == 0 && isdigit(UCH(*host))) {
		    char *next = NULL;

		    *portp = (int) strtol(host, &next, 10);
		    if (next != 0 && next != host && *next == '\0') {
			result = (host - 1);
			CTRACE((tfp, "HTParsePort %d\n", *portp));
		    }
		}
		break;
	    case '[':		/* for ipv6 */
		++brackets;
		break;
	    case ']':		/* for ipv6 */
		--brackets;
		break;
	    }
	}
    }
    return result;
}

#ifdef USE_IDNA
static int hex_decode(int ch)
{
    int result = -1;

    if (ch >= '0' && ch <= '9')
	result = (ch - '0');
    else if (ch >= 'a' && ch <= 'f')
	result = (ch - 'a') + 10;
    else if (ch >= 'A' && ch <= 'F')
	result = (ch - 'A') + 10;
    return result;
}

/*
 * Convert in-place the given hostname to IDNA form.  That requires up to 64
 * characters, and we've