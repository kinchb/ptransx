dnl $LynxId: configure.in,v 1.288 2017/04/30 16:02:52 tom Exp $
dnl
dnl Process this file with autoconf to produce a configure script.
dnl
dnl created jan/1997
dnl by T.E.Dickey <dickey@invisible-island.net>
dnl and Jim Spath <jspath@mail.bcpl.lib.md.us>
dnl
dnl ---------------------------------------------------------------------------
dnl Copyright 1997-2016,2017 by Thomas E. Dickey
dnl
dnl Permission to use, copy, modify, and distribute this software and its
dnl documentation for any purpose and without fee is hereby granted,
dnl provided that the above copyright notice appear in all copies and that
dnl both that copyright notice and this permission notice appear in
dnl supporting documentation, and that the name of the above listed
dnl copyright holder(s) not be used in advertising or publicity pertaining
dnl to distribution of the software without specific, written prior
dnl permission.
dnl
dnl THE ABOVE LISTED COPYRIGHT HOLDER(S) DISCLAIM ALL WARRANTIES WITH REGARD
dnl TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
dnl AND FITNESS, IN NO EVENT SHALL THE ABOVE LISTED COPYRIGHT HOLDER(S) BE
dnl LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
dnl WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
dnl ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
dnl OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
dnl ---------------------------------------------------------------------------
dnl
dnl ask PRCS to plug-in the project-version for the configure-script.
dnl $Format: "AC_REVISION($ProjectVersion$)"$
AC_REVISION(2.8.9dev.14)

# Save the original $CFLAGS so we can distinguish whether the user set those
# in the environment, or whether autoconf added -O and -g options:
ORIGINAL_CFLAGS="$CFLAGS"

# For autoconf 2.13, make sure we have no cache file at the beginning of this
# script.  That fixes problems with tests whose cached values change from one
# run to the next, as well as with tests that are order-dependent.
rm -f config.cache

AC_PREREQ(2.13.20020210)
AC_INIT(userdefs.h)

# autoconf 2.5x defaults to no cache file; we need the cache file's information
# for building the config page.  But start with it empty to avoid confusion by
# people who don't do a "make distclean" after applying patches.
cache_file=config.cache
rm -f config.cache; touch config.cache

CONFIG_H=lynx_cfg.h
AC_CONFIG_HEADER($CONFIG_H:config.hin)
AC_SUBST(CONFIG_H)

CF_CHECK_CACHE([AC_CANONICAL_SYSTEM])
AC_ARG_WITH(system-type,
[  --with-system-type=XXX  test: override derived host system-type],
[AC_MSG_WARN(overriding system type $host_os to $withval)
 host_os=$withval])

AC_ARG_PROGRAM

PACKAGE=lynx
dnl ask PRCS to plug-in the project-version for the packages.
# $Format: "VERSION=$ProjectVersion$"$
VERSION=2.8.9dev.14

AC_SUBST(PACKAGE)
AC_SUBST(VERSION)

AC_MSG_CHECKING(for DESTDIR)
CF_WITH_PATH(destdir,
[  --with-destdir=XXX      set DESTDIR destination for install],
DESTDIR,
[$DESTDIR],
[$DESTDIR])
AC_MSG_RESULT($DESTDIR)

AC_PREFIX_DEFAULT(/usr/local)

dnl --------------------------------------------------------------------------
dnl Checks for location of programs
dnl --------------------------------------------------------------------------

dnl Only add to this case statement