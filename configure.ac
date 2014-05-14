AC_PREREQ(2.61)
AC_INIT(KdV, 1.0, amd262@wildcats)
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AM_INIT_AUTOMAKE([foreign])

AC_PROG_CC_C99
if test -z "$CC" || test "$CC" = "no"; then
  AC_MSG_ERROR([No suitable C99 compiler found])
fi
AM_PROG_CC_C_O

AC_OUTPUT(Makefile src/Makefile)
