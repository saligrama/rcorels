#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([corels], [VERSION], [BUG-REPORT-ADDRESS])

## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
    # Use gsl-config for header and linker arguments
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
else
    AC_MSG_ERROR([gsl-config not found, is GSL1 installed?])
fi

AC_CONFIG_SRCDIR([src])
# AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])

AM_INIT_AUTOMAKE([-Wportability])

LT_INIT()


# Checks for programs.
AC_PROG_CC
AC_PROG_CXX

# Set language
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11

# Checks for libraries.
#AC_CHECK_LIB([c], [cxx], [printf], [], [
#echo "Standard C library is required for this program"
#exit -1])
AC_CHECK_LIB([m],[floor])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
LDFLAGS="$LDFLAGS -L/opt/local/lib -L/usr/local/lib -L/usr/lib64"
CPPFLAGS="$CPPFLAGS -I/opt/local/include -I/usr/local/include"
CXXFLAGS="$CXXFLAGS -I/opt/local/include -I/usr/local/include"
AC_CHECK_LIB([gmp], [__gmpz_clear], [], [
	echo "GMP (Gnu multiprecision) library is required for this program"
	exit -1])
AC_CHECK_LIB([gsl], [gsl_ran_poisson], [], [
	echo "GSL library is required for this program"
	exit -1])

AC_SUBST(GMP_CFLAGS, ${CFLAGS})
AC_SUBST(GMP_LIBS, ${LIBS})

AC_OUTPUT(src/Makevars)

# Checks for header files.
AC_CHECK_HEADERS([float.h gmp.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([gettimeofday memmove memset strdup strerror])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
