AC_INIT(jeweltreewriter, [1.00])
AC_CONFIG_SRCDIR([configure.ac])
# keep libtool/automake helper scripts here (otherwise libtoolize searches parent dirs)
AC_CONFIG_AUX_DIR([.])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow"

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
