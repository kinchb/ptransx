E_SYSTEM} in
  Darwin*)
      { $as_echo "$as_me:${as_lineno-$LINENO}: checking for main in -lcc_dynamic" >&5
$as_echo_n "checking for main in -lcc_dynamic... " >&6; }
if ${ac_cv_lib_cc_dynamic_main+:} false; then :
  $as_echo_n "(cached) " >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lcc_dynamic  $LIBS"
cat confdefs.h - <<_ACEOF >conftest.$ac_ext
/* end confdefs.h.  */


int
main ()
{
return main ();
  ;
  return 0;
}
_ACEOF
if ac_fn_c_try_link "$LINENO"; then :
  ac_cv_lib_cc_dynamic_main=yes
else
  ac_cv_lib_cc_dynamic_main=no
fi
rm -f core conftest.err conftest.$ac_objext \
    conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
{ $as_echo "$as_me:${as_lineno-$LINENO}: result: $ac_cv_lib_cc_dynamic_main" >&5
$as_echo "$ac_cv_lib_cc_dynamic_main" >&6; }
if test "x$ac_cv_lib_cc_dynamic_main" = xyes; then :
  LIBS="$LIBS -lcc_dynamic"
fi

    ;;
  OSF*)
      { $as_echo "$as_me:${as_lineno-$LINENO}: checking for vsnprintf in -ldb" >&5
$as_echo_n "checking for vsnprintf in -ldb... " >&6; }
if ${ac_cv_lib_db_vsnprintf+:} false; then :
  $as_echo_n "(cached) " >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-ldb  $LIBS"
cat confdefs.h - <<_ACEOF >conftest.$ac_ext
/* end confdefs.h.  */

/* Override any GCC internal prototype to avoid an error.
   Use char because int might match the return type of a GCC
   builtin and then its argument prototype would still apply.  */
#ifdef __cplusplus
extern "C"
#endif
char vsnprintf ();
int
main ()
{
return vsnprintf ();
  ;
  return 0;
}
_ACEOF
if ac_fn_c_try_link "$LINENO"; then :
  ac_cv_lib_db_vsnprintf=yes
else
  ac_cv_lib_db_vsnprintf=no
fi
rm -f core conftest.err conftest.$ac_objext \
    conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
{ $as_echo "$as_me:${as_lineno-$LINENO}: result: $ac_cv_lib_db_vsnprintf" >&5
$as_echo "$ac_cv_lib_db_vsnprintf" >&6; }
if test "x$ac_cv_lib_db_vsnprintf" = xyes; then :
  LIBS="$LIBS -ldb"
fi

    ;;
  SunOS*)
      { $as_echo "$as_me:${as_lineno-$LINENO}: checking for main in -lmvec" >&5
$as_echo_n "checking for main in -lmvec... " >&6; }
if ${ac_cv_lib_mvec_main+:} false; then :
  $as_echo_n "(cached) " >&6
else
  ac_check_lib_save_LIBS=$LIBS
LIBS="-lmvec  $LIBS"
cat confdefs.h - <<_ACEOF >conftest.$ac_ext
/* end confdefs.h.  */


int
main ()
{
return main ();
  ;
  return 0;
}
_ACEOF
if ac_fn_c_try_link "$LINENO"; then :
  ac_cv_lib_mvec_main=yes
else
  ac_cv_lib_mvec_main=no
fi
rm -f core conftest.err conftest.$ac_objext \
    conftest$ac_exeext conftest.$ac_ext
LIBS=$ac_check_lib_save_LIBS
fi
{ $as_echo "$as_me:${as_lineno-$LINENO}: result: $ac_cv_lib_mvec_main" >&5
$as_echo "$ac_cv_lib_mvec_main" >&6; }
if test "x$ac_cv_lib_mvec_main" = xyes; then :
  LIBS="$LIBS -lmvec"
fi

      { $as_echo "$as_me:${as_lineno-$LINENO}: checking for main in -lsunmath" >&5
$as_echo_n "checking for main in -lsunmath... " >&6; }
if ${ac_cv_lib_sunmath_main+:} false; then :
  $as_echo_n "(cached) " >&6
else
  ac_check