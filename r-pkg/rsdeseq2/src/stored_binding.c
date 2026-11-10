/* Read an already stored/forced binding. Never evaluate an active or delayed
 * binding to discover closure data. No numerical calculation happens here. */
#include <R.h>
#include <Rinternals.h>
#include <Rversion.h>
SEXP rsdeseq2_stored_binding(SEXP name, SEXP env, SEXP context) {
    if (TYPEOF(name) != STRSXP || XLENGTH(name) != 1 || STRING_ELT(name, 0) == NA_STRING ||
        TYPEOF(context) != STRSXP || XLENGTH(context) != 1 || STRING_ELT(context, 0) == NA_STRING ||
        TYPEOF(env) != ENVSXP) {
        Rf_error("stored binding inspection requires a name, environment and context");
    }
#if R_VERSION >= R_Version(4, 6, 0)
    SEXP symbol = Rf_installChar(STRING_ELT(name, 0));
    R_BindingType_t type = R_GetBindingType(symbol, env);
    if (type == R_BindingTypeActive) {
        Rf_error("active %s bindings are unsupported", Rf_translateCharUTF8(STRING_ELT(context, 0)));
    }
    if (type == R_BindingTypeDelayed || type == R_BindingTypeMissing) {
        Rf_error("%s constants must be stored primitive values; delayed expressions are unsupported", Rf_translateCharUTF8(STRING_ELT(context, 0)));
    }
    if (type != R_BindingTypeValue && type != R_BindingTypeForced) {
        Rf_error("unknown stored binding '%s'", Rf_translateCharUTF8(STRING_ELT(name, 0)));
    }
    /* Both accepted types have an existing value. R_getVar cannot need to
     * evaluate a promise or invoke an active-binding function in these cases. */
    return R_getVar(symbol, env, FALSE);
#else
    Rf_error("reading forced closure bindings requires R 4.6 or newer");
#endif
    return R_NilValue;
}
