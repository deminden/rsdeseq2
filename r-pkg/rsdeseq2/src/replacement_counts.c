#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include "validation.h"
#include "fit_abi.h"

extern const char *rsdeseq2_r_last_error(void);
static void checked(int status) {
    if (status) Rf_error("rsdeseq2 Rust backend: %s", rsdeseq2_r_last_error());
}
typedef struct {
    Rsdeseq2ReplacementCounts **owner;
    size_t genes, samples;
} ReplacementCopy;

/* Only the disposable sanitizer build enables error injection/counters. The
 * release bridge has neither this registered hook nor these branches. */
#ifdef RSDESEQ2_TEST_REPLACEMENT_UNWIND
static int fail_stage = 0;
static size_t owners_seen = 0, owners_dropped = 0;
SEXP rsdeseq2_test_replacement_transport(SEXP stage) {
    fail_stage = validate_integer_scalar(stage, 0, 4, "replacement test stage");
    SEXP result = PROTECT(Rf_allocVector(REALSXP, 3));
    REAL(result)[0] = (double)owners_seen;
    REAL(result)[1] = (double)owners_dropped;
    REAL(result)[2] = (double)(owners_seen - owners_dropped);
    UNPROTECT(1); return result;
}
static void inject_error(int stage) {
    if (fail_stage == stage) {
        fail_stage = 0;
        Rf_error("injected replacement transport error at stage %d", stage);
    }
}
#endif

static SEXP copy_replacement(void *data) {
    ReplacementCopy *copy = (ReplacementCopy *)data;
    size_t genes = 0, samples = 0;
    checked(rsdeseq2_r_replacement_shape(*copy->owner, &genes, &samples));
    if (genes != copy->genes || samples != copy->samples || genes > INT_MAX || samples > INT_MAX)
        Rf_error("replacement output dimensions do not match fitted counts");
    checked_product((R_xlen_t)genes, (R_xlen_t)samples);
#ifdef RSDESEQ2_TEST_REPLACEMENT_UNWIND
    inject_error(1);
#endif
    SEXP result = PROTECT(Rf_allocMatrix(REALSXP, (int)genes, (int)samples));
#ifdef RSDESEQ2_TEST_REPLACEMENT_UNWIND
    inject_error(2);
    if (fail_stage == 3) {
        fail_stage = 0;
        /* Copy validates shape before constructing a writable output slice. */
        checked(rsdeseq2_r_replacement_copy(*copy->owner, genes + 1, samples, REAL(result)));
    }
#endif
    checked(rsdeseq2_r_replacement_copy(*copy->owner, genes, samples, REAL(result)));
#ifdef RSDESEQ2_TEST_REPLACEMENT_UNWIND
    inject_error(4);
#endif
    UNPROTECT(1); return result;
}
static void drop_replacement(void *data, Rboolean jump) {
    (void)jump;
    ReplacementCopy *copy = (ReplacementCopy *)data;
    rsdeseq2_r_replacement_drop(copy->owner);
#ifdef RSDESEQ2_TEST_REPLACEMENT_UNWIND
    if (*copy->owner == NULL) ++owners_dropped;
#endif
}

/* The continuation is protected and allocated before native fitting, while
 * there is no native owner to leak if R allocation fails. The Rust fitting
 * frame has already returned before any R operation below can raise/longjmp. */
SEXP rsdeseq2_finish_replacement(Rsdeseq2ReplacementCounts **owner,
                               size_t genes, size_t samples, SEXP continuation) {
    if (*owner == NULL) return R_NilValue;
    ReplacementCopy copy = {owner, genes, samples};
#ifdef RSDESEQ2_TEST_REPLACEMENT_UNWIND
    ++owners_seen;
#endif
    return R_UnwindProtect(copy_replacement, &copy, drop_replacement, &copy, continuation);
}
