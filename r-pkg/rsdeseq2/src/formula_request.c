#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <stdint.h>
#include "validation.h"
#include "formula_transport.h"
typedef struct FormulaOutput FormulaOutput;
extern int rsdeseq2_r_formula_request_interruptible(const uint8_t *,size_t,FormulaOutput **,int (*)(void));
static void poll_interrupt(void *data){(void)data;R_CheckUserInterrupt();}
static int interrupted(void){return R_ToplevelExec(poll_interrupt,NULL)==FALSE;}
extern int rsdeseq2_r_formula_bytes(const FormulaOutput *,const uint8_t **,size_t *);
extern void rsdeseq2_r_formula_drop(FormulaOutput **);
extern const char *rsdeseq2_r_last_error(void);
static void checked(int status){if(status)Rf_error("rsdeseq2 Rust backend: %s",rsdeseq2_r_last_error());}
typedef struct { FormulaOutput *owner; } Output;
#ifdef RSDESEQ2_TEST_FORMULA_UNWIND
static int fail_stage=0;static size_t seen=0,dropped=0;
SEXP rsdeseq2_test_formula_transport(SEXP stage){fail_stage=validate_integer_scalar(stage,0,4,"formula test stage");SEXP result=PROTECT(Rf_allocVector(REALSXP,3));REAL(result)[0]=(double)seen;REAL(result)[1]=(double)dropped;REAL(result)[2]=(double)(seen-dropped);UNPROTECT(1);return result;}
void rsdeseq2_formula_inject(int stage){if(fail_stage==stage){fail_stage=0;Rf_error("injected formula transport error at stage %d",stage);}}
#endif
static SEXP decode(void *data){Output *output=(Output *)data;const uint8_t *bytes=NULL;size_t size=0;checked(rsdeseq2_r_formula_bytes(output->owner,&bytes,&size));
#ifdef RSDESEQ2_TEST_FORMULA_UNWIND
    rsdeseq2_formula_inject(1);
#endif
    SEXP result=PROTECT(rsdeseq2_formula_read(bytes,size));
#ifdef RSDESEQ2_TEST_FORMULA_UNWIND
    rsdeseq2_formula_inject(4);
#endif
    UNPROTECT(1);return result;
}
static void cleanup(void *data,Rboolean jump){(void)jump;Output *output=(Output *)data;rsdeseq2_r_formula_drop(&output->owner);
#ifdef RSDESEQ2_TEST_FORMULA_UNWIND
    ++dropped;
#endif
}
SEXP rsdeseq2_formula_request(SEXP request){
    /* Every R allocation/ALTREP access in encoding completes before ownership
     * leaves Rust. The continuation also exists before publishing any owner. */
    size_t size=rsdeseq2_formula_measure(request);uint8_t *bytes=(uint8_t *)R_alloc(size,1);rsdeseq2_formula_write(request,bytes,size);
    SEXP continuation=PROTECT(R_MakeUnwindCont());Output output={NULL};
    checked(rsdeseq2_r_formula_request_interruptible(bytes,size,&output.owner,interrupted));
#ifdef RSDESEQ2_TEST_FORMULA_UNWIND
    ++seen;
#endif
    SEXP result=R_UnwindProtect(decode,&output,cleanup,&output,continuation);UNPROTECT(1);return result;
}
