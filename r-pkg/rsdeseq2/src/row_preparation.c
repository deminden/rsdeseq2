#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <stddef.h>
#include "validation.h"
extern const char *rsdeseq2_r_last_error(void);
extern int rsdeseq2_r_validate_counts(size_t,size_t,const double *);
extern int rsdeseq2_r_normal_shrinkage_rows(size_t,size_t,const double *,const double *,const double *,const int32_t *,int32_t *,size_t *);
extern int rsdeseq2_r_offset_factors(size_t,size_t,const double *,const double *,const double *,double *);
extern int rsdeseq2_r_frozen_local_grid(size_t,size_t,const double *,const double *,const double *,double *);
static void checked(int status) {if(status)Rf_error("rsdeseq2 Rust backend: %s",rsdeseq2_r_last_error());}
SEXP rsdeseq2_frozen_local_grid(SEXP counts,SEXP sf,SEXP nf) {
    validate_matrix(counts,"counts");
    validate_normalization(counts,sf,nf,R_NilValue);
    SEXP grid=PROTECT(Rf_allocVector(REALSXP,999));
    checked(rsdeseq2_r_frozen_local_grid((size_t)matrix_rows(counts),(size_t)matrix_cols(counts),REAL(counts),
        sf==R_NilValue?NULL:REAL(sf),nf==R_NilValue?NULL:REAL(nf),REAL(grid)));
    UNPROTECT(1);return grid;
}
SEXP rsdeseq2_validate_counts(SEXP counts) {
    validate_matrix(counts,"counts");
    checked(rsdeseq2_r_validate_counts((size_t)matrix_rows(counts),(size_t)matrix_cols(counts),REAL(counts)));
    return counts;
}
SEXP rsdeseq2_normal_shrinkage_rows(SEXP counts,SEXP dispersions,SEXP fitted,SEXP excluded) {
    validate_matrix(counts,"counts");R_xlen_t genes=matrix_rows(counts);
    validate_vector(dispersions,REALSXP,genes,"dispersions");
    validate_vector(fitted,REALSXP,genes,"fitted dispersions");
    validate_optional_vector(excluded,LGLSXP,genes,"excluded rows");
    SEXP selected=PROTECT(Rf_allocVector(INTSXP,genes));size_t used=0;
    checked(rsdeseq2_r_normal_shrinkage_rows((size_t)genes,(size_t)matrix_cols(counts),REAL(counts),REAL(dispersions),REAL(fitted),
        excluded==R_NilValue?NULL:(const int32_t *)LOGICAL(excluded),(int32_t *)INTEGER(selected),&used));
    if(used>(size_t)genes)Rf_error("invalid selected row count");
    SEXP result=Rf_xlengthgets(selected,(R_xlen_t)used);UNPROTECT(1);return result;
}
SEXP rsdeseq2_offset_factors(SEXP counts,SEXP sf,SEXP nf,SEXP offsets) {
    validate_matrix(counts,"counts");R_xlen_t genes=matrix_rows(counts),samples=matrix_cols(counts);
    validate_optional_vector(sf,REALSXP,samples,"size factors");
    if(nf!=R_NilValue){validate_matrix(nf,"normalization factors");if(matrix_rows(nf)!=genes||matrix_cols(nf)!=samples)Rf_error("normalization factors dimensions must match counts");}
    validate_vector(offsets,REALSXP,samples,"formula offsets");
    SEXP result=PROTECT(Rf_allocMatrix(REALSXP,(int)genes,(int)samples));
    checked(rsdeseq2_r_offset_factors((size_t)genes,(size_t)samples,sf==R_NilValue?NULL:REAL(sf),nf==R_NilValue?NULL:REAL(nf),REAL(offsets),REAL(result)));
    SEXP names=Rf_getAttrib(counts,R_DimNamesSymbol);if(names!=R_NilValue)Rf_setAttrib(result,R_DimNamesSymbol,names);
    UNPROTECT(1);return result;
}
