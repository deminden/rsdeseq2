#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include "validation.h"

typedef struct {
    const char *name;
    int32_t kind, coding;
    size_t width;
    const double *values;
    const int32_t *codes;
    const char *const *levels, *const *names;
    size_t contrast_columns;
    const char *const *contrast_names;
} ModelVariable;
extern const char *rsdeseq2_r_last_error(void);
extern int rsdeseq2_r_average_contrast(size_t,size_t,const int32_t *,size_t,const int32_t *,double *);
extern int rsdeseq2_r_model_matrix_shape(size_t,size_t,const ModelVariable *,size_t,const int32_t *,int32_t,size_t,size_t *,size_t *,int32_t *);
extern int rsdeseq2_r_model_matrix_fill(size_t,size_t,const ModelVariable *,size_t,const int32_t *,int32_t,size_t,size_t,double *,int32_t *,uint8_t *,size_t);
static void checked(int status) {if(status)Rf_error("rsdeseq2 Rust backend: %s",rsdeseq2_r_last_error());}
static const char *string_at(SEXP value,R_xlen_t index) {
    if(TYPEOF(value)!=STRSXP||index>=XLENGTH(value)||STRING_ELT(value,index)==NA_STRING)Rf_error("model names and labels must be nonmissing strings");
    return Rf_translateCharUTF8(STRING_ELT(value,index));
}
static const char *const *strings(SEXP value,R_xlen_t length) {
    validate_vector(value,STRSXP,length,"model labels");
    const char **result=(const char **)R_alloc((size_t)length,sizeof(*result));
    for(R_xlen_t i=0;i<length;++i)result[i]=string_at(value,i);
    return result;
}
static SEXP field(SEXP value,const char *name) {return named_list_field(value,name,"model variable");}

SEXP rsdeseq2_average_contrast(SEXP columns,SEXP positive,SEXP negative) {
    int count=validate_integer_scalar(columns,1,INT_MAX,"average contrast columns");
    validate_vector(positive,INTSXP,XLENGTH(positive),"average contrast positive indices");
    validate_vector(negative,INTSXP,XLENGTH(negative),"average contrast negative indices");
    SEXP result=PROTECT(Rf_allocVector(REALSXP,count));
    checked(rsdeseq2_r_average_contrast((size_t)count,(size_t)XLENGTH(positive),(const int32_t *)INTEGER(positive),(size_t)XLENGTH(negative),(const int32_t *)INTEGER(negative),REAL(result)));
    UNPROTECT(1);return result;
}

SEXP rsdeseq2_model_matrix(SEXP samples,SEXP variables,SEXP terms,SEXP intercept,SEXP response) {
    int rows=validate_integer_scalar(samples,1,INT_MAX,"model samples");
    validate_boolean(intercept,"model intercept");
    if(TYPEOF(variables)!=VECSXP)Rf_error("model variables must be a list");
    R_xlen_t count=XLENGTH(variables);
    if(count>INT_MAX)Rf_error("model variables exceed matrix dimensions");
    int response_index=validate_integer_scalar(response,0,(double)count,"model response");
    SEXP dims=Rf_getAttrib(terms,R_DimSymbol);
    if(TYPEOF(terms)!=INTSXP||TYPEOF(dims)!=INTSXP||XLENGTH(dims)!=2||INTEGER(dims)[0]!=count||INTEGER(dims)[1]<0)Rf_error("model term codes must be an integer variable-by-term matrix");
    validate_vector(terms,INTSXP,checked_product(count,INTEGER(dims)[1]),"model terms");
    ModelVariable *input=(ModelVariable *)R_alloc((size_t)count,sizeof(*input));
    for(R_xlen_t i=0;i<count;++i) {
        SEXP variable=VECTOR_ELT(variables,i);ModelVariable *v=&input[i];memset(v,0,sizeof(*v));
        SEXP name=field(variable,"name");validate_vector(name,STRSXP,1,"model variable name");v->name=string_at(name,0);
        v->kind=validate_integer_scalar(field(variable,"kind"),0,1,"model variable kind");
        SEXP value=field(variable,"values"),names=field(variable,"columnNames");
        if(v->kind==0) {
            validate_matrix(value,"numeric model variable");
            if(matrix_rows(value)!=rows)Rf_error("numeric model variable rows must match samples");
            v->width=(size_t)matrix_cols(value);v->values=REAL(value);
            if(names!=R_NilValue)v->names=strings(names,(R_xlen_t)v->width);
        } else {
            SEXP levels=field(variable,"levels");
            if(TYPEOF(levels)!=STRSXP||XLENGTH(levels)<2)Rf_error("model factors require at least two declared levels");
            v->width=(size_t)XLENGTH(levels);v->levels=strings(levels,XLENGTH(levels));
            validate_vector(value,INTSXP,rows,"model factor codes");v->codes=(const int32_t *)INTEGER(value);
            v->coding=validate_integer_scalar(field(variable,"coding"),0,5,"factor coding");
            if(v->coding==5) {
                SEXP contrast=field(variable,"contrasts");validate_matrix(contrast,"factor contrast matrix");
                if(matrix_rows(contrast)!=(R_xlen_t)v->width)Rf_error("factor contrast rows must match declared levels");
                v->values=REAL(contrast);v->contrast_columns=(size_t)matrix_cols(contrast);
                if(names!=R_NilValue)v->contrast_names=strings(names,(R_xlen_t)v->contrast_columns);
            }
        }
    }
    size_t columns=0,name_bytes=0;int32_t response_dropped=0;
    checked(rsdeseq2_r_model_matrix_shape((size_t)rows,(size_t)count,input,(size_t)INTEGER(dims)[1],(const int32_t *)INTEGER(terms),LOGICAL(intercept)[0],(size_t)response_index,&columns,&name_bytes,&response_dropped));
    if(response_dropped){Rf_warning("the response appeared on the right-hand side and was dropped");Rf_warning("response-only model term has no assigned columns");}
    if(columns>INT_MAX||name_bytes>R_XLEN_T_MAX)Rf_error("model matrix output dimensions overflow");
    checked_product(rows,(R_xlen_t)columns);
    SEXP result=PROTECT(Rf_allocMatrix(REALSXP,rows,(int)columns));
    SEXP assign=PROTECT(Rf_allocVector(INTSXP,(R_xlen_t)columns));
    uint8_t *buffer=(uint8_t *)R_alloc(name_bytes?name_bytes:1,1);
    checked(rsdeseq2_r_model_matrix_fill((size_t)rows,(size_t)count,input,(size_t)INTEGER(dims)[1],(const int32_t *)INTEGER(terms),LOGICAL(intercept)[0],(size_t)response_index,columns,REAL(result),(int32_t *)INTEGER(assign),buffer,name_bytes));
    SEXP names=PROTECT(Rf_allocVector(STRSXP,(R_xlen_t)columns));
    const char *next=(const char *)buffer;
    for(size_t i=0;i<columns;++i){SET_STRING_ELT(names,(R_xlen_t)i,Rf_mkCharCE(next,CE_UTF8));next+=strlen(next)+1;}
    SEXP dimnames=PROTECT(Rf_allocVector(VECSXP,2));SET_VECTOR_ELT(dimnames,1,names);
    Rf_setAttrib(result,R_DimNamesSymbol,dimnames);Rf_setAttrib(result,Rf_install("assign"),assign);
    UNPROTECT(4);return result;
}
