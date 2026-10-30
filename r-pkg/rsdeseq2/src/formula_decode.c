#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "formula_transport.h"
typedef struct { const uint8_t *data; size_t at,size; } Reader;
static const uint8_t *wire_take(Reader *r,size_t n){if(n>r->size-r->at)Rf_error("truncated formula output");const uint8_t *p=r->data+r->at;r->at+=n;return p;}
static uint8_t wire_byte(Reader *r){return wire_take(r,1)[0];}
static int wire_flag(Reader *r){int v=wire_byte(r);if(v>1)Rf_error("invalid formula output flag");return v;}
static uint64_t wire_length(Reader *r){const uint8_t *b=wire_take(r,8);uint64_t n=0;for(int i=0;i<8;++i)n|=(uint64_t)b[i]<<(8*i);return n;}
static R_xlen_t wire_count(Reader *r,size_t unit){uint64_t n=wire_length(r);if(n>R_XLEN_T_MAX||n>(r->size-r->at)/unit)Rf_error("invalid formula output length");return (R_xlen_t)n;}
static SEXP wire_text(Reader *r){R_xlen_t n=wire_count(r,1);if(n>INT_MAX)Rf_error("formula text exceeds string length");return Rf_mkCharLenCE((const char *)wire_take(r,(size_t)n),(int)n,CE_UTF8);}
static SEXP wire_names(Reader *r){if(!wire_flag(r))return R_NilValue;R_xlen_t n=wire_count(r,8);SEXP result=PROTECT(Rf_allocVector(STRSXP,n));for(R_xlen_t i=0;i<n;++i)SET_STRING_ELT(result,i,wire_text(r));UNPROTECT(1);return result;}
static SEXP wire_value(Reader *r,int depth){
    if(depth>RSDESEQ2_FORMULA_MAX_DEPTH)Rf_error("formula output nesting exceeds 1024");
    int tag=wire_byte(r);R_xlen_t n=wire_count(r,tag==1||tag==2?8:1);
    SEXPTYPE type=tag==0?NILSXP:tag==1?REALSXP:tag==2?INTSXP:tag==3?LGLSXP:tag==4?STRSXP:tag==5?VECSXP:ANYSXP;
    if(type==ANYSXP||(tag==0&&n))Rf_error("invalid formula output type");
    SEXP out=PROTECT(tag==0?R_NilValue:Rf_allocVector(type,n));int protected=1;
#ifdef RSDESEQ2_TEST_FORMULA_UNWIND
    rsdeseq2_formula_inject(2);
#endif
    for(R_xlen_t i=0;i<n;++i){
        if(tag==1||tag==2){uint64_t bits=wire_length(r);double x;memcpy(&x,&bits,8);if(tag==1)REAL(out)[i]=x;else{if(ISNAN(x))INTEGER(out)[i]=NA_INTEGER;else if(!R_FINITE(x)||x<=INT_MIN||x>INT_MAX||x!=(int)x)Rf_error("invalid formula integer output");else INTEGER(out)[i]=(int)x;}}
        else if(tag==3){int v=wire_byte(r);if(v>2)Rf_error("invalid formula logical output");LOGICAL(out)[i]=v==2?NA_LOGICAL:v;}
        else if(tag==4)SET_STRING_ELT(out,i,wire_flag(r)?wire_text(r):NA_STRING);
        else if(tag==5){SEXP child=PROTECT(wire_value(r,depth+1));SET_VECTOR_ELT(out,i,child);UNPROTECT(1);}
    }
#ifdef RSDESEQ2_TEST_FORMULA_UNWIND
    rsdeseq2_formula_inject(3);
#endif
    uint64_t rows=0,cols=0;int matrix=wire_flag(r);
    if(matrix){
        rows=wire_length(r);cols=wire_length(r);
        if(rows>INT_MAX||cols>INT_MAX||rows*cols!=(uint64_t)n)Rf_error("invalid formula matrix output");
        SEXP reordered=PROTECT(Rf_allocVector(type,n));++protected;
        for(R_xlen_t i=0;i<n;++i){R_xlen_t at=i/(R_xlen_t)cols+(i%(R_xlen_t)cols)*(R_xlen_t)rows;
            if(tag==1)REAL(reordered)[at]=REAL(out)[i];else if(tag==2)INTEGER(reordered)[at]=INTEGER(out)[i];else if(tag==3)LOGICAL(reordered)[at]=LOGICAL(out)[i];else if(tag==4)SET_STRING_ELT(reordered,at,STRING_ELT(out,i));else if(tag==5)SET_VECTOR_ELT(reordered,at,VECTOR_ELT(out,i));else Rf_error("invalid NULL formula matrix");
        }
        out=reordered;SEXP dims=PROTECT(Rf_allocVector(INTSXP,2));INTEGER(dims)[0]=(int)rows;INTEGER(dims)[1]=(int)cols;Rf_setAttrib(out,R_DimSymbol,dims);UNPROTECT(1);
    }
    SEXP nm=PROTECT(wire_names(r));
    if(nm!=R_NilValue){
        if(matrix){
            if(XLENGTH(nm)!=n)Rf_error("formula element names do not match values");
            SEXP reordered=PROTECT(Rf_allocVector(STRSXP,n));
            for(R_xlen_t i=0;i<n;++i)SET_STRING_ELT(reordered,i/(R_xlen_t)cols+(i%(R_xlen_t)cols)*(R_xlen_t)rows,STRING_ELT(nm,i));
            Rf_setAttrib(out,R_NamesSymbol,reordered);UNPROTECT(1);
        }else Rf_setAttrib(out,R_NamesSymbol,nm);
    }
    UNPROTECT(1);
    SEXP dimnames=PROTECT(Rf_allocVector(VECSXP,2));int any=0;
    for(int i=0;i<2;++i){SEXP nm=PROTECT(wire_names(r));SET_VECTOR_ELT(dimnames,i,nm);any|=nm!=R_NilValue;UNPROTECT(1);}
    for(int i=0;i<2;++i){SEXP nm=PROTECT(wire_names(r));if(nm!=R_NilValue){SEXP labels=VECTOR_ELT(dimnames,i);if(labels==R_NilValue)Rf_error("formula label names require labels");Rf_setAttrib(labels,R_NamesSymbol,nm);}UNPROTECT(1);}
    if(any)Rf_setAttrib(out,R_DimNamesSymbol,dimnames);UNPROTECT(1);
    R_xlen_t attrs=wire_count(r,9);
    for(R_xlen_t i=0;i<attrs;++i){SEXP name=PROTECT(wire_text(r));SEXP attr=PROTECT(wire_value(r,depth+1));Rf_setAttrib(out,Rf_installChar(name),attr);UNPROTECT(2);}
    UNPROTECT(protected);return out;
}
SEXP rsdeseq2_formula_read(const uint8_t *data,size_t size){
    if(size>RSDESEQ2_FORMULA_MAX_BYTES||size<4||memcmp(data,"RFX2",4))Rf_error("invalid formula output header");
    Reader reader={data,4,size};SEXP result=PROTECT(wire_value(&reader,0));if(reader.at!=size)Rf_error("trailing formula output bytes");UNPROTECT(1);return result;
}
