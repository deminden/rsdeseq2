#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "formula_transport.h"
typedef struct { uint8_t *data; size_t at, capacity; } Writer;
static void wire_bytes(Writer *w,const void *v,size_t n) {
    if(n>w->capacity-w->at)Rf_error("formula transport exceeds size limit");
    if(w->data&&n)memcpy(w->data+w->at,v,n);w->at+=n;
}
static void wire_byte(Writer *w,uint8_t n){wire_bytes(w,&n,1);}
static void wire_length(Writer *w,uint64_t n){uint8_t b[8];for(int i=0;i<8;++i)b[i]=(uint8_t)(n>>(8*i));wire_bytes(w,b,8);}
static void wire_text(Writer *w,SEXP s) {
    if(TYPEOF(s)!=CHARSXP||s==NA_STRING)Rf_error("formula metadata requires nonmissing text");
    const char *v=Rf_translateCharUTF8(s);size_t n=strlen(v);wire_length(w,n);wire_bytes(w,v,n);
}
static void wire_names(Writer *w,SEXP s) {
    wire_byte(w,s!=R_NilValue);
    if(s==R_NilValue)return;
    if(TYPEOF(s)!=STRSXP)Rf_error("formula names must be text");
    wire_length(w,(uint64_t)XLENGTH(s));for(R_xlen_t i=0;i<XLENGTH(s);++i)wire_text(w,STRING_ELT(s,i));
}
static void wire_element_names(Writer *w,SEXP names,R_xlen_t rows,R_xlen_t cols,int matrix){
    if(!matrix||names==R_NilValue){wire_names(w,names);return;}
    if(TYPEOF(names)!=STRSXP||XLENGTH(names)!=rows*cols)Rf_error("formula names do not match values");
    wire_byte(w,1);wire_length(w,(uint64_t)XLENGTH(names));
    for(R_xlen_t i=0;i<XLENGTH(names);++i)wire_text(w,STRING_ELT(names,i/cols+(i%cols)*rows));
}
static const char *attributes[]={"class","levels","contrasts","scaled:center","scaled:scale","coefs","degree","units","tzone"};
static void wire_value(Writer *w,SEXP v,int depth) {
    if(depth>RSDESEQ2_FORMULA_MAX_DEPTH)Rf_error("formula transport nesting exceeds 1024");
    if(v==R_NilValue){wire_byte(w,0);wire_length(w,0);for(int i=0;i<6;++i)wire_byte(w,0);wire_length(w,0);return;}
    int tag;
    switch(TYPEOF(v)){case NILSXP:tag=0;break;case REALSXP:tag=1;break;case INTSXP:tag=2;break;case LGLSXP:tag=3;break;case STRSXP:tag=4;break;case VECSXP:tag=5;break;default:Rf_error("formula values must be primitive vectors or lists");return;}
    R_xlen_t n=XLENGTH(v),rows=0,cols=0;SEXP dims=Rf_getAttrib(v,R_DimSymbol);
    if(dims!=R_NilValue) {
        if(TYPEOF(dims)!=INTSXP||XLENGTH(dims)!=2||INTEGER(dims)[0]<0||INTEGER(dims)[1]<0)Rf_error("formula values require two-dimensional matrices");
        rows=INTEGER(dims)[0];cols=INTEGER(dims)[1];
        if((uint64_t)rows*(uint64_t)cols!=(uint64_t)n)Rf_error("formula dimensions do not match values");
    }
    wire_byte(w,(uint8_t)tag);wire_length(w,(uint64_t)n);
    /* Force ALTREP storage during C preparation, never from a Rust frame. */
    const double *real=tag==1?REAL(v):NULL;
    const int *integer=tag==2?INTEGER(v):tag==3?LOGICAL(v):NULL;
    for(R_xlen_t i=0;i<n;++i) {
        R_xlen_t at=dims==R_NilValue?i:(i/cols)+(i%cols)*rows;
        if(tag==1||tag==2){double x=tag==1?real[at]:integer[at]==NA_INTEGER?NA_REAL:(double)integer[at];uint64_t bits;memcpy(&bits,&x,8);wire_length(w,bits);}
        else if(tag==3)wire_byte(w,integer[at]==NA_LOGICAL?2:(uint8_t)(integer[at]!=0));
        else if(tag==4){SEXP s=STRING_ELT(v,at);wire_byte(w,s!=NA_STRING);if(s!=NA_STRING)wire_text(w,s);}
        else if(tag==5)wire_value(w,VECTOR_ELT(v,at),depth+1);
    }
    wire_byte(w,dims!=R_NilValue);if(dims!=R_NilValue){wire_length(w,(uint64_t)rows);wire_length(w,(uint64_t)cols);}
    wire_element_names(w,Rf_getAttrib(v,R_NamesSymbol),rows,cols,dims!=R_NilValue);
    SEXP dimnames=Rf_getAttrib(v,R_DimNamesSymbol);
    if(dimnames!=R_NilValue&&(TYPEOF(dimnames)!=VECSXP||XLENGTH(dimnames)!=2))Rf_error("formula dimnames require two entries");
    for(int i=0;i<2;++i)wire_names(w,dimnames==R_NilValue?R_NilValue:VECTOR_ELT(dimnames,i));
    for(int i=0;i<2;++i){SEXP labels=dimnames==R_NilValue?R_NilValue:VECTOR_ELT(dimnames,i);wire_names(w,labels==R_NilValue?R_NilValue:Rf_getAttrib(labels,R_NamesSymbol));}
    size_t count=0;
    for(size_t i=0;i<sizeof(attributes)/sizeof(*attributes);++i)if(Rf_getAttrib(v,Rf_install(attributes[i]))!=R_NilValue)++count;
    wire_length(w,count);
    for(size_t i=0;i<sizeof(attributes)/sizeof(*attributes);++i){SEXP attr=Rf_getAttrib(v,Rf_install(attributes[i]));if(attr!=R_NilValue){size_t n=strlen(attributes[i]);wire_length(w,n);wire_bytes(w,attributes[i],n);wire_value(w,attr,depth+1);}}
}
size_t rsdeseq2_formula_measure(SEXP input){Writer w={NULL,0,RSDESEQ2_FORMULA_MAX_BYTES};wire_bytes(&w,"RFX2",4);wire_value(&w,input,0);return w.at;}
void rsdeseq2_formula_write(SEXP input,uint8_t *output,size_t size){Writer w={output,0,size};wire_bytes(&w,"RFX2",4);wire_value(&w,input,0);if(w.at!=size)Rf_error("formula transport size changed during encoding");}
