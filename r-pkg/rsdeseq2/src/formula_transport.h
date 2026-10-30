#ifndef RSDESEQ2_FORMULA_TRANSPORT_H
#define RSDESEQ2_FORMULA_TRANSPORT_H
#include <Rinternals.h>
#include <stddef.h>
#include <stdint.h>
#define RSDESEQ2_FORMULA_MAX_BYTES (256u * 1024u * 1024u)
#define RSDESEQ2_FORMULA_MAX_DEPTH 1024
/* These codecs only translate primitive storage/metadata. The input pass is
 * complete before entering Rust; output allocation runs under unwind cleanup. */
size_t rsdeseq2_formula_measure(SEXP);
void rsdeseq2_formula_write(SEXP, uint8_t *, size_t);
SEXP rsdeseq2_formula_read(const uint8_t *, size_t);
#ifdef RSDESEQ2_TEST_FORMULA_UNWIND
void rsdeseq2_formula_inject(int);
#endif
#endif
