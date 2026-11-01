// Arithmetic schedules adapted from OpenBLAS 0.3.26.
// See OPENBLAS-LICENSE in this directory for its BSD license.
//! Linux x86 R uses the x87 extended accumulator for scalar reductions.
//! Keep that arithmetic in the Rust core, including the final rounding to f64.
use std::arch::asm;

#[derive(Clone, Copy)]
pub(super) struct Extended([u8; 16]);
impl Extended {
    pub(super) fn new(value: f64) -> Self {
        let mut out = Self([0; 16]);
        // SAFETY: both pointers address valid storage of the instruction's size;
        // each block balances its x87 pushes and pops and preserves control state.
        unsafe {
            asm!("fld qword ptr [{input}]", "fstp tbyte ptr [{output}]", input=in(reg) &value, output=in(reg) out.0.as_mut_ptr(), out("st(0)") _, options(nostack, preserves_flags));
        }
        out
    }
    pub(super) fn add(&mut self, value: Self) {
        unsafe {
            asm!("fld tbyte ptr [{input}]", "fld tbyte ptr [{other}]", "faddp st(1), st(0)", "fstp tbyte ptr [{input}]", input=in(reg) self.0.as_mut_ptr(), other=in(reg) value.0.as_ptr(), out("st(0)") _, out("st(1)") _, options(nostack, preserves_flags));
        }
    }
    pub(super) fn divide(&mut self, value: f64) {
        unsafe {
            asm!("fld tbyte ptr [{input}]", "fdiv qword ptr [{other}]", "fstp tbyte ptr [{input}]", input=in(reg) self.0.as_mut_ptr(), other=in(reg) &value, out("st(0)") _, options(nostack, preserves_flags));
        }
    }
    pub(super) fn difference(value: f64, other: Self) -> Self {
        let mut result = Self([0; 16]);
        unsafe {
            asm!("fld qword ptr [{value}]", "fld tbyte ptr [{other}]", "fsubp st(1), st(0)", "fstp tbyte ptr [{output}]", value=in(reg) &value, other=in(reg) other.0.as_ptr(), output=in(reg) result.0.as_mut_ptr(), out("st(0)") _, out("st(1)") _, options(nostack, preserves_flags));
        }
        result
    }
    fn square(value: f64) -> Self {
        let mut out = Self([0; 16]);
        unsafe {
            asm!("fld qword ptr [{input}]", "fmul st(0), st(0)", "fstp tbyte ptr [{output}]", input=in(reg) &value, output=in(reg) out.0.as_mut_ptr(), out("st(0)") _, options(nostack, preserves_flags));
        }
        out
    }
    fn sqrt(mut self) -> Self {
        unsafe {
            asm!("fld tbyte ptr [{input}]", "fsqrt", "fstp tbyte ptr [{input}]", input=in(reg) self.0.as_mut_ptr(), out("st(0)") _, options(nostack, preserves_flags));
        }
        self
    }
    pub(super) fn value(self) -> f64 {
        let mut out = 0.0;
        unsafe {
            asm!("fld tbyte ptr [{input}]", "fstp qword ptr [{output}]", input=in(reg) self.0.as_ptr(), output=in(reg) &mut out, out("st(0)") _, options(nostack, preserves_flags));
        }
        out
    }
}
pub(super) fn sum(values: impl Iterator<Item = f64>) -> Extended {
    let mut total = Extended::new(0.0);
    for value in values {
        total.add(Extended::new(value));
    }
    total
}
pub(super) fn mean(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let mut result = sum(values.iter().copied());
    result.divide(n);
    if result.value().is_finite() {
        let mut correction = Extended::new(0.0);
        for &value in values {
            correction.add(Extended::difference(value, result));
        }
        correction.divide(n);
        result.add(correction);
    }
    result.value()
}

pub(super) fn norm(values: &[f64]) -> f64 {
    let mut lanes = [Extended::new(0.0); 4];
    let end = values.len() / 8 * 8;
    for (i, &x) in values[..end].iter().enumerate() {
        lanes[i % 4].add(Extended::square(x));
    }
    for &x in &values[end..] {
        lanes[0].add(Extended::square(x));
    }
    let mut total = lanes[2];
    total.add(lanes[0]);
    total.add(lanes[1]);
    total.add(lanes[3]);
    total.sqrt().value()
}
