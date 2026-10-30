//! Versioned primitive transport. No executable objects or foreign pointers.
//! Lengths are little-endian u64, doubles retain IEEE bits (integer carriers
//! use canonical missing/zero representations), and all
//! matrices are row-major. Atomic tags: NULL=0, double=1, integer=2, logical=3,
//! text=4, list=5. Every value then carries dimensions/names/attributes.
use super::{Data, Value, invalid, validation::validate_value_with_limit};
use crate::errors::DeseqError;

/// Maximum transport payload, distinct from numerical API input sizes.
pub const MAX_BYTES: usize = 256 * 1024 * 1024;
/// Transport nodes include AST argument/list wrappers; their depth limit must
/// accommodate the evaluator's 128 semantic levels plus literal metadata.
pub const MAX_DEPTH: usize = 1024;
const MAGIC: &[u8; 4] = b"RFX2";

/// Serialize validated owned primitive data.
pub fn encode(value: &Value) -> Result<Vec<u8>, DeseqError> {
    validate_value_with_limit(value, 0, MAX_DEPTH)?;
    let mut w = Writer(MAGIC.to_vec());
    w.value(value, 0)?;
    Ok(w.0)
}
/// Decode and validate exactly one primitive value, rejecting trailing bytes.
pub fn decode(bytes: &[u8]) -> Result<Value, DeseqError> {
    if bytes.len() > MAX_BYTES || !bytes.starts_with(MAGIC) {
        return Err(invalid("invalid formula transport header or size"));
    }
    let mut r = Reader { bytes, at: 4 };
    let value = r.value(0)?;
    if r.at != bytes.len() {
        return Err(invalid("trailing formula transport bytes"));
    }
    validate_value_with_limit(&value, 0, MAX_DEPTH)?;
    Ok(value)
}
struct Writer(Vec<u8>);
impl Writer {
    fn put(&mut self, bytes: &[u8]) -> Result<(), DeseqError> {
        if self
            .0
            .len()
            .checked_add(bytes.len())
            .is_none_or(|n| n > MAX_BYTES)
        {
            return Err(invalid("formula transport exceeds size limit"));
        }
        self.0
            .try_reserve(bytes.len())
            .map_err(|e| invalid(e.to_string()))?;
        self.0.extend_from_slice(bytes);
        Ok(())
    }
    fn byte(&mut self, v: u8) -> Result<(), DeseqError> {
        self.put(&[v])
    }
    fn length(&mut self, n: usize) -> Result<(), DeseqError> {
        self.put(&(n as u64).to_le_bytes())
    }
    fn text(&mut self, s: &str) -> Result<(), DeseqError> {
        self.length(s.len())?;
        self.put(s.as_bytes())
    }
    fn names(&mut self, names: &Option<Vec<String>>) -> Result<(), DeseqError> {
        self.byte(u8::from(names.is_some()))?;
        if let Some(n) = names {
            self.length(n.len())?;
            for s in n {
                self.text(s)?;
            }
        }
        Ok(())
    }
    fn value(&mut self, v: &Value, depth: usize) -> Result<(), DeseqError> {
        if depth > MAX_DEPTH {
            return Err(invalid("formula transport nesting exceeds 1024"));
        }
        self.byte(match v.data {
            Data::Null => 0,
            Data::Number { integer: false, .. } => 1,
            Data::Number { integer: true, .. } => 2,
            Data::Logical(_) => 3,
            Data::Text(_) => 4,
            Data::List(_) => 5,
        })?;
        self.length(v.len())?;
        match &v.data {
            Data::Number { values, .. } => {
                for n in values {
                    self.put(&n.to_bits().to_le_bytes())?;
                }
            }
            Data::Logical(values) => {
                for n in values {
                    self.byte(n.map_or(2, u8::from))?;
                }
            }
            Data::Text(values) => {
                for n in values {
                    self.byte(u8::from(n.is_some()))?;
                    if let Some(n) = n {
                        self.text(n)?;
                    }
                }
            }
            Data::List(values) => {
                for n in values {
                    self.value(n, depth + 1)?;
                }
            }
            Data::Null => (),
        }
        self.byte(u8::from(v.dimensions.is_some()))?;
        if let Some([r, c]) = v.dimensions {
            self.length(r)?;
            self.length(c)?;
        }
        self.names(&v.names)?;
        self.names(&v.dimnames[0])?;
        self.names(&v.dimnames[1])?;
        self.names(&v.dimname_names[0])?;
        self.names(&v.dimname_names[1])?;
        self.length(v.attributes.len())?;
        for (name, value) in &v.attributes {
            self.text(name)?;
            self.value(value, depth + 1)?;
        }
        Ok(())
    }
}
struct Reader<'a> {
    bytes: &'a [u8],
    at: usize,
}
impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], DeseqError> {
        let end = self
            .at
            .checked_add(n)
            .filter(|e| *e <= self.bytes.len())
            .ok_or_else(|| invalid("truncated formula transport"))?;
        let result = &self.bytes[self.at..end];
        self.at = end;
        Ok(result)
    }
    fn byte(&mut self) -> Result<u8, DeseqError> {
        Ok(self.take(1)?[0])
    }
    fn flag(&mut self) -> Result<bool, DeseqError> {
        match self.byte()? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(invalid("invalid formula transport flag")),
        }
    }
    fn length(&mut self) -> Result<usize, DeseqError> {
        usize::try_from(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
            .map_err(|e| invalid(e.to_string()))
    }
    fn count(&mut self, unit: usize) -> Result<usize, DeseqError> {
        let n = self.length()?;
        if n.checked_mul(unit)
            .is_none_or(|n| n > self.bytes.len() - self.at)
        {
            return Err(invalid("invalid formula transport length"));
        }
        Ok(n)
    }
    fn text(&mut self) -> Result<String, DeseqError> {
        let n = self.count(1)?;
        String::from_utf8(self.take(n)?.to_vec()).map_err(|e| invalid(e.to_string()))
    }
    fn names(&mut self) -> Result<Option<Vec<String>>, DeseqError> {
        if !self.flag()? {
            return Ok(None);
        }
        let n = self.count(8)?;
        Ok(Some((0..n).map(|_| self.text()).collect::<Result<_, _>>()?))
    }
    fn value(&mut self, depth: usize) -> Result<Value, DeseqError> {
        if depth > MAX_DEPTH {
            return Err(invalid("formula transport nesting exceeds 1024"));
        }
        let tag = self.byte()?;
        let n = self.count(if matches!(tag, 1 | 2) {
            8
        } else if tag == 5 {
            23
        } else {
            1
        })?;
        let data = match tag {
            0 if n == 0 => Data::Null,
            1 | 2 => Data::Number {
                values: (0..n)
                    .map(|_| {
                        Ok(f64::from_bits(u64::from_le_bytes(
                            self.take(8)?.try_into().unwrap(),
                        )))
                    })
                    .collect::<Result<_, DeseqError>>()?,
                integer: tag == 2,
            },
            3 => Data::Logical(
                (0..n)
                    .map(|_| match self.byte()? {
                        0 => Ok(Some(false)),
                        1 => Ok(Some(true)),
                        2 => Ok(None),
                        _ => Err(invalid("invalid formula logical byte")),
                    })
                    .collect::<Result<_, _>>()?,
            ),
            4 => Data::Text(
                (0..n)
                    .map(|_| {
                        if self.flag()? {
                            self.text().map(Some)
                        } else {
                            Ok(None)
                        }
                    })
                    .collect::<Result<_, _>>()?,
            ),
            5 => Data::List(
                (0..n)
                    .map(|_| self.value(depth + 1))
                    .collect::<Result<_, _>>()?,
            ),
            _ => return Err(invalid("invalid formula transport value tag")),
        };
        let mut v = Value::new(data);
        if self.flag()? {
            v.dimensions = Some([self.length()?, self.length()?]);
        }
        v.names = self.names()?;
        v.dimnames = [self.names()?, self.names()?];
        v.dimname_names = [self.names()?, self.names()?];
        let n = self.count(9)?;
        v.attributes = (0..n)
            .map(|_| Ok((self.text()?, self.value(depth + 1)?)))
            .collect::<Result<_, DeseqError>>()?;
        Ok(v)
    }
}
