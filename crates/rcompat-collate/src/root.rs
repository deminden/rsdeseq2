//! The root collation: rules plus the data in [`crate::tables`].
//!
//! Each character of a string is replaced by its canonical decomposition
//! (without reordering combining marks, as R does), prevowels are moved
//! after the consonant they precede, contractions are matched, and every
//! remaining code point is mapped by the first rule that applies:
//!
//! 1. code points R sorts as unassigned: code point order after everything;
//! 2. stored characters (`EXPLICIT`, `VARIANTS`);
//! 3. atoms: one primary weight each, in R's order (`ATOM_RUNS`, the
//!    ideographs in `HAN_*`), common secondary, plain tertiary;
//! 4. marks: a secondary weight only (`MARKS`);
//! 5. completely ignorable code points: nothing;
//! 6. compatibility decompositions: the elements of the decomposition with the
//!    tertiary weight of the variant (`COMPAT`);
//! 7. uppercase letters: the elements of the lowercase letter with the
//!    uppercase tertiary weight;
//! 8. everything else (unassigned and private use): code point order after
//!    all atoms.

use std::collections::HashSet;
use std::sync::OnceLock;

use unicode_normalization::UnicodeNormalization;
use unicode_normalization::char::{
    canonical_combining_class, decompose_canonical, decompose_compatible,
};

use crate::tables;

/// Gap between consecutive root primary weights, so that tailorings can put
/// new weights between them.
pub(crate) const PRIMARY_STEP: u32 = 1 << 8;
/// Gap between consecutive secondary and tertiary codes.
pub(crate) const MINOR_STEP: u16 = 1 << 5;
/// The common secondary weight.
pub(crate) const COMMON: u16 = MINOR_STEP;
/// The tertiary weight of plain letters.
pub(crate) const PLAIN: u16 = tables::PLAIN * MINOR_STEP;
/// The tertiary weight of uppercase letters.
pub(crate) const UPPER: u16 = tables::UPPER * MINOR_STEP;
/// The first primary weight after every atom: unassigned code points follow
/// in code point order, then U+FFFD and U+FFFF.
const TAIL: u32 = (tables::ATOM_COUNT + 1) * PRIMARY_STEP;

/// Case of an element, used only by collations that sort uppercase first:
/// lowercase or uncased, mixed (a tailored contraction such as Danish `Aa`)
/// or uppercase.
pub(crate) const LOWER: u8 = 0;
pub(crate) const MIXED: u8 = 1;
pub(crate) const UPPERCASE: u8 = 2;

/// One collation element. A weight of zero means "no weight at this level".
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Element {
    pub(crate) primary: u32,
    pub(crate) secondary: u16,
    pub(crate) tertiary: u16,
    pub(crate) case: u8,
}

impl Element {
    pub(crate) const fn new(primary: u32, secondary: u16, tertiary: u16) -> Self {
        Self {
            primary,
            secondary,
            tertiary,
            case: LOWER,
        }
    }

    pub(crate) const fn with_case(self, case: u8) -> Self {
        Self { case, ..self }
    }
}

/// Case of a character for the elements it maps to directly: uppercase when
/// it is an uppercase letter or has a lowercase form.
pub(crate) fn char_case(c: u32) -> u8 {
    let upper = char::from_u32(c).is_some_and(char::is_uppercase) || lowercase(c).is_some();
    if upper { UPPERCASE } else { LOWER }
}

/// Symbols of capital letters that R's uppercase-first order shows to be
/// uppercase (the peseta sign "Pts", U+214D "A/S" and the squared words of
/// capitals U+1F18A..1F19A and U+1F1A6..1F1AD, which have no decomposition)
/// or lowercase (the superscript signs SM, TM, MC, MD and MR, whose
/// decompositions are capitals).
const UPPERCASE_SYMBOLS: [(u32, u32); 4] = [
    (0x20A7, 0x20A7),
    (0x214D, 0x214D),
    (0x1F18A, 0x1F19A),
    (0x1F1A6, 0x1F1AD),
];
const LOWERCASE_SYMBOLS: [(u32, u32); 3] = [(0x2120, 0x2120), (0x2122, 0x2122), (0x1F16A, 0x1F16C)];

/// Small kana: lowercase under uppercase-first collations, where the other
/// kana and the halfwidth forms are uppercase.
const SMALL_KANA: [(u32, u32); 28] = [
    (0x3041, 0x3041),
    (0x3043, 0x3043),
    (0x3045, 0x3045),
    (0x3047, 0x3047),
    (0x3049, 0x3049),
    (0x3063, 0x3063),
    (0x3083, 0x3083),
    (0x3085, 0x3085),
    (0x3087, 0x3087),
    (0x308E, 0x308E),
    (0x3095, 0x3096),
    (0x30A1, 0x30A1),
    (0x30A3, 0x30A3),
    (0x30A5, 0x30A5),
    (0x30A7, 0x30A7),
    (0x30A9, 0x30A9),
    (0x30C3, 0x30C3),
    (0x30E3, 0x30E3),
    (0x30E5, 0x30E5),
    (0x30E7, 0x30E7),
    (0x30EE, 0x30EE),
    (0x30F5, 0x30F6),
    (0x31F0, 0x31FF),
    (0xFF67, 0xFF6F),
    (0x1B132, 0x1B132),
    (0x1B150, 0x1B152),
    (0x1B155, 0x1B155),
    (0x1B164, 0x1B167),
];

/// Kana and halfwidth forms whose case R's uppercase-first order shows to be
/// uppercase without a lowercase form.
const CASED_FORMS: [(u32, u32); 7] = [
    (0x3041, 0x3096),
    (0x30A1, 0x30FA),
    (0x31F0, 0x31FF),
    (0xFF61, 0xFF9F),
    (0xFFA0, 0xFFDC),
    (0xFFE8, 0xFFEE),
    (0x1B000, 0x1B16F),
];

/// Case of every element of a kana, a halfwidth form or a symbol listed
/// above, which the rules of [`char_case`] and of the decompositions do not
/// give.
fn form_case(c: u32) -> Option<u8> {
    // Every listed code point is at or after the peseta sign.
    if c < 0x20A7 {
        return None;
    }
    if in_ranges(&UPPERCASE_SYMBOLS, c) {
        return Some(UPPERCASE);
    }
    if in_ranges(&LOWERCASE_SYMBOLS, c) {
        return Some(LOWER);
    }
    if !in_ranges(&CASED_FORMS, c) {
        return None;
    }
    Some(if in_ranges(&SMALL_KANA, c) {
        LOWER
    } else {
        UPPERCASE
    })
}

/// Lookup structures built once from the tables.
pub(crate) struct Root {
    /// (first code point, length, index of the first atom), by code point.
    atoms: Vec<(u32, u32, u32)>,
    /// (code point, secondary, tertiary), by code point.
    marks: Vec<(u32, u16, u16)>,
    /// The secondary weight above every mark.
    top: u16,
    /// Contractions: (sequence, character), by sequence.
    contractions: Vec<(Vec<u32>, u32)>,
    longest_contraction: usize,
    /// Elements of the ASCII characters when they stand alone.
    ascii: [Option<Element>; 128],
}

pub(crate) fn root() -> &'static Root {
    static ROOT: OnceLock<Root> = OnceLock::new();
    ROOT.get_or_init(Root::build)
}

fn in_ranges(ranges: &[(u32, u32)], c: u32) -> bool {
    let i = ranges.partition_point(|r| r.0 <= c);
    i > 0 && c <= ranges[i - 1].1
}

fn tail(c: u32) -> Element {
    let offset = match c {
        0xFFFD => 0x11_0000,
        0xFFFF => 0x11_0001,
        _ => c,
    };
    Element::new(TAIL + offset, COMMON, PLAIN)
}

impl Root {
    fn build() -> Self {
        let mut atoms: Vec<(u32, u32, u32)> = Vec::with_capacity(tables::ATOM_RUNS.len() + 25_000);
        let mut index = 0u32;
        for run in tables::ATOM_RUNS {
            let len = (run & 0x7FF) + 1;
            atoms.push((run >> 11, len, index));
            index += len;
        }
        // Ideographs: runs continue where their block stopped, skipping the
        // characters listed in R's order in HAN_MOVED.
        let mut moved = HashSet::with_capacity(tables::HAN_MOVED.len());
        let mut next_moved = 0;
        for entry in tables::HAN_RUNS {
            if entry & 0x8000 != 0 {
                let (first, _) = tables::HAN_BLOCKS[usize::from(entry & 0x1F)];
                moved.insert(first + u32::from(tables::HAN_MOVED[next_moved]));
                next_moved += 1;
            }
        }
        let mut cursor: Vec<u32> = tables::HAN_BLOCKS.iter().map(|b| b.0).collect();
        next_moved = 0;
        for entry in tables::HAN_RUNS {
            let block = usize::from(entry & 0x1F);
            if entry & 0x8000 != 0 {
                let c = tables::HAN_BLOCKS[block].0 + u32::from(tables::HAN_MOVED[next_moved]);
                next_moved += 1;
                atoms.push((c, 1, index));
                index += 1;
                continue;
            }
            let block = usize::from(entry >> 10);
            let mut left = u32::from(entry & 0x3FF) + 1;
            while left > 0 {
                while moved.contains(&cursor[block]) {
                    cursor[block] += 1;
                }
                let start = cursor[block];
                let mut len = 0;
                while len < left && !moved.contains(&(start + len)) {
                    len += 1;
                }
                atoms.push((start, len, index));
                index += len;
                cursor[block] = start + len;
                left -= len;
            }
        }
        atoms.sort_unstable_by_key(|a| a.0);
        let mut marks = Vec::with_capacity(tables::MARKS.len());
        let mut secondary = (2 + tables::OWN_SECONDARIES) * MINOR_STEP;
        for (i, entry) in tables::MARKS.iter().enumerate() {
            if i > 0 && entry & 0x80 != 0 {
                secondary += MINOR_STEP;
            }
            // The tertiary code has 7 bits; the mask makes the narrowing exact.
            let tertiary = (entry & 0x7F) as u16 * MINOR_STEP;
            marks.push((entry >> 8, secondary, tertiary));
        }
        marks.sort_unstable_by_key(|m| m.0);
        let top = secondary + MINOR_STEP;
        let mut contractions: Vec<(Vec<u32>, u32)> = tables::CONTRACTING
            .iter()
            .filter_map(|&c| {
                let ch = char::from_u32(c)?;
                Some((ch.to_string().nfd().map(u32::from).collect(), c))
            })
            .collect();
        contractions.extend(
            tables::PAIR_CONTRACTIONS
                .iter()
                .map(|&(a, b, c)| (vec![a, b], c)),
        );
        contractions.sort_unstable();
        let longest_contraction = contractions.iter().map(|c| c.0.len()).max().unwrap_or(1);
        let mut root = Self {
            atoms,
            marks,
            top,
            contractions,
            longest_contraction,
            ascii: [None; 128],
        };
        for c in 0..128u32 {
            let mut out = Vec::new();
            root.push_char(c, &mut out);
            root.ascii[c as usize] = out.first().copied();
        }
        root
    }

    pub(crate) fn atom_primary(&self, c: u32) -> Option<u32> {
        let i = self.atoms.partition_point(|a| a.0 <= c);
        let (first, len, index) = *self.atoms.get(i.checked_sub(1)?)?;
        (c < first + len).then(|| (index + c - first + 1) * PRIMARY_STEP)
    }

    fn mark(&self, c: u32) -> Option<(u16, u16)> {
        let i = self.marks.binary_search_by_key(&c, |m| m.0).ok()?;
        Some((self.marks[i].1, self.marks[i].2))
    }

    /// Element of an encoded table entry: `reference << 11 | secondary << 6 |
    /// tertiary`; reference 0 stands for the secondary weight above every
    /// mark.
    fn decode(&self, value: u32) -> Element {
        let reference = value >> 11;
        // Both fields are masked to fit.
        let secondary = ((value >> 6) & 0x1F) as u16;
        let tertiary = (value & 0x3F) as u16 * MINOR_STEP;
        match self.atom_primary(reference) {
            Some(primary) => Element::new(primary, (1 + secondary) * MINOR_STEP, tertiary),
            None if reference == 0 => Element::new(0, self.top, tertiary),
            None => {
                let (s, _) = self.mark(reference).unwrap_or((COMMON, PLAIN));
                Element::new(0, s, tertiary)
            }
        }
    }

    fn stored(&self, c: u32, out: &mut Vec<Element>) -> bool {
        let i = tables::EXPLICIT.partition_point(|e| e.0 <= c);
        if i > 0 && i < tables::EXPLICIT.len() && tables::EXPLICIT[i - 1].0 == c {
            let start = usize::from(tables::EXPLICIT[i - 1].1);
            let end = usize::from(tables::EXPLICIT[i].1);
            let case = char_case(c);
            out.extend(tables::ELEMENTS[start..end].iter().map(|&v| {
                let e = self.decode(v);
                if e.primary != 0 { e.with_case(case) } else { e }
            }));
            return true;
        }
        let i = tables::VARIANTS.partition_point(|v| v.0 >> 11 <= c);
        if let Some(&(run, element)) = i.checked_sub(1).and_then(|i| tables::VARIANTS.get(i)) {
            let offset = c - (run >> 11);
            if offset <= run & 0x7FF {
                let e = self.decode(element + (offset << 11));
                out.push(if e.primary != 0 {
                    e.with_case(char_case(c))
                } else {
                    e
                });
                return true;
            }
        }
        false
    }

    /// Push the elements of one code point of a decomposed string.
    pub(crate) fn push_char(&self, c: u32, out: &mut Vec<Element>) {
        let start = out.len();
        self.push_char_elements(c, out);
        if let Some(case) = form_case(c) {
            for e in &mut out[start..] {
                e.case = case;
            }
        }
    }

    fn push_char_elements(&self, c: u32, out: &mut Vec<Element>) {
        if c < 128
            && let Some(e) = self.ascii.get(c as usize).copied().flatten()
        {
            out.push(e);
            return;
        }
        if in_ranges(&tables::UNASSIGNED, c) {
            out.push(tail(c));
        } else if self.stored(c, out) {
        } else if let Some(primary) = self.atom_primary(c) {
            out.push(Element::new(primary, COMMON, PLAIN).with_case(char_case(c)));
        } else if let Some((s, t)) = self.mark(c) {
            out.push(Element::new(0, s, t));
        } else if in_ranges(&tables::IGNORABLE, c) {
        } else if let Some(parts) = compatibility(c) {
            let start = out.len();
            self.push_decomposed(&parts, out);
            let i = tables::COMPAT.partition_point(|e| e >> 10 <= c);
            let entry = i.checked_sub(1).map_or(0, |i| tables::COMPAT[i]);
            // Five-bit fields.
            let lower = ((entry >> 5) & 0x1F) as u16 * MINOR_STEP;
            let upper = (entry & 0x1F) as u16 * MINOR_STEP;
            for e in &mut out[start..] {
                if e.primary != 0 && e.tertiary == PLAIN && lower != 0 {
                    e.tertiary = lower;
                } else if e.primary != 0 && e.tertiary == UPPER && upper != 0 {
                    e.tertiary = upper;
                }
            }
        } else if let Some(lower) = lowercase(c) {
            let start = out.len();
            self.push_char(lower, out);
            for e in &mut out[start..] {
                if e.primary != 0 {
                    e.case = UPPERCASE;
                    if e.tertiary == PLAIN {
                        e.tertiary = UPPER;
                    }
                }
            }
        } else {
            out.push(tail(c));
        }
    }

    /// Elements of a sequence already in canonical decomposed form.
    pub(crate) fn push_decomposed(&self, seq: &[u32], out: &mut Vec<Element>) {
        let seq = reorder_prevowels(seq);
        let mut i = 0;
        while i < seq.len() {
            if let Some((len, target)) = self.contraction(&seq[i..]) {
                self.push_char(target, out);
                i += len;
            } else {
                self.push_char(seq[i], out);
                i += 1;
            }
        }
    }

    /// The longest contraction at the start of `seq`.
    pub(crate) fn contraction(&self, seq: &[u32]) -> Option<(usize, u32)> {
        let first = *seq.first()?;
        let lo = self.contractions.partition_point(|c| c.0[0] < first);
        if self.contractions.get(lo).is_none_or(|c| c.0[0] != first) {
            return None;
        }
        let longest = self.longest_contraction.min(seq.len());
        (2..=longest).rev().find_map(|len| {
            let key = &seq[..len];
            self.contractions[lo..]
                .binary_search_by(|c| c.0.as_slice().cmp(key))
                .ok()
                .map(|j| (len, self.contractions[lo + j].1))
        })
    }

    /// Elements of an ASCII byte, `None` for completely ignorable ones.
    pub(crate) fn ascii(&self, b: u8) -> Option<Element> {
        self.ascii.get(usize::from(b)).copied().flatten()
    }
}

/// A prevowel followed by one of its consonants sorts after the consonant.
fn reorder_prevowels(seq: &[u32]) -> std::borrow::Cow<'_, [u32]> {
    let is_pre = |c: u32| tables::PREVOWELS.iter().any(|p| p.0 == c);
    if !seq.iter().any(|&c| is_pre(c)) {
        return std::borrow::Cow::Borrowed(seq);
    }
    let mut out = seq.to_vec();
    let mut i = 0;
    while i + 1 < out.len() {
        let (x, y) = (out[i], out[i + 1]);
        if tables::PREVOWELS
            .iter()
            .any(|p| p.0 == x && p.1 <= y && y <= p.2)
        {
            out.swap(i, i + 1);
            i += 2;
        } else {
            i += 1;
        }
    }
    std::borrow::Cow::Owned(out)
}

/// The compatibility decomposition of a code point, in canonical order, when
/// it differs from the code point.
fn compatibility(c: u32) -> Option<Vec<u32>> {
    let ch = char::from_u32(c)?;
    let mut parts = Vec::new();
    decompose_compatible(ch, |x| parts.push(x));
    if parts.len() == 1 && parts[0] == ch {
        return None;
    }
    canonical_order(&mut parts);
    Some(parts.into_iter().map(u32::from).collect())
}

/// Sort each run of combining marks by combining class (stable).
fn canonical_order(chars: &mut [char]) {
    let mut start = 0;
    while start < chars.len() {
        if canonical_combining_class(chars[start]) == 0 {
            start += 1;
            continue;
        }
        let mut end = start;
        while end < chars.len() && canonical_combining_class(chars[end]) != 0 {
            end += 1;
        }
        chars[start..end].sort_by_key(|&c| canonical_combining_class(c));
        start = end;
    }
}

/// The lowercase letter of an uppercase one, when it is a single character.
fn lowercase(c: u32) -> Option<u32> {
    let ch = char::from_u32(c)?;
    let mut lower = ch.to_lowercase();
    let l = lower.next()?;
    (lower.next().is_none() && l != ch).then_some(u32::from(l))
}

/// The code points of `text` with each character replaced by its canonical
/// decomposition. Like R (ICU with normalization off), combining marks are
/// not reordered: `"a\u{301}\u{325}"` and `"a\u{325}\u{301}"` differ.
pub(crate) fn decompose(text: &str) -> Vec<u32> {
    let mut out = Vec::with_capacity(text.len());
    for ch in text.chars() {
        decompose_canonical(ch, |x| out.push(u32::from(x)));
    }
    out
}

/// The canonical decomposition of `text` with combining marks in canonical
/// order (NFD), for collations that normalize.
pub(crate) fn decompose_ordered(text: &str) -> Vec<u32> {
    text.nfd().map(u32::from).collect()
}
