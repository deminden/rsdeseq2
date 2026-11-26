//! Tailorings: rules applied on top of the root collation.
//!
//! A tailoring is a short list of rules in the usual notation of collation
//! rules: `&x < y` puts `y` right after `x` as a new base letter, `<<` as an
//! accent variant, `<<<` as a case or letter variant and `=` as equal;
//! `&[before 1]x < y` puts `y` right before `x`. A rule item of several
//! characters is a contraction (Czech `ch`), and `x/y` sorts `x` as if `y`
//! followed it. `&[before 2]x << y` puts `y` right before the accent of
//! `x`. A tailoring can also move a script to another place (`reorder`),
//! switch off root contractions (`suppress`), sort uppercase before lowercase
//! (`upper_first`, Danish), compare accents from the end of the string
//! (`backwards`, Canadian French) and put combining marks in canonical order
//! first (`normalize`, Greek). A tailored item takes the case of its
//! characters: uppercase, lowercase or, for a contraction like Danish `Aa`,
//! mixed. Each tailoring's rules and options were read off R's order for
//! that locale; see `docs/provenance/specs/collation.md`.

mod rules;

use std::collections::BTreeMap;
use std::sync::OnceLock;

use crate::root::{
    COMMON, Element, LOWER, MIXED, PLAIN, PRIMARY_STEP, UPPERCASE, char_case, decompose, root,
};

pub(crate) use rules::DEFINITIONS;

/// One tailoring as written in [`rules::DEFINITIONS`].
pub(crate) struct Definition {
    /// The name ICU reports for the collation.
    pub(crate) name: &'static str,
    /// The rules.
    pub(crate) rules: &'static str,
    /// Scripts moved, in order, to the given place: (first atom, first atom
    /// after the script, atom before which the script goes).
    pub(crate) reorder: &'static [(u32, u32, u32)],
    /// Root contractions not used (by the character they stand for).
    pub(crate) suppress: &'static [u32],
    /// Uppercase sorts before lowercase at the tertiary level.
    pub(crate) upper_first: bool,
    /// Secondary weights are compared from the end of the string.
    pub(crate) backwards: bool,
    /// Combining marks are put in canonical order before collation.
    pub(crate) normalize: bool,
}

/// A tailored collation: an index into [`DEFINITIONS`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Tailoring(u8);

struct Compiled {
    /// Tailored sequences (canonically decomposed) and their elements.
    mappings: BTreeMap<Vec<u32>, Vec<Element>>,
    longest: usize,
    touches_ascii: bool,
    /// Primary ranges moved: (first, last, shift as a signed offset).
    moves: Vec<(u32, u32, i64)>,
}

impl Tailoring {
    /// The tailoring with the given ICU name.
    pub(crate) fn find(name: &str) -> Option<Self> {
        let i = DEFINITIONS.iter().position(|d| d.name == name)?;
        u8::try_from(i).ok().map(Self)
    }

    fn definition(self) -> &'static Definition {
        &DEFINITIONS[usize::from(self.0)]
    }

    fn compiled(self) -> &'static Compiled {
        static COMPILED: OnceLock<Vec<OnceLock<Compiled>>> = OnceLock::new();
        let all = COMPILED.get_or_init(|| DEFINITIONS.iter().map(|_| OnceLock::new()).collect());
        all[usize::from(self.0)].get_or_init(|| compile(self.definition()))
    }

    pub(crate) fn name(self) -> &'static str {
        self.definition().name
    }

    pub(crate) fn upper_first(self) -> bool {
        self.definition().upper_first
    }

    pub(crate) fn backwards(self) -> bool {
        self.definition().backwards
    }

    pub(crate) fn normalize(self) -> bool {
        self.definition().normalize
    }

    /// Whether the tailoring changes how some ASCII string sorts against
    /// another (then the ASCII fast path does not apply).
    pub(crate) fn touches_ascii(self) -> bool {
        self.compiled().touches_ascii
    }

    pub(crate) fn push_decomposed(self, seq: &[u32], out: &mut Vec<Element>) {
        let c = self.compiled();
        let start = out.len();
        let def = self.definition();
        push_with(&c.mappings, c.longest, def.suppress, seq, out);
        for e in &mut out[start..] {
            for &(first, last, shift) in &c.moves {
                if e.primary >= first && e.primary <= last {
                    // Shifts keep weights positive and within u32 by construction.
                    e.primary = u32::try_from(i64::from(e.primary) + shift).unwrap_or(e.primary);
                    break;
                }
            }
        }
    }
}

/// Elements of a decomposed sequence with tailored mappings taking
/// precedence over the root collation. A tailored letter written as a base
/// and a mark also matches when other marks of a lower combining class come
/// between them (`"c\u{327}\u{30C}"` is `"č"` with a cedilla).
fn push_with(
    mappings: &BTreeMap<Vec<u32>, Vec<Element>>,
    longest: usize,
    suppress: &[u32],
    seq: &[u32],
    out: &mut Vec<Element>,
) {
    let root = root();
    let mut seq = seq.to_vec();
    let mut i = 0;
    while i < seq.len() {
        let max = longest.min(seq.len() - i);
        let mut tailored = (1..=max)
            .rev()
            .find_map(|len| mappings.get(&seq[i..i + len]).map(|e| (len, e)));
        if tailored.as_ref().is_none_or(|t| t.0 == 1)
            && let Some((j, els)) = discontiguous(mappings, &seq, i)
        {
            // Move the matched mark next to its base and match normally.
            let mark = seq.remove(j);
            seq.insert(i + 1, mark);
            tailored = Some((2, els));
        }
        let rooted = root
            .contraction(&seq[i..])
            .filter(|(_, target)| !suppress.contains(target));
        match (tailored, rooted) {
            (Some((len, els)), Some((rlen, _))) if len >= rlen => {
                out.extend_from_slice(els);
                i += len;
            }
            (Some((len, els)), None) => {
                out.extend_from_slice(els);
                i += len;
            }
            (_, Some((len, target))) => {
                root.push_char(target, out);
                i += len;
            }
            (None, None) => {
                root.push_char(seq[i], out);
                i += 1;
            }
        }
    }
}

/// A mapping of a base and one mark where the mark comes later in the run of
/// marks after the base, not blocked by a mark of the same or a higher
/// combining class: (position of the mark, elements).
fn discontiguous<'a>(
    mappings: &'a BTreeMap<Vec<u32>, Vec<Element>>,
    seq: &[u32],
    i: usize,
) -> Option<(usize, &'a Vec<Element>)> {
    let mut highest = 0u8;
    for j in i + 2..seq.len() {
        let previous = combining_class(seq[j - 1]);
        if previous == 0 {
            return None;
        }
        highest = highest.max(previous);
        let class = combining_class(seq[j]);
        if class == 0 {
            return None;
        }
        if class > highest
            && let Some(els) = mappings.get(&[seq[i], seq[j]][..])
        {
            return Some((j, els));
        }
    }
    None
}

fn combining_class(c: u32) -> u8 {
    char::from_u32(c).map_or(0, unicode_normalization::char::canonical_combining_class)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Relation {
    Primary,
    Secondary,
    Tertiary,
    Equal,
}

/// A relation, the item it places and the item's extension.
type Item = (Relation, Vec<u32>, Vec<u32>);
/// A reset (`[before n]` level and text) and the items that follow it.
type Rule = (Option<u8>, Vec<u32>, Vec<Item>);

/// Parse `rules` into resets, each followed by (relation, item, extension).
fn parse(rules: &str) -> Vec<Rule> {
    let mut out: Vec<Rule> = Vec::new();
    let mut tokens = rules.split_whitespace().peekable();
    while let Some(token) = tokens.next() {
        if let Some(rest) = token.strip_prefix('&') {
            let (before, text) = if rest == "[before" {
                let level = tokens
                    .next()
                    .and_then(|l| l.trim_end_matches(']').parse().ok());
                (level, tokens.next().unwrap_or(""))
            } else {
                (None, rest)
            };
            out.push((before, chars(text), Vec::new()));
            continue;
        }
        let relation = match token {
            "<" => Relation::Primary,
            "<<" => Relation::Secondary,
            "<<<" => Relation::Tertiary,
            "=" => Relation::Equal,
            _ => continue,
        };
        let item = tokens.next().unwrap_or("");
        let (text, extension) = item.split_once('/').unwrap_or((item, ""));
        if let Some(last) = out.last_mut() {
            last.2.push((relation, chars(text), chars(extension)));
        }
    }
    out
}

/// Code points of a rule item; `\uXXXX` stands for an invisible character.
fn chars(text: &str) -> Vec<u32> {
    let mut out = Vec::new();
    let mut rest = text;
    while !rest.is_empty() {
        if let Some(hex) = rest.strip_prefix("\\u") {
            let code = hex.get(..4).and_then(|h| u32::from_str_radix(h, 16).ok());
            if let Some(code) = code {
                out.push(code);
                rest = &hex[4..];
                continue;
            }
        }
        let mut it = rest.chars();
        if let Some(c) = it.next() {
            out.push(u32::from(c));
        }
        rest = it.as_str();
    }
    let text: String = out.iter().filter_map(|&c| char::from_u32(c)).collect();
    decompose(&text)
}

fn compile(def: &Definition) -> Compiled {
    let mut mappings: BTreeMap<Vec<u32>, Vec<Element>> = BTreeMap::new();
    // Inserted primary weights after each root primary, in order; finalised
    // as root weight + position + 1.
    let mut inserted: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    let mut next_id = 0usize;
    // Elements may carry provisional primaries: ids above `PROVISIONAL`.
    const PROVISIONAL: u32 = 0xF000_0000;
    let mut minor_used: BTreeMap<(u32, u16, u16, u8), u16> = BTreeMap::new();
    for (before, reset, items) in parse(def.rules) {
        let mut longest = mappings.keys().map(Vec::len).max().unwrap_or(1);
        let mut current = Vec::new();
        push_with(&mappings, longest, def.suppress, &reset, &mut current);
        let mut before = before;
        for (relation, item, extension) in items {
            let mut new = current.clone();
            let case = item_case(&item);
            let last = new.iter().rposition(|e| e.primary != 0 || e.secondary != 0);
            match (relation, last) {
                (Relation::Equal, _) | (_, None) => {}
                (Relation::Primary, Some(k)) => {
                    let anchor = new[k].primary;
                    let id = next_id;
                    next_id += 1;
                    insert_primary(
                        &mut inserted,
                        anchor,
                        id,
                        before.take() == Some(1),
                        PROVISIONAL,
                    );
                    new.truncate(k + 1);
                    new[k] = Element::new(PROVISIONAL + id as u32, COMMON, PLAIN).with_case(case);
                }
                (Relation::Secondary, Some(k)) => {
                    let e = new[k];
                    let n = minor_used
                        .entry((e.primary, e.secondary, e.tertiary, 2))
                        .or_insert(0);
                    *n += 1;
                    new.truncate(k + 1);
                    let s = if before.take() == Some(2) {
                        e.secondary - *n
                    } else {
                        e.secondary + *n
                    };
                    new[k] = Element::new(e.primary, s, PLAIN).with_case(case);
                }
                (Relation::Tertiary, Some(k)) => {
                    let e = new[k];
                    let n = minor_used
                        .entry((e.primary, e.secondary, e.tertiary, 3))
                        .or_insert(0);
                    *n += 1;
                    new.truncate(k + 1);
                    new[k] = Element::new(e.primary, e.secondary, e.tertiary + *n).with_case(case);
                }
            }
            // The next item is placed relative to this one without its
            // extension.
            current = new.clone();
            if !extension.is_empty() {
                push_with(&mappings, longest, def.suppress, &extension, &mut new);
            }
            mappings.insert(item.clone(), new);
            longest = longest.max(item.len());
        }
    }
    // Final primary weights of the insertions.
    let mut weight = vec![0u32; next_id];
    for (&anchor, ids) in &inserted {
        for (pos, &id) in ids.iter().enumerate() {
            weight[id] = anchor + u32::try_from(pos).unwrap_or(0) + 1;
        }
    }
    let resolve = |p: u32| -> u32 {
        if p >= PROVISIONAL {
            weight.get((p - PROVISIONAL) as usize).copied().unwrap_or(p)
        } else {
            p
        }
    };
    for els in mappings.values_mut() {
        for e in els.iter_mut() {
            e.primary = resolve(e.primary);
        }
    }
    // Only a mapping of an ASCII-only sequence can change how two ASCII
    // strings compare; moving a script does not.
    let touches_ascii =
        def.upper_first || def.backwards || mappings.keys().any(|k| k.iter().all(|&c| c < 128));
    let root = root();
    let mut moves = Vec::new();
    for &(first, end, target) in def.reorder {
        let (Some(lo), Some(hi), Some(to)) = (
            root.atom_primary(first),
            root.atom_primary(end),
            root.atom_primary(target),
        ) else {
            continue;
        };
        // The script [lo, hi) moves to `to`; what lay in [to, lo) moves up.
        let len = i64::from(hi) - i64::from(lo);
        moves.push((lo, hi - 1, i64::from(to) - i64::from(lo)));
        moves.push((to, lo - 1, len));
    }
    let longest = mappings.keys().map(Vec::len).max().unwrap_or(1);
    Compiled {
        mappings,
        longest,
        touches_ascii,
        moves,
    }
}

/// Case of a rule item: uppercase or lowercase when all its cased characters
/// are, mixed otherwise.
fn item_case(item: &[u32]) -> u8 {
    let upper = item.iter().filter(|&&c| char_case(c) == UPPERCASE).count();
    let lower = item
        .iter()
        .filter(|&&c| char::from_u32(c).is_some_and(char::is_lowercase))
        .count();
    match (upper, lower) {
        (0, _) => LOWER,
        (_, 0) => UPPERCASE,
        _ => MIXED,
    }
}

/// Record a new primary right after (or, with `before`, right before) the
/// weight `anchor`, which is a root weight or a provisional one.
fn insert_primary(
    inserted: &mut BTreeMap<u32, Vec<usize>>,
    anchor: u32,
    id: usize,
    before: bool,
    provisional: u32,
) {
    if anchor >= provisional {
        let target = (anchor - provisional) as usize;
        for ids in inserted.values_mut() {
            if let Some(pos) = ids.iter().position(|&x| x == target) {
                ids.insert(if before { pos } else { pos + 1 }, id);
                return;
            }
        }
        return;
    }
    let base = anchor - anchor % PRIMARY_STEP;
    if before && anchor.is_multiple_of(PRIMARY_STEP) {
        inserted.entry(base - PRIMARY_STEP).or_default().push(id);
    } else {
        let ids = inserted.entry(base).or_default();
        let pos = usize::try_from(anchor % PRIMARY_STEP).unwrap_or(0);
        ids.insert(pos.min(ids.len()), id);
    }
}
