//! The tailorings, as rules read off R's order for each locale
//! (`scripts/provenance/probe_collation_tailorings.R`, checked against every
//! code point and the expansions with `probe_collation_supported.R`; the
//! derivation is described in `docs/provenance/specs/collation.md`).
//!
//! `&[before 1] ǀ` is the end of the Latin letters: the first letter after
//! them in the root order is Greek alpha. Danish sorts uppercase first,
//! Canadian French compares accents from the end of the string, and Greek
//! puts combining marks in canonical order and its script before Latin.

use super::Definition;

/// First Cyrillic letter, first letter after the Cyrillic group (Glagolitic),
/// first Greek letter (the first letter after the Latin group), first letter
/// after the Greek group (Coptic alfa) and the first Latin letter.
const CYRILLIC: u32 = 0x0430;
const AFTER_CYRILLIC: u32 = 0x2C30;
const GREEK: u32 = 0x03B1;
const AFTER_GREEK: u32 = 0x2C81;
const LATIN: u32 = 0x0061;

const fn plain(name: &'static str, rules: &'static str) -> Definition {
    Definition {
        name,
        rules,
        reorder: &[],
        suppress: &[],
        upper_first: false,
        backwards: false,
        normalize: false,
    }
}

const CYRILLIC_FIRST: &[(u32, u32, u32)] = &[(CYRILLIC, AFTER_CYRILLIC, LATIN)];
const CYRILLIC_AFTER_LATIN: &[(u32, u32, u32)] = &[(CYRILLIC, AFTER_CYRILLIC, GREEK)];

const fn with(
    name: &'static str,
    rules: &'static str,
    reorder: &'static [(u32, u32, u32)],
) -> Definition {
    Definition {
        name,
        rules,
        reorder,
        suppress: &[],
        upper_first: false,
        backwards: false,
        normalize: false,
    }
}

const TURKISH: &str = "&C < ç <<< Ç &G < ğ <<< Ğ &[before 1] i < ı <<< I &i <<< İ \
                       &O < ö <<< Ö &S < ş <<< Ş &U < ü <<< Ü";
const CROATIAN: &str = "&C < č <<< Č < ć <<< Ć &D < dž <<< ǆ <<< Dž <<< ǅ <<< DŽ <<< Ǆ < đ <<< Đ \
                        &L < lj <<< ǉ <<< Lj <<< ǈ <<< LJ <<< Ǉ &N < nj <<< ǌ <<< Nj <<< ǋ <<< NJ <<< Ǌ \
                        &S < š <<< Š &Z < ž <<< Ž";
const NORWEGIAN: &str = "&D << đ <<< Đ << ð <<< Ð &t <<< þ/h &T <<< Þ/h &y << ü <<< Ü << ű <<< Ű \
                         &[before 1] ǀ < æ <<< Æ << ä <<< Ä << ę <<< Ę < ø <<< Ø << ö <<< Ö << ő <<< Ő \
                         << œ <<< Œ < å <<< Å << aa <<< Aa <<< AA";

/// All tailorings, by ICU name.
pub(crate) static DEFINITIONS: [Definition; 29] = [
    plain(
        "sv",
        "&D << đ <<< Đ << ð <<< Ð &t <<< þ/h &T <<< Þ/h &y << ü <<< Ü << ű <<< Ű \
         &[before 1] ǀ < å <<< Å < ä <<< Ä << æ <<< Æ << ę <<< Ę < ö <<< Ö << ø <<< Ø << ő <<< Ő \
         << œ <<< Œ << ô <<< Ô",
    ),
    plain(
        "fi",
        "&n\\u0335 << ŋ <<< Ŋ &t\\u0335 << ŧ <<< Ŧ &g\\u0335 << ǥ <<< Ǥ &Y << ü <<< Ü &z\\u0335 << ʒ <<< Ʒ \
         &[before 1] ǀ < å <<< Å < ä <<< Ä << æ <<< Æ < ö <<< Ö << ø <<< Ø",
    ),
    plain("no", NORWEGIAN),
    plain("es", "&N < ñ <<< Ñ"),
    plain("gl", "&N < ñ <<< Ñ"),
    plain(
        "pl",
        "&A < ą <<< Ą &C < ć <<< Ć &E < ę <<< Ę &L < ł <<< Ł &N < ń <<< Ń &O < ó <<< Ó \
         &S < ś <<< Ś &Z < ź <<< Ź < ż <<< Ż",
    ),
    plain(
        "cs",
        "&C < č <<< Č &H < ch <<< cH <<< Ch <<< CH &R < ř <<< Ř &S < š <<< Š &Z < ž <<< Ž",
    ),
    plain(
        "sk",
        "&A < ä <<< Ä &C < č <<< Č &H < ch <<< cH <<< Ch <<< CH &O < ô <<< Ô &R < ř <<< Ř \
         &S < š <<< Š &Z < ž <<< Ž",
    ),
    plain(
        "hu",
        "&C < cs <<< Cs <<< CS &D < dz <<< Dz <<< DZ < dzs <<< Dzs <<< DZS &G < gy <<< Gy <<< GY \
         &L < ly <<< Ly <<< LY &N < ny <<< Ny <<< NY &O < ö <<< Ö << ő <<< Ő &S < sz <<< Sz <<< SZ \
         &T < ty <<< Ty <<< TY &U < ü <<< Ü << ű <<< Ű &Z < zs <<< Zs <<< ZS \
         &cs <<< ccs/cs <<< Ccs/cs <<< CCS/CS &dz <<< ddz/dz <<< Ddz/dz <<< DDZ/DZ \
         &dzs <<< ddzs/dzs <<< Ddzs/dzs <<< DDZS/DZS &gy <<< ggy/gy <<< Ggy/gy <<< GGY/GY \
         &ly <<< lly/ly <<< Lly/ly <<< LLY/LY &ny <<< nny/ny <<< Nny/ny <<< NNY/NY \
         &sz <<< ssz/sz <<< Ssz/sz <<< SSZ/SZ &ty <<< tty/ty <<< Tty/ty <<< TTY/TY \
         &zs <<< zzs/zs <<< Zzs/zs <<< ZZS/ZS",
    ),
    with("hr", CROATIAN, CYRILLIC_AFTER_LATIN),
    with("bs", CROATIAN, CYRILLIC_AFTER_LATIN),
    plain(
        "sl",
        "&C < č <<< Č < ć <<< Ć &D < đ <<< Đ &S < š <<< Š &Z < ž <<< Ž",
    ),
    plain(
        "ro",
        "&A < ă <<< Ă < â <<< Â &I < î <<< Î &S < ş = ș <<< Ş = Ș &T < ţ = ț <<< Ţ = Ț",
    ),
    plain(
        "lt",
        "&A << ą <<< Ą &C < č <<< Č &E << ę <<< Ę << ė <<< Ė &I << į <<< Į << y <<< Y \
         &S < š <<< Š &U << ų <<< Ų << ū <<< Ū &Z < ž <<< Ž",
    ),
    plain(
        "lv",
        "&A < ā <<< Ā &C < č <<< Č &E < ē <<< Ē &G < ģ <<< Ģ &I < y <<< Y < ī <<< Ī \
         &K < ķ <<< Ķ &L < ļ <<< Ļ &N < ņ <<< Ņ &O < ō <<< Ō &R < ŗ <<< Ŗ &S < š <<< Š \
         &U < ū <<< Ū &Z < ž <<< Ž",
    ),
    plain(
        "et",
        "&[before 1] t < š <<< Š < z <<< Z < ž <<< Ž &[before 1] x < õ <<< Õ < ä <<< Ä < ö <<< Ö < ü <<< Ü",
    ),
    plain("tr", TURKISH),
    with(
        "az",
        "&C < ç <<< Ç &E < ə <<< Ə &G < ğ <<< Ğ &H < x <<< X &[before 1] i < ı <<< I \
         &i <<< İ &K < q <<< Q &O < ö <<< Ö &S < ş <<< Ş &U < ü <<< Ü &Z < w <<< W",
        CYRILLIC_AFTER_LATIN,
    ),
    with("ru", "", CYRILLIC_FIRST),
    with("uk", "&г < ґ <<< Ґ &ꙇ < ї <<< ꙶ <<< Ї", CYRILLIC_FIRST),
    with("be", "&е < ё <<< Ё &у < ў <<< Ў", CYRILLIC_FIRST),
    with("bg", "", CYRILLIC_FIRST),
    Definition {
        name: "sr",
        rules: "",
        reorder: CYRILLIC_FIRST,
        suppress: &[0x0419, 0x0439],
        upper_first: false,
        backwards: false,
        normalize: false,
    },
    Definition {
        name: "mk",
        rules: "&ԃ < ѓ <<< Ѓ &[before 1] у < ќ <<< Ќ",
        reorder: CYRILLIC_FIRST,
        suppress: &[0x0419, 0x0439],
        upper_first: false,
        backwards: false,
        normalize: false,
    },
    plain(
        "de@collation=phonebook",
        "&AE << ä <<< Ä &OE << ö <<< Ö &UE << ü <<< Ü",
    ),
    Definition {
        name: "da",
        rules: "&D << đ <<< Đ << ð <<< Ð &oe << œ <<< Œ &t <<< þ/h &T <<< Þ/h &y << ü <<< Ü << ű <<< Ű \
                &[before 1] ǀ < æ <<< Æ << ä <<< Ä < ø <<< Ø << ö <<< Ö << ő <<< Ő \
                < å <<< Å <<< aa <<< Aa <<< AA",
        reorder: &[],
        suppress: &[],
        upper_first: true,
        backwards: false,
        normalize: false,
    },
    plain(
        "is",
        "&[before 1] b < á <<< Á &d << đ <<< Đ < ð <<< Ð &[before 1] f < é <<< É &[before 1] j < í <<< Í \
         &[before 1] p < ó <<< Ó &[before 1] v < ú <<< Ú &[before 1] z < ý <<< Ý \
         &[before 1] ǀ < æ <<< Æ << ä <<< Ä < ö <<< Ö << ø <<< Ø < å <<< Å",
    ),
    Definition {
        name: "el",
        rules: "",
        reorder: &[(GREEK, AFTER_GREEK, LATIN)],
        suppress: &[],
        upper_first: false,
        backwards: false,
        normalize: true,
    },
    Definition {
        name: "fr_CA",
        rules: "",
        reorder: &[],
        suppress: &[],
        upper_first: false,
        backwards: true,
        normalize: false,
    },
];
