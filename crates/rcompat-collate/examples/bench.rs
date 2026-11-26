//! Time comparisons of short strings (no benchmark harness needed).
//!
//! ```sh
//! cargo run --release -p rcompat-collate --example bench
//! ```

use std::hint::black_box;
use std::time::Instant;

use rcompat_collate::Collator;

fn time(name: &str, collator: &Collator, a: &str, b: &str) {
    let n = 1_000_000;
    let start = Instant::now();
    for _ in 0..n {
        black_box(collator.compare(black_box(a), black_box(b)));
    }
    let ns = start.elapsed().as_nanos() as f64 / f64::from(n);
    println!("{name:<28} {ns:8.1} ns per comparison");
}

fn main() {
    let root = Collator::root();
    time("root, ASCII", &root, "treatment_B", "treatment_a");
    time("root, accented", &root, "\u{e9}tat_1", "etat_2");
    time(
        "root, Cyrillic",
        &root,
        "\u{43a}\u{43e}\u{43d}\u{442}",
        "\u{43a}\u{43e}\u{442}",
    );
    let sv = Collator::new("sv").unwrap_or_default();
    time("sv, ASCII", &sv, "treatment_B", "treatment_a");
    time("sv, accented", &sv, "\u{e5}sa", "\u{e4}ng");
    time("C", &Collator::c(), "treatment_B", "treatment_a");
}
