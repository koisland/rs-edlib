# rs-edlib
Rust implementation of [edlib](https://github.com/Martinsos/edlib/tree/master).

### Usage

```bash
cargo add --git https://github.com/koisland/rs-edlib
```

```rust
use rs_edlib::{Alignment, AlignConfig};

let query: &str = "hello";
let target: &str = "world!";
let align_res = Alignment::run(
    AlignConfig::default(),
    query,
    target
).unwrap();
println!("edit_distance('{query}', '{target}') = {:?}", align_res.edit_distance)
```

### Why?

1. For learning purposes. Interested after reading about it in [`Racon`](https://genome.cshlp.org/content/27/5/737.full.pdf).
2. Need to get better at reading C/C++ code.
3. [`edlib-rs`](https://github.com/jean-pierreBoth/edlib-rs/tree/master) cmake build script is broken and the crate is just a thin wrapper for the C++ lib.

Some thoughts on original during rewrite:
* Great documentation!
    * Honestly most readable C/C++ code I've seen so far.
    * If I had one complaint, would be nice if some variables like `k` were described. Had to go to original [Myers 1999](https://dl.acm.org/doi/pdf/10.1145/316542.316550) paper to understand.
* It would be nice if a small example/diagram was provided demonstrating the banded block approach. Supplementary doc lacking imo.
    * This could be solved with a unittest.
* Unit tests for original implementation.
* Pointer arithmetic makes things harder than they need to be. Indices would be easier to read.
* Separating functions into separate files rather than a single file.


### References
* https://github.com/Martinsos/edlib/tree/master
* For more alignment info:
    * https://curiouscoding.nl/posts/astarpa2/
    * [`astar-pairwise-aligner`](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner)
