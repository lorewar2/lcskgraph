mod lcskgraphdp;
mod lcskgraphefficient;
mod bittree;
use lcskgraphdp::*;
use petgraph::dot::Dot;
fn main() {
    // test run the lcsk++ incomplete code
    let x = b"AAAAAAA".to_vec();
    let y = b"AABBBAA".to_vec();
    let z = b"AABCBAA".to_vec();
    let mut aligner = Aligner::new(2, -2, -2, &x);
    // z differs from x in 3 locations
    aligner.global(&y).add_to_graph();
    aligner.global(&z).add_to_graph();
    let output_graph = aligner.graph();
    println!("{:?}", Dot::new(&output_graph.map(|_, n| (*n) as char, |_, e| *e)));
}