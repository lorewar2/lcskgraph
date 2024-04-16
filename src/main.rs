mod lcskgraphdp;
mod lcskgraphefficient;
mod bit_tree;
use lcskgraphdp::*;
use petgraph::dot::Dot;
use petgraph::graph::{NodeIndex};

use crate::lcskgraphefficient::{simple_dfs_all_paths, find_kmer_matches};
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
    let mut all_paths: Vec<Vec<usize>> = vec![];
    let mut all_sequences: Vec<Vec<u8>> = vec![];
    simple_dfs_all_paths(output_graph, 0, vec![], vec![], &mut all_paths, &mut all_sequences);
    for path_index in 0..all_paths.len() {
        for node_index in 0..all_paths[path_index].len() {
            print!("{} {}, ", all_paths[path_index][node_index], all_sequences[path_index][node_index] as char);
        }
        println!();
    }
    let query = [65, 65, 65, 65];
    let result = find_kmer_matches(&query, &all_sequences, &all_paths, 3);
    for re in result {
        println!("{} {} {:?}", re.0, re.1, re.2);
    }
}