mod poa;
mod lcskgraphefficient;
mod bit_tree;
use poa::*;
use petgraph::dot::Dot;
use petgraph::visit::Topo;
use crate::lcskgraphefficient::{simple_dfs_all_paths, find_kmer_matches, lcskpp_graph, convert_topological_indices_to_ascending_indices};

fn main() {
    // test run the lcsk++ incomplete code
    let x = b"CTATAGAGTA".to_vec();
    let y = b"ATTATG".to_vec();
    let aligner = Aligner::new(2, -2, -2, &x, 0, 0, 1);
    let output_graph = aligner.graph();
    println!("{:?}", Dot::new(&output_graph.map(|_, n| (*n) as char, |_, e| *e)));
    let mut all_paths: Vec<Vec<usize>> = vec![];
    let mut all_sequences: Vec<Vec<u8>> = vec![];
    simple_dfs_all_paths(output_graph, 0, vec![], vec![], &mut all_paths, &mut all_sequences);
    let (kmer_pos_vec, kmer_path_vec) = find_kmer_matches(&y, &all_sequences, &all_paths, 2);
    for re in &kmer_pos_vec {
        println!("{} {}", re.0, re.1);
    }
    // get topology ordering
    let mut topo = Topo::new(&output_graph);
    // go through the nodes topologically
    let mut topo_indices = vec![];
    while let Some(node) = topo.next(&output_graph) {
        topo_indices.push(node.index());
    }
    let (converted_paths , converted_matches) = convert_topological_indices_to_ascending_indices(&topo_indices, &all_paths, &kmer_pos_vec);
    
    println!("{:?}", topo_indices);
    println!("{:?}", converted_paths);
    for path_index in 0..all_paths.len() {
        for node_index in 0..all_paths[path_index].len() {
            print!("{}=={} {}, ", all_paths[path_index][node_index], converted_paths[path_index][node_index], all_sequences[path_index][node_index] as char);
        }
        println!();
    }
    println!("{:?}", converted_matches);
    lcskpp_graph(converted_matches, kmer_path_vec, &converted_paths, 2);
}