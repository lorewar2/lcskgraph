mod poa;
mod lcskgraphefficient;
mod lcskgraphdp;
mod bit_tree;
use std::collections::HashMap;
use poa::*;
use petgraph::dot::Dot;
use petgraph::visit::Topo;
use crate::lcskgraphefficient::{simple_dfs_all_paths, find_kmer_matches, lcskpp_graph};
use lcskgraphdp::Aligner as aligner2;
use rand::{Rng, SeedableRng, rngs::StdRng};

fn main() {
    // test run the lcsk++ incomplete code
    //let x = b"ATAGTAAAATATATG".to_vec(); // test case 1 which was fixed
    //let y = b"ATTATG".to_vec();
    // test case 2
    //let y = b"ATTATAAAG".to_vec();
    //let x = b"ATAGTAAAATATATG".to_vec();
    //let x = b"CTATAGAGTA".to_vec();
    //let y = b"ATTATG".to_vec();
    for seed in 0..100000 {
        let string_vec = get_random_sequences_from_generator(10, 2, seed);
        let x = string_vec[0].as_bytes().to_vec();
        let y = string_vec[1].as_bytes().to_vec();
        let aligner = Aligner::new(2, -2, -2, &x, 0, 0, 1);
        let output_graph = aligner.graph();
        //println!("{:?}", Dot::new(&output_graph.map(|_, n| (*n) as char, |_, e| *e)));
        let mut all_paths: Vec<Vec<usize>> = vec![];
        let mut all_sequences: Vec<Vec<u8>> = vec![];

        // get topology ordering
        let mut topo = Topo::new(&output_graph);
        // go through the nodes topologically // make a hashmap with node_index as key and incrementing indices as value
        let mut topo_indices = vec![];
        let mut topo_map = HashMap::new();
        let mut incrementing_index: usize = 0;
        while let Some(node) = topo.next(&output_graph) {
            topo_indices.push(node.index());
            topo_map.insert(node.index(), incrementing_index);
            incrementing_index += 1;
        }
        simple_dfs_all_paths(output_graph, 0, vec![], vec![], &mut all_paths, &mut all_sequences, &topo_map);
        let (kmer_pos_vec, kmers_plus_k, kmer_path_vec, kmers_previous_node_in_paths) = find_kmer_matches(&y, &all_sequences, &all_paths, 2);
        let k_score = lcskpp_graph(kmer_pos_vec, kmers_plus_k, kmer_path_vec, kmers_previous_node_in_paths, &all_paths, 2);
        // test fulldplcsk++ 
        let mut aligner2 = aligner2::new(0, 0, 0, &x);
        aligner2.global(&y, 2);
        let dp_score = aligner2.traceback.get_score();
        println!("efficient_score: {} dp_score: {}", k_score, dp_score);
        assert!(k_score == dp_score as u32);
    }
}

pub fn get_random_sequences_from_generator(sequence_length: usize, num_of_sequences: usize, seed: u64) -> Vec<String> {
    let mut rng = StdRng::seed_from_u64(seed);
    //vector to save all the sequences 
    let mut randomvec: Vec<String> = vec![];
    //generate the first sequence of random bases of length sequence_length
    let mut firstseq: Vec<char> = vec![];
    for _ in 0..sequence_length {
        firstseq.push(match rng.gen_range(0..4) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'X'
        });
    }
    //randomvec.push(firstseq.iter().collect::<String>());
    //loop for 10 
    for _ in 0..num_of_sequences {
        //clone the sequence
        let mut mutseq = firstseq.clone();
        //mutate the all the bases with 0.05 chance
        for i in 0..mutseq.len() {
            match rng.gen_range(0..20) {
                0 => {
                    mutseq[i] = match rng.gen_range(0..4) {
                        0 => 'A',
                        1 => 'C',
                        2 => 'G',
                        3 => 'T',
                        _ => 'X'
                    }
                },
                _ => {}
            }
        }
        //put indels at location with chance 0.1 
        for i in 0..mutseq.len() {
            let mean_value: f64 = 1.5; //2.0 before
            //get length of the indel geometric distributed mean value 1.5
            let indel_length: usize  = ((1.0 - rng.gen::<f64>()).ln() / (1.00 - (1.00 / mean_value) as f64).ln()).ceil() as usize;
            match rng.gen_range(0..20) {
                //insertion of elements
                0 => {
                    if i + indel_length < mutseq.len() {
                        for _ in 0..indel_length{
                            mutseq.insert(i + 1, mutseq[i]);
                        }
                    }
                },
                //deletion of elements
                1 => {
                    if i + indel_length < mutseq.len() {
                        for _ in 0..indel_length{
                            mutseq.remove(i);
                        }
                    }
                }
                _ => {}
            }
        }
        println!("{:?}", mutseq.iter().collect::<String>());
        //insert to vector
        randomvec.push(mutseq.iter().collect::<String>());
    }
    randomvec
}
