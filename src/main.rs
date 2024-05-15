mod poa;
mod poa_lcsk_banded;
mod lcskgraphefficient;
mod lcskgraphdp;
mod bit_tree;
use std::collections::HashMap;
use poa::*;
use petgraph::{dot::Dot, Direction::Incoming};
use petgraph::visit::Topo;
use crate::lcskgraphefficient::{divide_poa_graph_get_paths, find_sequence_in_graph, dfs_get_sequence_paths, better_find_kmer_matches, find_kmer_matches, find_kmer_matches_for_divided, lcskpp_graph, lcskpp_graph_for_divided, simple_dfs_all_paths};
use lcskgraphdp::Aligner as aligner2;
use poa_lcsk_banded::Aligner as aligner_banded;
use petgraph::graph::NodeIndex;
use rand::{Rng, SeedableRng, rngs::StdRng};

const CUT_THRESHOLD: usize = 5; //cut when number of nodes exceed this threshold
const KMER: usize = 10;
fn main() {
    let seed = 0;
    {
        println!("seed {}", seed);
        let mut string_vec = get_random_sequences_from_generator(50, 3, seed);
        let x = string_vec[0].as_bytes().to_vec();
        let y = string_vec.pop().unwrap().as_bytes().to_vec();
        
        let mut aligner = Aligner::new(2, -2, -2, &x, 0, 0, 1);
        for index in 1..string_vec.len() {
            aligner.global(&string_vec[index].as_bytes().to_vec()).add_to_graph();
        }
        
        let output_graph = aligner.graph();
        println!("{:?}", Dot::new(&output_graph.map(|_, n| (*n) as char, |_, e| *e)));
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
        //println!("Finding graph IDs");
        //dfs_get_sequence_paths(0,  string_vec.clone(), output_graph, topo_indices[0], vec![], vec![], &mut all_paths, &mut all_sequences, &topo_map);
        for sequence in string_vec.clone() {
            println!("{:?}", sequence);
            let mut error_index = 0;
            loop {
                let (error_occured, temp_path, temp_sequence) = find_sequence_in_graph (sequence.as_bytes().to_vec().clone(), output_graph, &topo_indices, &topo_map, error_index);
                if error_index > 10 {
                    println!("WHAT {} {:?}", sequence, temp_path);
                    break;
                }
                if !error_occured {
                    println!("{:?}", temp_path);
                    all_paths.push(temp_path);
                    all_sequences.push(temp_sequence);
                    break;
                }
                error_index += 1;
                
            }
        }
        //println!("Finding kmers");
        let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = better_find_kmer_matches(&y, &all_sequences, &all_paths, KMER);
        //println!("LCSKgraph");
        println!("{:?}", kmer_pos_vec);
        let (lcsk_path, k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, all_paths.len(), KMER, kmer_graph_path, &topo_indices);

        let mut aligner = aligner_banded::new(2, -2, -2, &x, 0, 0, 1);
        for index in 1..string_vec.len() {
            aligner.global(&string_vec[index].as_bytes().to_vec(), &lcsk_path, 10).add_to_graph();
        }
        //println!("{} {}", all_paths.len(), k_score);
        /*//println!("Getting paths by dividing..");
        //let (all_all_paths, all_all_sequences, max_paths) = divide_poa_graph_get_paths (output_graph, &topo_indices, 2, CUT_THRESHOLD, &topo_map);
        //let (kmer_pos_vec, kmers_plus_k, kmer_path_vec, kmers_previous_node_in_paths) = find_kmer_matches_for_divided(&y, &all_all_sequences, &all_all_paths, KMER);
        
        //println!("{:?}", kmers_plus_k);
        //let k_score = lcskpp_graph_for_divided(kmer_pos_vec, kmers_plus_k, kmer_path_vec, kmers_previous_node_in_paths, max_paths, KMER);
        //let (kmer_pos_vec, kmers_plus_k, kmer_path_vec, kmers_previous_node_in_paths) = find_kmer_matches(&y, &all_sequences, &all_paths, KMER);
        let mut all_paths: Vec<Vec<usize>> = vec![];
        let mut all_sequences: Vec<Vec<u8>> = vec![];
        // have to use all the nodes with no incoming nodes as starts... what a pain!!
        let mut start_nodes = vec![];
        for topo_index in &topo_indices {
            if output_graph.neighbors_directed(NodeIndex::new(*topo_index), Incoming).count() == 0 {
                start_nodes.push(topo_index);
            }
            else {
                break;
            }
        }
        for start_node in start_nodes {
            simple_dfs_all_paths(output_graph, *start_node, vec![], vec![], &mut all_paths, &mut all_sequences, &topo_map);
        }
        
        //println!("{}", all_paths.len());
        let (kmer_pos_vec, kmers_plus_k, kmer_path_vec, kmers_previous_node_in_paths, kmer_path_index) = better_find_kmer_matches(&y, &all_sequences, &all_paths, KMER);
        println!("{:?}", kmer_pos_vec);
        let k_score = lcskpp_graph( kmer_pos_vec,  kmer_path_vec, kmers_previous_node_in_paths, all_paths.len(), KMER, kmer_path_index);
        println!("all paths {:?}", all_paths);
        for path in all_paths {
            for node_id in path {
                print!("{} ", topo_indices[node_id]);
            }
            println!("");
        }
        let k_old_score = k_score;
        // test fulldplcsk++ 
        //let mut aligner2 = aligner2::new(0, 0, 0, &x);
        //aligner2.global(&y, KMER);
        //let dp_score = aligner2.traceback.get_score();
        //println!("efficient_score: {} dp_score: {}", k_score, dp_score);
        if k_new_score != k_old_score {
            println!("{} {}", k_new_score, k_old_score);
            mismatch_count += 1;
        }
        else {
            println!("Matched");
            match_count += 1;
        }
        assert!(k_new_score == k_old_score);
    }
    println!("{} {}", match_count, mismatch_count);*/
    }
}

fn get_random_sequences_from_generator(sequence_length: usize, num_of_sequences: usize, seed: u64) -> Vec<String> {
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
