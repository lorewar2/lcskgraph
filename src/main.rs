mod poa;
mod lcskgraphefficient;
mod lcskgraphdp;
mod bit_tree;
use std::{collections::HashMap};
use poa::*;
use petgraph::visit::Topo;
use crate::lcskgraphefficient::{find_sequence_in_graph, better_find_kmer_matches, lcskpp_graph};
use rand::{Rng, SeedableRng, rngs::StdRng};
use std::time::Instant;
use rust_htslib::{bam, bam::Read};
use std::env;

//const CUT_THRESHOLD: usize = 5; //cut when number of nodes exceed this threshold
const KMER: usize = 4;
const SEQ_LEN: usize = 100;
const NUM_OF_ITER: u64 = 10;
const BAND_SIZE: usize = 10;

fn main() {
    // get the arguments 1. num of threads 2. read bam
    // s for synthetic test p for pacbio test s for subread processing
    // kmer size band size 
    let mut kmer_size: usize = 10;
    let mut band_size: usize = 200;
    let mut sequence_length: usize = 10000;
    let mut bam_path: String = "".to_string();

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match &arg[..] {
            "-h" | "--help" => help(),
            "--version" => {
                println!("{} {}", "LCSKGRAPH++", "1.0.0");
            }
            "-v" | "--verbose" => {
                println!("Verbose mode is not supported yet.");
            }
            "-c" | "--config" => {
                if let Some(arg_config) = args.next() {
                    //config = arg_config;
                } else {
                    panic!("No value specified for parameter --config.");
                }
            }
            _ => {
                if arg.starts_with('-') {
                    println!("Unkown argument {}", arg);
                } else {
                    println!("Unkown positional argument {}", arg);
                }
            }
        }
    }
   //run_pacbio_data();
   //run_synthetic_data();
}

fn help() {
    println!("Usage: lcskgraph [OPTIONS]");
    println!("Options:\n -h, --help\tPrint help");
    println!("\n--version\tPrint version information");
    println!("\n-t <0,1,2>\tk for lcsk");
    println!("\n-k <N>\tk for lcsk");
    println!("\n-b <N>\tBand size for POA");
    println!("\n-l <N>\tSequence length for synthetic data");
    println!("\n-i <N>\tThe sequence length for synthetic data");
    println!("\n-f <PATH>\tSpecify the path for pacbio subread bam file");
}

fn make_read_file_from_subread_bam () {

}

fn run_pacbio_data() {
    println!("Processing Pacbio Data");
    let read_file_dir = "data/sample_pacbio.bam";
    // get data from bam file
    let mut bam = bam::Reader::from_path(&read_file_dir).unwrap();
    let mut read_set = vec![];
    let mut current_set = "".to_string();
    let mut temp_read = vec![];
    for record_option in bam.records().into_iter() {
        match record_option {                                                                                                       
            Ok(x) => {                                                                                                       
                let record = x;                                                                                                         
                let record_set = String::from_utf8(record.qname().to_vec()).unwrap().split("/").collect::<Vec<&str>>()[1].to_string();
                if current_set == "".to_string() {
                    println!("Start here");
                    current_set = record_set;
                    temp_read.push(String::from_utf8(record.seq().as_bytes()).unwrap());
                    println!()
                }
                else if current_set == record_set {

                    println!("Just adding read");
                    temp_read.push(String::from_utf8(record.seq().as_bytes()).unwrap());
                }
                else {
                    println!("Read set complete onto the next");
                    current_set = record_set;
                    read_set.push(temp_read);
                    temp_read = vec![String::from_utf8(record.seq().as_bytes()).unwrap()];
                    if read_set.len() >= 10 {
                        println!("Got 10 sets exiting!");
                        break;
                    }
                }
            },
            Err(_) => {
                break;
            }
        }
    }
    let mut new_read_set = vec![];
    // reverse the reads and stuff
    for reads in &read_set {
        let mut temp_reads = vec![];
        temp_reads.push(reads[0].clone());
        temp_reads.push(reads[2].clone());
        temp_reads.push(reads[4].clone());
        println!("{:?}", temp_reads.len());
        new_read_set.push(temp_reads);

    }
    for (index, reads) in new_read_set.iter().enumerate() {
        println!("index {}", index);
        let string_vec = reads.clone();
        lcsk_test_pipeline(string_vec);
    }
    let mut lcsk_stuff_sum = (0, 0, 0); // score time memory
    let mut poa_stuff_sum = (0, 0, 0);
    for (index, reads) in new_read_set.iter().enumerate() {
        let string_vec = reads.clone();
        let results = lcsk_test_pipeline(string_vec);
        // print current seed results
        println!("Read number {}\nNormal poa\n\tScore: {}\n\tTime: {}meus\n\tMemory_usage: {}KB\nLcsk poa\n\tScore: {}\n\tTime: {}meus\n\tMemory usage: {}KB\n", index, results.1, results.3, results.5, results.0, results.2, results.4);
        lcsk_stuff_sum = (lcsk_stuff_sum.0 + results.0, lcsk_stuff_sum.1 + results.2, lcsk_stuff_sum.2 + results.4);
        poa_stuff_sum = (poa_stuff_sum.0 + results.1, poa_stuff_sum.1 + results.3, poa_stuff_sum.2 + results.5);
    }
    println!("=======================\nAverage \nNormal poa\n\tScore: {}\n\tTime: {}meus\n\tMemory_usage: {}KB\nLcsk poa\n\tScore: {}\n\tTime: {}meus\n\tMemory usage: {}KB\n=======================\n", poa_stuff_sum.0 / NUM_OF_ITER as usize, poa_stuff_sum.1 / NUM_OF_ITER as usize, poa_stuff_sum.2/ NUM_OF_ITER as usize, lcsk_stuff_sum.0/ NUM_OF_ITER as usize, lcsk_stuff_sum.1/ NUM_OF_ITER as usize, lcsk_stuff_sum.2/ NUM_OF_ITER as usize);
    //io::stdin().read_line(&mut String::new()).unwrap();
}

fn run_synthetic_data() {
    println!("Processing Synthetic Data");
    let mut lcsk_stuff_sum = (0, 0, 0); // score time memory
    let mut poa_stuff_sum = (0, 0, 0);
    for seed in 0..NUM_OF_ITER {
        let string_vec = get_random_sequences_from_generator(SEQ_LEN, 3, seed);
        let results = lcsk_test_pipeline(string_vec);
        // print current seed results
        println!("Seed {}\nNormal poa\n\tScore: {}\n\tTime: {}meus\n\tMemory_usage: {}KB\nLcsk poa\n\tScore: {}\n\tTime: {}meus\n\tMemory usage: {}KB\n", seed, results.1, results.3, results.5, results.0, results.2, results.4);
        lcsk_stuff_sum = (lcsk_stuff_sum.0 + results.0, lcsk_stuff_sum.1 + results.2, lcsk_stuff_sum.2 + results.4);
        poa_stuff_sum = (poa_stuff_sum.0 + results.1, poa_stuff_sum.1 + results.3, poa_stuff_sum.2 + results.5);
    }
    println!("=======================\nAverage \nNormal poa\n\tScore: {}\n\tTime: {}meus\n\tMemory_usage: {}KB\nLcsk poa\n\tScore: {}\n\tTime: {}meus\n\tMemory usage: {}KB\n=======================\n", poa_stuff_sum.0 / NUM_OF_ITER as usize, poa_stuff_sum.1 / NUM_OF_ITER as usize, poa_stuff_sum.2/ NUM_OF_ITER as usize, lcsk_stuff_sum.0/ NUM_OF_ITER as usize, lcsk_stuff_sum.1/ NUM_OF_ITER as usize, lcsk_stuff_sum.2/ NUM_OF_ITER as usize);
    //io::stdin().read_line(&mut String::new()).unwrap();
}

fn lcsk_test_pipeline(reads: Vec<String>) -> (usize, usize, usize, usize, usize, usize){
    //for evaluating varibles
    let lcsk_poa_score;
    let normal_poa_score;
    let lcsk_poa_time;
    let normal_poa_time;
    let lcsk_poa_memory;
    let normal_poa_memory;

    let mut string_vec = reads.clone();
    let x = string_vec[0].as_bytes().to_vec();
    let y = string_vec.pop().unwrap().as_bytes().to_vec();

    let mut aligner = Aligner::new(2, -2, -2, &x, 0, 0, 1);
    for index in 1..string_vec.len() {
        aligner.global(&string_vec[index].as_bytes().to_vec()).add_to_graph();
    }
    
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
    

    //println!("Finding graph IDs");
    //dfs_get_sequence_paths(0,  string_vec.clone(), output_graph, topo_indices[0], vec![], vec![], &mut all_paths, &mut all_sequences, &topo_map);
    for sequence in string_vec.clone() {
        //println!("{:?}", sequence);
        let mut error_index = 0;
        loop {
            let (error_occured, temp_path, temp_sequence) = find_sequence_in_graph (sequence.as_bytes().to_vec().clone(), output_graph, &topo_indices, &topo_map, error_index);
            if error_index > 10 {
                //println!("WHAT {} {:?}", sequence, temp_path);
                break;
            }
            if !error_occured {
                //println!("{:?}", temp_path);
                all_paths.push(temp_path);
                all_sequences.push(temp_sequence);
                break;
            }
            error_index += 1;
            
        }
    }
    let now = Instant::now();
    //println!("Finding kmers");
    let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = better_find_kmer_matches(&y, &all_sequences, &all_paths, KMER);
    //println!("LCSKgraph");
    //println!("{:?}", kmer_pos_vec);
    let (lcsk_path, _k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, all_paths.len(), KMER, kmer_graph_path, &topo_indices);
    
    //let output_graph = aligner.graph();
    //println!("{:?}", Dot::new(&output_graph.map(|_, n| (*n) as char, |_, e| *e)));
    lcsk_poa_score = aligner.semiglobal_banded(&y, &lcsk_path, BAND_SIZE).alignment().score as usize;
    let elapsed = now.elapsed();
    lcsk_poa_memory = aligner.poa.memory_usage as usize;
    lcsk_poa_time = elapsed.as_micros() as usize;
    //println!("Elapsed: {:.2?}", elapsed);
    //println!("score {}", lcsk_poa_score);
   // second try   
    let mut aligner = Aligner::new(2, -2, -2, &x, 0, 0, 1);
    for index in 1..string_vec.len() {
        aligner.global(&string_vec[index].as_bytes().to_vec()).add_to_graph();
    }
    let now = Instant::now();
    normal_poa_score = aligner.semiglobal(&y).alignment().score as usize;
    let elapsed = now.elapsed();
    normal_poa_memory = aligner.poa.memory_usage as usize;
    //println!("score {}", normal_poa_score);
    normal_poa_time = elapsed.as_micros() as usize;
    //println!("Elapsed: {:.2?}", elapsed);
    return (lcsk_poa_score, normal_poa_score, lcsk_poa_time, normal_poa_time, lcsk_poa_memory, normal_poa_memory)
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