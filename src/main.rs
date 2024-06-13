#![allow(dead_code)]
mod poa;
mod lcskgraphefficient;
mod bit_tree;
use std::{collections::HashMap, process::exit};
use poa::*;
use petgraph::visit::Topo;
use crate::lcskgraphefficient::{find_sequence_in_graph, better_find_kmer_matches, lcskpp_graph, anchoring_lcsk_path_for_threading};
use rand::{Rng, SeedableRng, rngs::StdRng};
use std::time::Instant;
use rust_htslib::{bam, bam::Read};
use std::env;
use petgraph::dot::Dot;

const KMER_SIZE: usize = 10;
const SEQ_LEN: usize = 10000;
const BAND_SIZE: usize = 100;
const SUB_SECTION_LEN: usize = 1000;
const MIDDLE_SECTION_LEN: usize = 3162;

fn main() {
    arg_runner();
}

fn lcsk_test_pipeline(reads: Vec<String>, kmer_size: usize, band_size: usize) -> (usize, usize, usize, usize, usize, usize){
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

    let mut aligner = Aligner::new(2, -2, -2, &x, 0, 0, band_size as i32);
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
    let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = better_find_kmer_matches(&y, &all_sequences, &all_paths, kmer_size);
    //println!("LCSKgraph");
    //println!("{:?}", kmer_pos_vec);
    let (lcsk_path, lcsk_path_unconverted, _k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, all_paths.len(), kmer_size, kmer_graph_path, &topo_indices);
    println!("time for lcsk++ {:?}", now.elapsed());
    println!("lcsk++ path length {}", lcsk_path.len());
    // find the anchors here and display to test
    let anchors = anchoring_lcsk_path_for_threading(&lcsk_path_unconverted, &lcsk_path, 2, output_graph, 2);
    println!("{:?}", anchors);
    // get start and end from anchors and do poa for each section
    
    //let output_graph = aligner.graph();
    println!("{:?}", Dot::new(&output_graph.map(|_, n| (*n) as char, |_, e| *e)));
    lcsk_poa_score = aligner.semiglobal_banded(&y, &lcsk_path, band_size).alignment().score as usize;
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

fn lcsk_only_start (seed: usize) {
    let mut string_vec = get_random_sequences_from_generator(SEQ_LEN, 3, seed);
    let x = string_vec[0].as_bytes().to_vec();
    let y = string_vec.pop().unwrap().as_bytes().to_vec();

    let mut aligner = Aligner::new(2, -2, -2, &x, 0, 0, BAND_SIZE as i32);
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
    // find the kmers in top section 
    //println!("Finding kmers");

    // find matches in top section (full horizontal query, subsection length verticle graph)
    let (top_path, top_score) = get_lcsk_path_from_section (&y, &all_sequences, &all_paths, KMER_SIZE, &topo_indices, SEQ_LEN, SUB_SECTION_LEN);

    // find matches in left section (full verticle graph, subsection length query)
    let (left_path, left_score) = get_lcsk_path_from_section (&y, &all_sequences, &all_paths, KMER_SIZE, &topo_indices, SUB_SECTION_LEN, SUB_SECTION_LEN);
    // find matches in middle section (middle section length graph, middle section length query)
    let (middle_path, middle_score) = get_lcsk_path_from_section (&y, &all_sequences, &all_paths, KMER_SIZE, &topo_indices, MIDDLE_SECTION_LEN, MIDDLE_SECTION_LEN);
    // get the scores from all three and use the one with highest score
    println!("top {} {:?}", top_score, top_path[0]);
    println!("left{} {:?}", left_score, left_path[0]);
    println!("middle{} {:?}", middle_score, middle_path[0]);

    // and only use the start of kmer
    //let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = better_find_kmer_matches(&y, &all_sequences, &all_paths, KMER_SIZE);
    //println!("LCSKgraph");
    //println!("{:?}", kmer_pos_vec);
    //let (lcsk_path, _k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, all_paths.len(), KMER_SIZE, kmer_graph_path, &topo_indices);
    //println!("new with multi window {:?}", lcsk_path[0]);
    let elapsed = now.elapsed();
    println!("new TIME TAKEN {:?}", elapsed);
}

fn get_lcsk_path_from_section (query: &[u8], graph_sequences: &Vec<Vec<u8>>, graph_ids: &Vec<Vec<usize>>, k: usize, topo_map: &Vec<usize>, required_query_len: usize, required_graph_len: usize) -> (Vec<(usize, usize)>, u32){
    // cut the query
    let new_query;
    if required_query_len < query.len() {
        new_query = query[0..required_query_len].to_vec();
    }
    else {
        new_query = query.to_vec();
    }
    // cut the sequences and paths
    let mut new_sequences = vec![];
    let mut new_graph_ids = vec![];
    for (index, sequence) in graph_sequences.iter().enumerate() {
        let new_graph_sequence;
        let new_graph_id;
        if required_graph_len < sequence.len() {
            new_graph_sequence = sequence[0..required_graph_len].to_vec();
            new_graph_id = graph_ids[index][0..required_graph_len].to_vec();
        }
        else {
            new_graph_sequence = sequence.to_vec();
            new_graph_id = graph_ids[index].to_vec();
        }
        new_sequences.push(new_graph_sequence);
        new_graph_ids.push(new_graph_id);
    }
    // get the results
    let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = better_find_kmer_matches(&new_query, &new_sequences, &new_graph_ids, k);

    let (lcsk_path,lcsk_path_unconverted, k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, new_sequences.len(), k, kmer_graph_path, &topo_map);

    return (lcsk_path, k_new_score);
}

fn arg_runner() {
    // get the arguments 1. num of threads 2. read bam
    // s for synthetic test p for pacbio test s for subread processing
    // kmer size band size 
    let mut synthetic_data = false;
    let mut benchmarking = false; // only valid for pacbio data
    let mut kmer_size: usize = 10;
    let mut band_size: usize = 200;
    let mut sequence_length: usize = 10000;
    let mut number_of_iter: usize = 0; // required for benchmarking
    let mut input_bam_path: String = "data/sample_pacbio.bam".to_string();
    let mut output_bam_path: String = "output/read_pacbio.bam".to_string();

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match &arg[..] {
            "-h" | "--help" => help(),
            "--version" => {
                println!("LCSKGRAPH++ 1.0.0");
            }
            "-s" => {
                println!("Synthetic data mode");
                synthetic_data = true;
            }
            "-p" => {
                println!("Pacbio data mode");
                synthetic_data = false;
                if let Some(in_string) = args.next() {
                    input_bam_path = in_string;
                } else {
                    panic!("No input bam file specified!");
                }
                if let Some(out_string) = args.next() {
                    output_bam_path = out_string;
                } else {
                    panic!("No output bam file specified!");
                }
            }
            "-t" => {
                println!("Benchmarking..");
                benchmarking = true;
                if let Some(i_string) = args.next() {
                    number_of_iter = match i_string.parse::<usize>() {
                        Ok(x) => {x},
                        Err(_) => {panic!("Invalid value for -i")},
                    }
                } else {
                    panic!("No value specified for parameter -i");
                }
            }
            "-k" => {
                if let Some(k_string) = args.next() {
                    kmer_size = match k_string.parse::<usize>() {
                        Ok(x) => {x},
                        Err(_) => {panic!("Invalid value for -k")},
                    }
                } else {
                    panic!("No value specified for parameter -k.");
                }
            }
            "-b" => {
                if let Some(b_string) = args.next() {
                    band_size = match b_string.parse::<usize>() {
                        Ok(x) => {x},
                        Err(_) => {panic!("Invalid value for -b")},
                    }
                } else {
                    panic!("No value specified for parameter -b");
                }
            }
            "-l" => {
                if let Some(l_string) = args.next() {
                    sequence_length = match l_string.parse::<usize>() {
                        Ok(x) => {x},
                        Err(_) => {panic!("Invalid value for -l")},
                    }
                } else {
                    panic!("No value specified for parameter -l");
                }
            }
            _ => {
                if arg.starts_with('-') {
                    println!("Unknown argument {}", arg);
                } else {
                    println!("Unknown positional argument {}", arg);
                }
            }
        }
    }
    // check for synthetic and run
    if synthetic_data == true {
        if number_of_iter == 0 {
            panic!("Number of iterations can't be 0 for synthetic data benchmarking!");
        }
        if benchmarking == false {
            panic!("Synthetic data can only be run in benchmarking mode!");
        }
        println!("Benchmarking Synthetic data for {} iterations, sequences of length: {} using k: {} band size: {}", number_of_iter, sequence_length, kmer_size, band_size);
        run_synthetic_data_benchmark(kmer_size, sequence_length, number_of_iter, band_size);
    }
    // check for pacbio and run
    if synthetic_data == false {
        if benchmarking == true {
            if number_of_iter == 0 {
                panic!("Number of iterations can't be 0 for pacbio data benchmarking!");
            }
            println!("Benchmarking Pacbio data {} for {} iterations, using k: {} band size: {}", input_bam_path, number_of_iter, kmer_size, band_size);
            run_pacbio_data_benchmark(kmer_size, number_of_iter, band_size, input_bam_path);
        }
        else {
            make_read_file_from_subread_bam(kmer_size, band_size, input_bam_path, output_bam_path);
            println!("SS");
        }
    }
}

fn help() {
    println!("Usage: lcskgraph [OPTIONS]");
    println!("Options:\n -h, --help\tPrint help");
    println!(" --version\tPrint version information");
    println!(" -s \t\tSynthetic data, no output");
    println!(" -p <IN> <OUT>\tPacbio data, input subread bam file path, output bam file path");
    println!(" -t <N>\t\tBenchmarking, Number of iterations, default 10");
    println!(" -k <N>\t\tk for lcsk, default 10");
    println!(" -b <N>\t\tBand size for POA, default 200");
    println!(" -l <N>\t\tSequence length for synthetic data, default 10,000");
    println!(" -i <N>\t\tNumber of iterations for synthetic data or pacbio, default 0");
    // exit the program
    exit(0x0100);
}

fn make_read_file_from_subread_bam (_kmer_size: usize, _band_size: usize, _input_path: String, _output_path: String) {

}

fn run_pacbio_data_benchmark (kmer_size: usize, num_of_iter: usize, band_size: usize, input_path: String) {
    //println!("Processing Pacbio Data");
    let read_file_dir = input_path;
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
                    //println!("Start here");
                    current_set = record_set;
                    temp_read.push(String::from_utf8(record.seq().as_bytes()).unwrap());
                    //println!()
                }
                else if current_set == record_set {
                    //println!("Just adding read");
                    temp_read.push(String::from_utf8(record.seq().as_bytes()).unwrap());
                }
                else {
                    //println!("Read set complete onto the next");
                    current_set = record_set;
                    read_set.push(temp_read);
                    temp_read = vec![String::from_utf8(record.seq().as_bytes()).unwrap()];
                    if read_set.len() >= num_of_iter {
                        println!("Got {} sets exiting!", num_of_iter);
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
        //println!("{:?}", temp_reads.len());
        new_read_set.push(temp_reads);

    }
    let mut lcsk_stuff_sum = (0, 0, 0); // score time memory
    let mut poa_stuff_sum = (0, 0, 0);
    for (index, reads) in new_read_set.iter().enumerate() {
        println!("Progress {:.2}%", ((index * 100) as f32 / num_of_iter as f32));
        let string_vec = reads.clone();
        let results = lcsk_test_pipeline(string_vec, kmer_size, band_size);
        // print current seed results
        println!("Read number {}\nNormal poa \tScore: {} \tTime: {}meus \tMemory_usage: {}KB\nLcsk poa \tScore: {} \tTime: {}meus \tMemory usage: {}KB", index, results.1, results.3, results.5, results.0, results.2, results.4);
        lcsk_stuff_sum = (lcsk_stuff_sum.0 + results.0, lcsk_stuff_sum.1 + results.2, lcsk_stuff_sum.2 + results.4);
        poa_stuff_sum = (poa_stuff_sum.0 + results.1, poa_stuff_sum.1 + results.3, poa_stuff_sum.2 + results.5);
    }
    println!("=======================\nSummary Average of {} runs \nNormal poa\n\tScore: {}\n\tTime: {}meus\n\tMemory_usage: {}KB\nLcsk poa\n\tScore: {}\n\tTime: {}meus\n\tMemory usage: {}KB\n=======================\n", num_of_iter, poa_stuff_sum.0 / num_of_iter, poa_stuff_sum.1 / num_of_iter, poa_stuff_sum.2/ num_of_iter, lcsk_stuff_sum.0/ num_of_iter, lcsk_stuff_sum.1/ num_of_iter, lcsk_stuff_sum.2/ num_of_iter);
    //io::stdin().read_line(&mut String::new()).unwrap();
}

fn run_synthetic_data_benchmark (kmer_size: usize, sequence_length: usize, num_of_iter: usize, band_size: usize) {
    println!("Processing Synthetic Data");
    let mut lcsk_stuff_sum = (0, 0, 0); // score time memory
    let mut poa_stuff_sum = (0, 0, 0);
    for seed in 0..num_of_iter {
        println!("Progress {:.2}%", ((seed * 100) as f32 / num_of_iter as f32));
        let string_vec = get_random_sequences_from_generator(sequence_length, 3, seed);
        let results = lcsk_test_pipeline(string_vec, kmer_size, band_size);
        // print current seed results
        println!("Seed {}\nNormal poa \tScore: {} \tTime: {}meus \tMemory_usage: {}KB\nLcsk poa \tScore: {} \tTime: {}meus \tMemory usage: {}KB", seed, results.1, results.3, results.5, results.0, results.2, results.4);
        lcsk_stuff_sum = (lcsk_stuff_sum.0 + results.0, lcsk_stuff_sum.1 + results.2, lcsk_stuff_sum.2 + results.4);
        poa_stuff_sum = (poa_stuff_sum.0 + results.1, poa_stuff_sum.1 + results.3, poa_stuff_sum.2 + results.5);
    }
    println!("=======================\nSummary Average of {} runs \nNormal poa\n\tScore: {}\n\tTime: {}meus\n\tMemory_usage: {}KB\nLcsk poa\n\tScore: {}\n\tTime: {}meus\n\tMemory usage: {}KB\n=======================\n", num_of_iter, poa_stuff_sum.0 / num_of_iter, poa_stuff_sum.1 / num_of_iter, poa_stuff_sum.2/ num_of_iter, lcsk_stuff_sum.0/ num_of_iter, lcsk_stuff_sum.1/ num_of_iter, lcsk_stuff_sum.2/ num_of_iter);
    //io::stdin().read_line(&mut String::new()).unwrap();
}

fn get_random_sequences_from_generator(sequence_length: usize, num_of_sequences: usize, seed: usize) -> Vec<String> {
    let mut rng = StdRng::seed_from_u64(seed as u64);
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
        //println!("{:?}", mutseq.iter().collect::<String>());
        //insert to vector
        randomvec.push(mutseq.iter().collect::<String>());
    }
    randomvec
}