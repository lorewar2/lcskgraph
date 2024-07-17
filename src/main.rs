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
use std::thread;
use std::fs::OpenOptions;
use std::io::prelude::*;

const KMER_SIZE: usize = 10;
const SEQ_LEN: usize = 10000;
const BAND_SIZE: usize = 10000;
const SUB_SECTION_LEN: usize = 1000;
const MIDDLE_SECTION_LEN: usize = 3162;

fn main() {
    arg_runner();
}

fn process_the_reads_get_consensus_and_save_in_fa (input_name: &String, input_reads: Vec<String>, output_fa: &String, kmer_size: usize, band_size: usize) {
    let mut lcsk_aligner = Aligner::new(2, -2, -2, &input_reads[0].as_bytes().to_vec(), 0, 0, band_size as i32);
    let mut all_paths: Vec<Vec<usize>> = vec![];
    let mut all_sequences: Vec<Vec<u8>> = vec![];
    // run lcskpoa with the strings 
    for index in 0..input_reads.len() {
        if index == 0 {
            continue;
        }
        let output_graph = lcsk_aligner.graph();
        // initialize topology order and stuff
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
        // find the previous read path in graph
        let mut error_index = 0;
        loop {
            let (error_occured, temp_path, temp_sequence) = find_sequence_in_graph (input_reads[index - 1].as_bytes().to_vec().clone(), output_graph, &topo_indices, &topo_map, error_index);
            if error_index > 10 {
                break;
            }
            if !error_occured {
                all_paths.push(temp_path);
                all_sequences.push(temp_sequence);
                break;
            }
            error_index += 1;   
        }
        // get the lcsk path
        let query = &input_reads[index].as_bytes().to_vec();
        let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = better_find_kmer_matches(&query, &all_sequences, &all_paths, kmer_size);
        let (lcsk_path, _lcsk_path_unconverted, _k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, all_paths.len(), kmer_size, kmer_graph_path, &topo_indices);
        // poa
        lcsk_aligner.custom_banded(query, &lcsk_path, band_size).add_to_graph();
    }
    // get consensus 
    let consensus_u8 = lcsk_aligner.consensus();
    // write it to fa file
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .create(true)
        .open(output_fa)
        .unwrap();
    if let Err(e) = writeln!(file, ">{}", input_name) {
        eprintln!("Couldn't write to file: {}", e);
    }
    if let Err(e) = writeln!(file, "{}", String::from_utf8(consensus_u8).unwrap()) {
        eprintln!("Couldn't write to file: {}", e);
    }
}

fn make_output_file_from_subread_bam (kmer_size: usize, band_size: usize, input_path: String, output_path: String) {
    let read_file_dir = input_path;
    // get all the data from bam file
    let mut bam = bam::Reader::from_path(&read_file_dir).unwrap();
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
                    
                    // process here
                    process_the_reads_get_consensus_and_save_in_fa (&current_set, temp_read, &output_path, kmer_size, band_size);
                    temp_read = vec![String::from_utf8(record.seq().as_bytes()).unwrap()];
                    current_set = record_set;
                }
            },
            Err(_) => {
                break;
            }
        }
    }
}

fn lcsk_test_pipeline(reads: Vec<String>, kmer_size: usize, band_size: usize, cut_limit: usize) -> ((usize, usize, usize), (usize, usize, usize), (usize, usize, usize)){
    //for evaluating varibles
    let normal_stat: (usize, usize, usize); // score time memory
    let lcsk_stat: (usize, usize, usize);
    let threaded_stat: (usize, usize, usize);
    let only_lcsk_time;
    let mut children = vec![];
    // Make the aligner using full poa and first 2 reads
    let mut string_vec = reads.clone();
    let x = string_vec[0].as_bytes().to_vec();
    let y = string_vec.pop().unwrap().as_bytes().to_vec();
    let mut original_aligner = Aligner::new(2, -2, -2, &x, 0, 0, band_size as i32);
    for index in 1..string_vec.len() {
        original_aligner.global(&string_vec[index].as_bytes().to_vec()).add_to_graph();
    }
    let output_graph = original_aligner.graph();
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
    // find the path in graph (find graph node ids)
    for sequence in string_vec.clone() {
        let mut error_index = 0;
        loop {
            let (error_occured, temp_path, temp_sequence) = find_sequence_in_graph (sequence.as_bytes().to_vec().clone(), output_graph, &topo_indices, &topo_map, error_index);
            if error_index > 10 {
                break;
            }
            if !error_occured {
                all_paths.push(temp_path);
                all_sequences.push(temp_sequence);
                break;
            }
            error_index += 1;
            
        }
    }
    // LCSK PATH
    let now = Instant::now();
    let (kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, kmer_graph_path) = better_find_kmer_matches(&y, &all_sequences, &all_paths, kmer_size);
    let (lcsk_path, lcsk_path_unconverted, _k_new_score) = lcskpp_graph(kmer_pos_vec, kmer_path_vec, kmers_previous_node_in_paths, all_paths.len(), kmer_size, kmer_graph_path, &topo_indices);
    only_lcsk_time = now.elapsed().as_micros() as usize;
    //println!("time for lcsk++ {:?}", only_lcsk_time);
    // NORMAL LCSK POA
    let mut lcsk_aligner = original_aligner.clone();
    let now = Instant::now();
    let score = lcsk_aligner.custom_banded(&y, &lcsk_path, band_size).alignment().score;
    let time = now.elapsed().as_micros() as usize;
    let memory = lcsk_aligner.poa.memory_usage as usize;
    lcsk_stat = (score as usize, time + only_lcsk_time, memory);
    // THREADED LCSK POA
    let now = Instant::now();
    let mut total_section_score = 0;
    let mut total_section_memory = 0;
    if lcsk_path.len() > 0 {
        // find the anchors and graph sections (TODO intergrate query section finding in this and sections lcsk path)
        let (anchors, section_graphs, _node_tracker, section_queries, section_lcsks) = anchoring_lcsk_path_for_threading(&lcsk_path_unconverted, &lcsk_path, 2, output_graph, cut_limit,  y.len(), topo_indices, &y);
        for anchor_index in 0..anchors.len() - 1 {
            let section_query = section_queries[anchor_index].clone();
            let section_lcsk = section_lcsks[anchor_index].clone();
            let section_graph = section_graphs[anchor_index].clone();
            //println!("section query len {} section graph len {}", section_query.len(), section_graph.node_count());
            children.push(thread::spawn(move || {
                let mut aligner = Aligner::empty(2, -2, -2, 0, 0, band_size as i32);
                let score = aligner.custom_banded_threaded(&section_query, &section_lcsk, band_size, section_graph).alignment().score;
                let memory = aligner.poa.memory_usage as usize;
                (score, memory)
            }));
        }
        // get tall the section scores and add them up
        for _child_index in 0..children.len() {
            let result = children.pop().unwrap().join().unwrap();
            //println!("section child {} score {} memory {}", child_index, result.0, result.1);
            total_section_score += result.0;
            total_section_memory += result.1;
        }
        //println!("total section score {}", total_section_score);
    }
    let elapsed = now.elapsed();
    threaded_stat = (total_section_score as usize, elapsed.as_micros() as usize + only_lcsk_time, total_section_memory as usize);
    // NORMAL POA
    let now = Instant::now();
    let score = original_aligner.semiglobal(&y).alignment().score as usize;
    let elapsed = now.elapsed();
    let memory = original_aligner.poa.memory_usage as usize;
    let time = elapsed.as_micros() as usize;
    normal_stat = (score, memory, time);
    return (normal_stat, lcsk_stat, threaded_stat);
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
    let mut cut_limit: usize = 1000;

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
            "-c" => {
                if let Some(c_string) = args.next() {
                    cut_limit = match c_string.parse::<usize>() {
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
        run_synthetic_data_benchmark(kmer_size, sequence_length, number_of_iter, band_size, cut_limit);
    }
    // check for pacbio and run
    if synthetic_data == false {
        if benchmarking == true {
            if number_of_iter == 0 {
                panic!("Number of iterations can't be 0 for pacbio data benchmarking!");
            }
            println!("Benchmarking Pacbio data {} for {} iterations, using k: {} band size: {}", input_bam_path, number_of_iter, kmer_size, band_size);
            run_pacbio_data_benchmark(kmer_size, number_of_iter, band_size, input_bam_path, cut_limit);
        }
        else {
            make_output_file_from_subread_bam(kmer_size, band_size, input_bam_path, output_bam_path);
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
    println!(" -c <N>\t\t Cut limit for sectioning the graph and query, off by default");
    // exit the program
    exit(0x0100);
}

fn run_pacbio_data_benchmark (kmer_size: usize, num_of_iter: usize, band_size: usize, input_path: String, cut_limit: usize) {
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
    let mut normal_sum = (0, 0, 0);
    let mut threaded_sum = (0, 0, 0);
    let mut lcsk_sum = (0, 0, 0);
    for (index, reads) in new_read_set.iter().enumerate() {
        println!("Progress {:.2}%", ((index * 100) as f32 / num_of_iter as f32));
        let string_vec = reads.clone();
        let (normal, lcsk, threaded) = lcsk_test_pipeline(string_vec, kmer_size, band_size, cut_limit);
        // print current seed results
        normal_sum = (normal_sum.0 + normal.0, normal_sum.1 + normal.1, normal_sum.2 + normal.2);
        lcsk_sum = (lcsk_sum.0 + lcsk.0, lcsk_sum.1 + lcsk.1, lcsk_sum.2 + lcsk.2);
        threaded_sum = (threaded_sum.0 + threaded.0, threaded_sum.1 + threaded.1, threaded_sum.2 + threaded.2);
        println!("Read Number {}", index + 1);
        println!("Normal poa \t\tScore: {} \tTime: {}meus \tMemory_usage: {}KB", normal.0, normal.1, normal.2);
        println!("LCSK poa \t\tScore: {} \tTime: {}meus \tMemory_usage: {}KB", lcsk.0, lcsk.1, lcsk.2);
        println!("Threaded LCSK poa \tScore: {} \tTime: {}meus \tMemory_usage: {}KB", threaded.0, threaded.1, threaded.2);
    }
    println!("=======================\nSummary Average of {} runs", num_of_iter);
    println!("Normal poa \t\tScore: {} \tTime: {}meus \tMemory_usage: {}KB", normal_sum.0 / num_of_iter, normal_sum.1 / num_of_iter, normal_sum.2 / num_of_iter);
    println!("LCSK poa \t\tScore: {} \tTime: {}meus \tMemory_usage: {}KB", lcsk_sum.0 / num_of_iter, lcsk_sum.1 / num_of_iter, lcsk_sum.2 / num_of_iter);
    println!("Threaded LCSK poa \tScore: {} \tTime: {}meus \tMemory_usage: {}KB", threaded_sum.0 / num_of_iter, threaded_sum.1 / num_of_iter, threaded_sum.2 / num_of_iter);
}

fn run_synthetic_data_benchmark (kmer_size: usize, sequence_length: usize, num_of_iter: usize, band_size: usize, cut_limit: usize) {
    println!("Processing Synthetic Data");
    let mut normal_sum = (0, 0, 0);
    let mut threaded_sum = (0, 0, 0);
    let mut lcsk_sum = (0, 0, 0);
    for seed in 0..num_of_iter {
        println!("Progress {:.2}%", ((seed * 100) as f32 / num_of_iter as f32));
        let string_vec = get_random_sequences_from_generator(sequence_length, 3, seed);
        let (normal, lcsk, threaded) = lcsk_test_pipeline(string_vec, kmer_size, band_size, cut_limit);
        // print current seed results
        normal_sum = (normal_sum.0 + normal.0, normal_sum.1 + normal.1, normal_sum.2 + normal.2);
        lcsk_sum = (lcsk_sum.0 + lcsk.0, lcsk_sum.1 + lcsk.1, lcsk_sum.2 + lcsk.2);
        threaded_sum = (threaded_sum.0 + threaded.0, threaded_sum.1 + threaded.1, threaded_sum.2 + threaded.2);
        println!("Seed {}", seed);
        println!("Normal poa \t\tScore: {} \tTime: {}meus \tMemory_usage: {}KB", normal.0, normal.1, normal.2);
        println!("LCSK poa \t\tScore: {} \tTime: {}meus \tMemory_usage: {}KB", lcsk.0, lcsk.1, lcsk.2);
        println!("Threaded LCSK poa \tScore: {} \tTime: {}meus \tMemory_usage: {}KB", threaded.0, threaded.1, threaded.2);
    }
    println!("=======================\nSummary Average of {} runs", num_of_iter);
    println!("Normal poa \t\tScore: {} \tTime: {}meus \tMemory_usage: {}KB", normal_sum.0 / num_of_iter, normal_sum.1 / num_of_iter, normal_sum.2 / num_of_iter);
    println!("LCSK poa \t\tScore: {} \tTime: {}meus \tMemory_usage: {}KB", lcsk_sum.0 / num_of_iter, lcsk_sum.1 / num_of_iter, lcsk_sum.2 / num_of_iter);
    println!("Threaded LCSK poa \tScore: {} \tTime: {}meus \tMemory_usage: {}KB", threaded_sum.0 / num_of_iter, threaded_sum.1 / num_of_iter, threaded_sum.2 / num_of_iter);
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