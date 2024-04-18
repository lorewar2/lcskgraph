use crate::bit_tree::MaxBitTree;
use fxhash::FxHasher;
use std::cmp::{max};
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use petgraph::graph::{self, NodeIndex};
use petgraph::{Directed, Graph};

pub type POAGraph = Graph<u8, i32, Directed, usize>;
pub type HashMapFx<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

pub fn lcskpp_graph(kmer_pos_vec: Vec<(u32, u32)>, kmer_path_vec: Vec<Vec<usize>>, paths: &Vec<Vec<usize>>, k: usize) {
    // return nothing if empty
    if kmer_pos_vec.is_empty() {
        return;
    }

    let k = k as u32;

    let mut events: Vec<(u32, u32, u32, u32, Vec<usize>)> = Vec::new();
    let mut max_ns = vec![0; paths.len()];
    // generate the required events
    for (idx, &(x, y)) in kmer_pos_vec.iter().enumerate() {
        // in path find y and add k to the index and that is the value, get one path just use the first one
        let mut y_plusk_index = 0;
        
        if let Some(index_path) = paths[kmer_path_vec[idx][0]].iter().position(|r| r == &(y as usize)) {
            y_plusk_index = paths[kmer_path_vec[idx][0]][index_path + k as usize] as u32;
        }
        else {
            print!("SERIOUS ERROR Y IN PATH NOT FOUND");
        }
        // add the previous one as well
        events.push((x, y, y_plusk_index, (idx + kmer_pos_vec.len()) as u32, kmer_path_vec[idx].clone()));
        events.push((x + k, y, y_plusk_index, idx as u32, kmer_path_vec[idx].clone()));
        for path in kmer_path_vec[idx].clone() {
            max_ns[path] = max(max_ns[path], x + k);
            max_ns[path] = max(max_ns[path], y_plusk_index);
        }
    }
    // sorting is okay with topologically converted indices
    events.sort_unstable();
    for event in &events {
        println!("{} {} {} {} {:?}", event.0, event.1, event.2, event.3, event.4);
    }
    let mut max_bit_tree_path = vec![];
    // generate empty fenwick trees
    for (index, n) in max_ns.iter().enumerate() {
        let max_col_dp: MaxBitTree<(u32, u32)> = MaxBitTree::new(*n as usize);
        max_bit_tree_path.push(max_col_dp);
    }
    let mut dp: Vec<(u32, i32)> = Vec::with_capacity(events.len());
    let mut best_dp = (k, 0);

    dp.resize(events.len(), (0, 0));
    for ev in events {
        // p is the match index
        let p = (ev.2 % kmer_pos_vec.len() as u32) as usize;
        // the the graph indices in this case is j
        let j = ev.1;
        // is start if higher than this
        let is_start = ev.2 >= (kmer_pos_vec.len() as u32);

        if is_start {
            dp[p] = (k, -1);
            // go through the paths available in this event, and get the max corrosponding value and pos
            for path in ev.4 {
                let (temp_value, temp_position) = max_bit_tree_path[path].get(j as usize);
                if (temp_value > dp[p].0) && (temp_value > 0) {
                    dp[p] = (k + temp_value, temp_position as i32);
                    best_dp = max(best_dp, (dp[p].0, p as i32));
                }
            } // done until here
        } else {
            // See if this kmer continues a different kmer
            // this has to be checked by path if enough nodes are available from end ev.1
            if ev.0 > k && ev.1 > k {
                // check the diagonally prev slot for a kmer
                // need the events path and find the k - 1 from path to find p start of continueing so event need the path index as well
                if let Ok(cont_idx) = kmer_pos_vec.binary_search(&(ev.0 - k - 1, ev.1 - k - 1)) {
                    let prev_score = dp[cont_idx].0;
                    let candidate = (prev_score + 1, cont_idx as i32);
                    dp[p] = max(dp[p], candidate);
                    best_dp = max(best_dp, (dp[p].0, p as i32));
                }
            }
            // add this end max col dp find which path fenwick tree and set it
            for path in ev.4 {
                max_bit_tree_path[path].set(ev.1 as usize, (dp[p].0, p as u32));
            }
        }
    }

    let mut traceback = Vec::new();
    let (best_score, mut prev_match) = best_dp;
    while prev_match >= 0 {
        traceback.push(prev_match as usize);
        prev_match = dp[prev_match as usize].1;
    }
    traceback.reverse();
    //path :traceback;
    //score: best_score,
    //dp_vector: dp, */
}

pub fn convert_topological_indices_to_ascending_indices (topo_indices: &Vec<usize>, paths: &Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let mut converted_paths: Vec<Vec<usize>> = vec![];
    for path in paths {
        let mut temp_converted_path = vec![];
        let mut path_index = 0;
        let mut ascending_index = 0;
        for topo_index in topo_indices {
            println!("{} {}", topo_index, path[path_index]);
            if *topo_index == path[path_index] {
                temp_converted_path.push(ascending_index);
                path_index += 1;
            }
            ascending_index += 1;
        }
        converted_paths.push(temp_converted_path);
    }
    converted_paths
}

pub fn find_kmer_matches(query: &[u8], graph_sequences: &Vec<Vec<u8>>, graph_ids: &Vec<Vec<usize>>, k: usize) -> (Vec<(u32, u32)>, Vec<Vec<usize>>) {
    // hash the query
    let set = hash_kmers(query, k);
    // go through the paths and get the indices of path and make a list with query index, graph index, paths
    let mut kmers_result_vec: Vec<(u32, u32)> = vec![];
    let mut kmers_paths: Vec<Vec<usize>> = vec![];
    for (index, seq) in graph_sequences.iter().enumerate() {
        let matches = find_kmer_matches_seq1_hashed(&set, seq, k);
        // go through the matches and see if they are in the already made list
        for a_match in matches {
            // first get the graph index for the match
            let graph_index = graph_ids[index][a_match.1 as usize] as u32;
            // if it is in the result vec, add path if indices match and path different
            if let Some(index_result) = kmers_result_vec.iter().position(|r| (r.0, r.1) == (a_match.0, graph_index)) {
                kmers_paths[index_result].push(index);
            }
            // not in the result vec just add
            else {
                kmers_result_vec.push((a_match.0, graph_index));
                kmers_paths.push(vec![index]);
            }
        }
    }
    (kmers_result_vec, kmers_paths)
}

pub fn simple_dfs_all_paths (
    graph: &POAGraph, 
    start: usize,
    current_path: Vec<usize>,
    current_sequence: Vec<u8>,
    all_paths: &mut Vec<Vec<usize>>,
    all_sequences: &mut Vec<Vec<u8>>
) {
    let mut path = current_path.clone();
    let mut sequence = current_sequence.clone();
    path.push(start);
    sequence.push(graph.raw_nodes()[start].weight);
    // if no neighbours, last node reached
    if graph.neighbors(NodeIndex::new(start)).count() == 0 {
        all_paths.push(path.clone());
        all_sequences.push(sequence.clone());
    } else {
        // go through neighbours recursively
        for neighbor in graph.neighbors(NodeIndex::new(start)) {
            simple_dfs_all_paths(graph, neighbor.index(), path.clone(), sequence.clone(), all_paths, all_sequences);
        }
    }
}

pub fn hash_kmers(seq: &[u8], k: usize) -> HashMapFx<&[u8], Vec<u32>> {
    let slc = seq;
    let mut set: HashMapFx<&[u8], Vec<u32>> = HashMapFx::default();
    for i in 0..(slc.len() + 1).saturating_sub(k) {
        set.entry(&slc[i..i + k])
            .or_insert_with(Vec::new)
            .push(i as u32);
    }
    set
}

pub fn find_kmer_matches_seq1_hashed(
    seq1_set: &HashMapFx<&[u8], Vec<u32>>,
    seq2: &[u8],
    k: usize,
) -> Vec<(u32, u32)> {
    let mut matches = Vec::new();

    for i in 0..(seq2.len() + 1).saturating_sub(k) {
        let slc = &seq2[i..i + k];
        if let Some(matches1) = seq1_set.get(slc) {
            for pos1 in matches1 {
                matches.push((*pos1, i as u32));
            }
        }
    }
    matches.sort_unstable();
    matches
}