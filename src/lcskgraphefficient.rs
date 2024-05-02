use crate::bit_tree::MaxBitTree;
use fxhash::FxHasher;
use std::cmp::max;
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use petgraph::Outgoing;
use petgraph::graph::NodeIndex;
use petgraph::{Directed, Graph};

pub type POAGraph = Graph<u8, i32, Directed, usize>;
pub type HashMapFx<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

pub fn lcskpp_graph_for_divided (
    kmer_pos_vec: Vec<(u32, u32, usize)>,
    kmers_plus_k: Vec<u32>,
    kmer_path_vec: Vec<Vec<(usize, usize)>>,
    kmers_previous_node_in_paths: Vec<Vec<u32>>,
    num_paths: usize,
    k: usize
    ) {
    // return nothing if empty
    if kmer_pos_vec.is_empty() {
        return;
    }

    let k = k as u32;

    let mut events: Vec<(u32, u32, u32, u32, Vec<(usize, usize)>, Vec<u32>)> = Vec::new();
    let mut max_ns = vec![0; num_paths];
    // generate the required events
    for (idx, &(x, y, section)) in kmer_pos_vec.iter().enumerate() {
        let y_plus_k = kmers_plus_k[idx];
        // add the previous one as well
        // IF END SAVE y - K - 1 NODE IF START SAVE y - 1 NODE (-1 if not available)
        events.push((x, y, y_plus_k, (idx + kmer_pos_vec.len()) as u32, kmer_path_vec[idx].clone(), kmers_previous_node_in_paths[idx].clone()));
        events.push((x + k - 1, y, y_plus_k, idx as u32, kmer_path_vec[idx].clone(), kmers_previous_node_in_paths[idx].clone()));
        for path in kmer_path_vec[idx].clone() {
            max_ns[path.0] = max(max_ns[path.0], x + k - 1);
            max_ns[path.0] = max(max_ns[path.0], y_plus_k);
        }
    }
    // sorting is okay with topologically converted indices
    events.sort_unstable();
    for event in events {
        println!("{:?}", event);
    }
    
    /*let mut max_bit_tree_path = vec![];
    // generate empty fenwick trees
    for (_index, n) in max_ns.iter().enumerate() {
        let max_col_dp: MaxBitTree<(u32, u32)> = MaxBitTree::new(*n as usize);
        max_bit_tree_path.push(max_col_dp);
    }
    let mut dp: Vec<(u32, i32)> = Vec::with_capacity(events.len());
    let mut best_dp = (k, 0, 0); // score, coloumn, path
    // update trees only after the next query point is hit
    let mut tree_update_required_level = 0;
    let mut tree_update_vec: Vec<(u32, (u32, u32), u32, Vec<usize>)> = vec![(0, (0, 0), 0, vec![])]; //current x, value,  index, paths (trees)
    dp.resize(events.len(), (0, 0));
    for ev in events {
        if tree_update_required_level != ev.0 && tree_update_vec.len() != 0 {
            // update trees here
            for tree_info in &tree_update_vec {
                for path in &tree_info.3 {
                    max_bit_tree_path[*path].set(tree_info.2 as usize, tree_info.1);
                }
            }
            //tree_update_vec = vec![];
        }
        // p is the match index
        let p = (ev.3 % kmer_pos_vec.len() as u32) as usize;
        //println!("x:{} y_start:{} y_end:{} p:{} idx:{} paths:{:?}", ev.0, ev.1, ev.2, p, ev.3, ev.4);
        // is start if higher than this
        let is_start = ev.3 >= (kmer_pos_vec.len() as u32);
        if is_start {
            //print!("IS START \n");
            dp[p] = (k, -1);
            // go through the paths available in this event, and get the max corrosponding value and pos
            for path_index in 0..ev.4.len() {
                let path = ev.4[path_index];
                let prev_node = ev.5[path_index];
                if prev_node != u32::MAX {
                    //println!("prev node value = {}", prev_node);
                    let (temp_value, temp_position) = max_bit_tree_path[path].get(prev_node as usize);
                    //println!("temp value from fenwick tree {}", temp_value);
                    if (temp_value + k > dp[p].0) && (temp_value > 0) {
                        dp[p] = (k + temp_value, temp_position as i32);
                        best_dp = max(best_dp, (dp[p].0, p as i32, path));
                        //println!("best_dp {}", best_dp.0);
                    }
                }
            }
        } else {
            //print!("IS END \n");
            // See if this kmer continues a different kmer
            if ev.0 >= k {
                for path_index in 0..ev.4.len() {
                    let path = ev.4[path_index];
                    let prev_node = ev.5[path_index];
                    if prev_node != u32::MAX {
                        if let Ok(cont_idx) = kmer_pos_vec.binary_search(&(ev.0 - k, prev_node)) {
                            //println!("!!!!!!!!!!");
                            let prev_score = dp[cont_idx].0;
                            let candidate = (prev_score + 1, cont_idx as i32);
                            dp[p] = max(dp[p], candidate);
                            best_dp = max(best_dp, (dp[p].0, p as i32, path));
                            //println!("candidate location {} {}", ev.0 - k, prev_node);
                            //println!("cont value from candidate p {} score {}", cont_idx, candidate.0);
                            //println!("best_dp {}", best_dp.0);
                        }
                    }
                }
            }
            // set all trees which have this match as this // maybe update this in the next iteration to prevent query overlapping
            //  update required, current x, value,  index, paths (trees)
            tree_update_vec.push((ev.0, (dp[p].0, p as u32), ev.2, ev.4.clone()));
            tree_update_required_level = ev.0;
        }
    }
    let mut traceback = Vec::new();
    let (best_score, mut prev_match, mut path) = best_dp;
    //println!("BEST SCORE: {} PREV_MATCH: {}", best_score, prev_match);
    //println!("PATHS");
    while prev_match >= 0 {
        //println!("{} ", prev_match);
        traceback.push(prev_match as usize);
        prev_match = dp[prev_match as usize].1;
    }
    //traceback.reverse();
    //path :traceback;
    //score: best_score,
    //dp_vector: dp,
    best_score*/
}

pub fn find_kmer_matches_for_divided (
    query: &[u8],
    graph_sequences: &Vec<Vec<Vec<u8>>>,
    graph_ids: &Vec<Vec<Vec<usize>>>,
    k: usize
    ) -> (Vec<(u32, u32, usize)>, Vec<u32>, Vec<Vec<(usize, usize)>>, Vec<Vec<u32>>) {
    // hash the query
    let set = hash_kmers(query, k);
    // go through the paths and get the indices of path and make a list with query index, graph index, paths
    let mut kmers_result_vec: Vec<(u32, u32, usize)> = vec![]; //
    let mut kmers_plus_k: Vec<u32> = vec![];
    let mut kmers_paths: Vec<Vec<(usize, usize)>> = vec![];
    let mut kmers_previous_node_in_paths: Vec<Vec<u32>> = vec![];
    for (top_index, section) in graph_sequences.iter().enumerate() {
        for (bottom_index, seq) in section.iter().enumerate() {
            let matches = find_kmer_matches_seq1_hashed(&set, seq, k);
            // go through the matches and see if they are in the already made list
            for a_match in matches {
                // first get the graph index for the match
                let graph_index = graph_ids[top_index][bottom_index][a_match.1 as usize] as u32;
                let graph_index_minus_1;
                let graph_index_plus_k = graph_ids[top_index][bottom_index][a_match.1 as usize - 1 + k] as u32;
                if a_match.1 > 0 {
                    graph_index_minus_1 = graph_ids[top_index][bottom_index][a_match.1 as usize - 1] as u32;
                }
                else {
                    graph_index_minus_1 = u32::MAX;
                }
                // if it is in the result vec, add path if indices match and path different
                if let Some(index_result) = kmers_result_vec.iter().position(|r| (r.0, r.1, r.2) == (a_match.0, graph_index, top_index)) {
                    kmers_paths[index_result].push((top_index, bottom_index));
                    kmers_previous_node_in_paths[index_result].push(graph_index_minus_1);
                }
                // not in the result vec just add
                else {
                    kmers_result_vec.push((a_match.0, graph_index, top_index));
                    kmers_plus_k.push(graph_index_plus_k);
                    kmers_paths.push(vec![(top_index, bottom_index)]);
                    kmers_previous_node_in_paths.push(vec![graph_index_minus_1]);
                }
            }
        }
    }
    
    (kmers_result_vec, kmers_plus_k, kmers_paths, kmers_previous_node_in_paths)
}

pub fn divide_poa_graph_get_paths (output_graph: &POAGraph, topo_indices: &Vec<usize>, total_num_sequences: usize, cut_threshold: usize, topo_map: &HashMap<usize, usize>) -> (Vec<Vec<Vec<usize>>>, Vec<Vec<Vec<u8>>>, usize) {
    let mut cut_start_end: Vec<(usize, usize)>= vec![];
    let mut current_topo_indices_index: usize = 0;
    let mut current_cut_limit: usize = cut_threshold;
    let mut current_start: usize = topo_indices[current_topo_indices_index];
    // add nodes to the processed nodes until we reach the cut limit
    while topo_indices.len() > current_topo_indices_index + 1 {
        if current_topo_indices_index + 1 < topo_indices.len() {
            if current_topo_indices_index >= current_cut_limit {
                println!("Try cutting...");
                let cut_result = try_to_make_the_cut(output_graph, topo_indices[current_topo_indices_index], total_num_sequences);
                // make a cut if successful update the current_cut_limit
                if cut_result {
                    println!("Cut successful! :) at {}", topo_indices[current_topo_indices_index]);
                    // this node becomes the end node make a entry in the vec
                    cut_start_end.push((current_start, topo_indices[current_topo_indices_index]));
                    current_start = topo_indices[current_topo_indices_index + 1];
                    current_cut_limit += cut_threshold;
                }
                // else just add entry to the processed nodes
                else {
                    println!("Cut failed! :( at {}", topo_indices[current_topo_indices_index]);
                    current_cut_limit += 1;
                }
            }
            println!("index {} node {} cut limit {}", current_topo_indices_index, topo_indices[current_topo_indices_index], current_cut_limit);
        }
        else {
            break;
        }
        current_topo_indices_index += 1;
    }
    cut_start_end.push((current_start, topo_indices[current_topo_indices_index]));
    println!("{:?}", cut_start_end);
    let mut all_all_paths: Vec<Vec<Vec<usize>>> = vec![];
    let mut all_all_sequences: Vec<Vec<Vec<u8>>> = vec![];
    let mut max_number_of_paths_per_section: usize = 0;
    for current_start_end in cut_start_end {
        let (start, end) = current_start_end;
        let mut all_paths: Vec<Vec<usize>> = vec![];
        let mut all_sequences: Vec<Vec<u8>> = vec![];
        simple_dfs_with_start_end(output_graph, start, end, vec![], vec![], &mut all_paths, &mut all_sequences, topo_map);
        println!("{:?}", all_paths);
        if max_number_of_paths_per_section < all_paths.len() {
            max_number_of_paths_per_section = all_paths.len()
        }
        all_all_paths.push(all_paths);
        all_all_sequences.push(all_sequences);
    }
    return (all_all_paths, all_all_sequences, max_number_of_paths_per_section);
}

pub fn try_to_make_the_cut(output_graph: &POAGraph, topo_index: usize, total_num_sequences: usize) -> bool {
    // for the selected node, see if there is only one outgoing edge
    let mut neighbour_nodes = output_graph.neighbors_directed(NodeIndex::new(topo_index), Outgoing);
    let mut number_of_outgoing_node = 0;
    let mut number_of_seq_match = false;
    while let Some(neighbour_node) = neighbour_nodes.next() {
        number_of_outgoing_node += 1;
        if number_of_outgoing_node == 1 {
            let mut edges = output_graph.edges_connecting(NodeIndex::new(topo_index), neighbour_node);
            let mut weight: i32 = 0;
            while let Some(edge) = edges.next() {
                weight += edge.weight().clone();
            }
            if weight == total_num_sequences as i32 {
                number_of_seq_match = true;
            }
        }
        if number_of_outgoing_node >= 2 {
            break;
        }
    }
    // if cut is successful return true
    if (number_of_outgoing_node == 1) && (number_of_seq_match) {
        return true;
    }
    else {
        return false;
    }
}

pub fn simple_dfs_with_start_end (
    graph: &POAGraph, 
    start: usize,
    end: usize,
    current_path: Vec<usize>,
    current_sequence: Vec<u8>,
    all_paths: &mut Vec<Vec<usize>>,
    all_sequences: &mut Vec<Vec<u8>>,
    topo_map: &HashMap<usize, usize>
) {
    let mut path = current_path.clone();
    let mut sequence = current_sequence.clone();
    path.push(*topo_map.get(&start).unwrap());
    sequence.push(graph.raw_nodes()[start].weight);
    //println!("current node {}", start);
    // if no neighbours, last node reached
    if (graph.neighbors(NodeIndex::new(start)).count() == 0) || (start == end) {
        all_paths.push(path.clone());
        all_sequences.push(sequence.clone());
    } else {
        // go through neighbours recursively
        for neighbor in graph.neighbors(NodeIndex::new(start)) {
            simple_dfs_with_start_end(graph, neighbor.index(), end, path.clone(), sequence.clone(), all_paths, all_sequences, topo_map);
        }
    }
}

pub fn lcskpp_graph(kmer_pos_vec: Vec<(u32, u32)>, kmers_plus_k: Vec<u32>, kmer_path_vec: Vec<Vec<usize>>, kmers_previous_node_in_paths: Vec<Vec<u32>>, num_of_paths: usize, k: usize) -> u32 {
    // return nothing if empty
    if kmer_pos_vec.is_empty() {
        return 0;
    }

    let k = k as u32;

    let mut events: Vec<(u32, u32, u32, u32, Vec<usize>, Vec<u32>)> = Vec::new();
    let mut max_ns = vec![0; num_of_paths];
    // generate the required events
    for (idx, &(x, y)) in kmer_pos_vec.iter().enumerate() {
        let y_plus_k = kmers_plus_k[idx];
        // add the previous one as well
        // IF END SAVE y - K - 1 NODE IF START SAVE y - 1 NODE (-1 if not available)
        events.push((x, y, y_plus_k, (idx + kmer_pos_vec.len()) as u32, kmer_path_vec[idx].clone(), kmers_previous_node_in_paths[idx].clone()));
        events.push((x + k - 1, y, y_plus_k, idx as u32, kmer_path_vec[idx].clone(), kmers_previous_node_in_paths[idx].clone()));
        for path in kmer_path_vec[idx].clone() {
            max_ns[path] = max(max_ns[path], x + k - 1);
            max_ns[path] = max(max_ns[path], y_plus_k);
        }
    }
    // sorting is okay with topologically converted indices
    events.sort_unstable();
    let mut max_bit_tree_path = vec![];
    // generate empty fenwick trees
    for (_index, n) in max_ns.iter().enumerate() {
        let max_col_dp: MaxBitTree<(u32, u32)> = MaxBitTree::new(*n as usize);
        max_bit_tree_path.push(max_col_dp);
    }
    let mut dp: Vec<(u32, i32)> = Vec::with_capacity(events.len());
    let mut best_dp = (k, 0, 0); // score, coloumn, path
    // update trees only after the next query point is hit
    let mut tree_update_required_level = 0;
    let mut tree_update_vec: Vec<(u32, (u32, u32), u32, Vec<usize>)> = vec![(0, (0, 0), 0, vec![])]; //current x, value,  index, paths (trees)
    dp.resize(events.len(), (0, 0));
    for ev in events {
        if tree_update_required_level != ev.0 && tree_update_vec.len() != 0 {
            // update trees here
            for tree_info in &tree_update_vec {
                for path in &tree_info.3 {
                    max_bit_tree_path[*path].set(tree_info.2 as usize, tree_info.1);
                }
            }
            //tree_update_vec = vec![];
        }
        // p is the match index
        let p = (ev.3 % kmer_pos_vec.len() as u32) as usize;
        //println!("x:{} y_start:{} y_end:{} p:{} idx:{} paths:{:?}", ev.0, ev.1, ev.2, p, ev.3, ev.4);
        // is start if higher than this
        let is_start = ev.3 >= (kmer_pos_vec.len() as u32);
        if is_start {
            //print!("IS START \n");
            dp[p] = (k, -1);
            // go through the paths available in this event, and get the max corrosponding value and pos
            for path_index in 0..ev.4.len() {
                let path = ev.4[path_index];
                let prev_node = ev.5[path_index];
                if prev_node != u32::MAX {
                    //println!("prev node value = {}", prev_node);
                    let (temp_value, temp_position) = max_bit_tree_path[path].get(prev_node as usize);
                    //println!("temp value from fenwick tree {}", temp_value);
                    if (temp_value + k > dp[p].0) && (temp_value > 0) {
                        dp[p] = (k + temp_value, temp_position as i32);
                        best_dp = max(best_dp, (dp[p].0, p as i32, path));
                        //println!("best_dp {}", best_dp.0);
                    }
                }
            }
        } else {
            //print!("IS END \n");
            // See if this kmer continues a different kmer
            if ev.0 >= k {
                for path_index in 0..ev.4.len() {
                    let path = ev.4[path_index];
                    let prev_node = ev.5[path_index];
                    if prev_node != u32::MAX {
                        if let Ok(cont_idx) = kmer_pos_vec.binary_search(&(ev.0 - k, prev_node)) {
                            //println!("!!!!!!!!!!");
                            let prev_score = dp[cont_idx].0;
                            let candidate = (prev_score + 1, cont_idx as i32);
                            dp[p] = max(dp[p], candidate);
                            best_dp = max(best_dp, (dp[p].0, p as i32, path));
                            //println!("candidate location {} {}", ev.0 - k, prev_node);
                            //println!("cont value from candidate p {} score {}", cont_idx, candidate.0);
                            //println!("best_dp {}", best_dp.0);
                        }
                    }
                }
            }
            // set all trees which have this match as this // maybe update this in the next iteration to prevent query overlapping
            //  update required, current x, value,  index, paths (trees)
            tree_update_vec.push((ev.0, (dp[p].0, p as u32), ev.2, ev.4.clone()));
            tree_update_required_level = ev.0;
        }
    }
    let mut traceback = Vec::new();
    let (best_score, mut prev_match, mut path) = best_dp;
    //println!("BEST SCORE: {} PREV_MATCH: {}", best_score, prev_match);
    //println!("PATHS");
    while prev_match >= 0 {
        //println!("{} ", prev_match);
        traceback.push(prev_match as usize);
        prev_match = dp[prev_match as usize].1;
    }
    //traceback.reverse();
    //path :traceback;
    //score: best_score,
    //dp_vector: dp,
    best_score
}

pub fn find_kmer_matches(query: &[u8], graph_sequences: &Vec<Vec<u8>>, graph_ids: &Vec<Vec<usize>>, k: usize) -> (Vec<(u32, u32)>, Vec<u32>, Vec<Vec<usize>>, Vec<Vec<u32>>) {
    // hash the query
    let set = hash_kmers(query, k);
    // go through the paths and get the indices of path and make a list with query index, graph index, paths
    let mut kmers_result_vec: Vec<(u32, u32)> = vec![];
    let mut kmers_plus_k: Vec<u32> = vec![];
    let mut kmers_paths: Vec<Vec<usize>> = vec![];
    let mut kmers_previous_node_in_paths: Vec<Vec<u32>> = vec![];
    for (index, seq) in graph_sequences.iter().enumerate() {
        let matches = find_kmer_matches_seq1_hashed(&set, seq, k);
        // go through the matches and see if they are in the already made list
        for a_match in matches {
            // first get the graph index for the match
            let graph_index = graph_ids[index][a_match.1 as usize] as u32;
            let graph_index_minus_1;
            let graph_index_plus_k = graph_ids[index][a_match.1 as usize - 1 + k] as u32;
            if a_match.1 > 0 {
                graph_index_minus_1 = graph_ids[index][a_match.1 as usize - 1] as u32;
            }
            else {
                graph_index_minus_1 = u32::MAX;
            }
            // if it is in the result vec, add path if indices match and path different
            if let Some(index_result) = kmers_result_vec.iter().position(|r| (r.0, r.1) == (a_match.0, graph_index)) {
                kmers_paths[index_result].push(index);
                kmers_previous_node_in_paths[index_result].push(graph_index_minus_1);
            }
            // not in the result vec just add
            else {
                kmers_result_vec.push((a_match.0, graph_index));
                kmers_plus_k.push(graph_index_plus_k);
                kmers_paths.push(vec![index]);
                kmers_previous_node_in_paths.push(vec![graph_index_minus_1]);
            }
        }
    }
    (kmers_result_vec, kmers_plus_k, kmers_paths, kmers_previous_node_in_paths)
}

pub fn simple_dfs_all_paths (
    graph: &POAGraph, 
    start: usize,
    current_path: Vec<usize>,
    current_sequence: Vec<u8>,
    all_paths: &mut Vec<Vec<usize>>,
    all_sequences: &mut Vec<Vec<u8>>,
    topo_map: &HashMap<usize, usize>
) {
    let mut path = current_path.clone();
    let mut sequence = current_sequence.clone();
    path.push(*topo_map.get(&start).unwrap());
    sequence.push(graph.raw_nodes()[start].weight);
    // if no neighbours, last node reached
    if graph.neighbors(NodeIndex::new(start)).count() == 0 {
        all_paths.push(path.clone());
        all_sequences.push(sequence.clone());
    } else {
        // go through neighbours recursively
        for neighbor in graph.neighbors(NodeIndex::new(start)) {
            simple_dfs_all_paths(graph, neighbor.index(), path.clone(), sequence.clone(), all_paths, all_sequences, topo_map);
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