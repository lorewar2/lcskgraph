use crate::bittree::MaxBitTree;
use fxhash::FxHasher;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
pub type HashMapFx<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

pub fn simple_dfs() {
    // first do this one and check the time complexty and accuracy and stuff
}

pub fn dfs_with_max_link_cuts() {

}

pub fn dfs_with_bfs_check_cuts() {

}

pub fn lcskpp(matches: &[(u32, u32)], k: usize) {
    if matches.is_empty() {
        return;
    }

    let k = k as u32;

    // incoming matches must be sorted to let us find the predecessor kmers by binary search.
    for i in 1..matches.len() {
        assert!(matches[i - 1] < matches[i]);
    }

    let mut events: Vec<(u32, u32, u32)> = Vec::new();
    let mut n = 0;

    for (idx, &(x, y)) in matches.iter().enumerate() {
        events.push((x, y, (idx + matches.len()) as u32));
        events.push((x + k, y + k, idx as u32));

        n = max(n, x + k);
        n = max(n, y + k);
    }
    events.sort_unstable();

    let mut max_col_dp: MaxBitTree<(u32, u32)> = MaxBitTree::new(n as usize);
    let mut dp: Vec<(u32, i32)> = Vec::with_capacity(events.len());
    let mut best_dp = (k, 0);

    dp.resize(events.len(), (0, 0));

    for ev in events {
        let p = (ev.2 % matches.len() as u32) as usize;
        let j = ev.1;
        let is_start = ev.2 >= (matches.len() as u32);

        if is_start {
            dp[p] = (k, -1);
            let (best_value, best_position) = max_col_dp.get(j as usize);
            if best_value > 0 {
                dp[p] = (k + best_value, best_position as i32);
                best_dp = max(best_dp, (dp[p].0, p as i32));
            }
        } else {
            // See if this kmer continues a different kmer
            if ev.0 > k && ev.1 > k {
                if let Ok(cont_idx) = matches.binary_search(&(ev.0 - k - 1, ev.1 - k - 1)) {
                    let prev_score = dp[cont_idx].0;
                    let candidate = (prev_score + 1, cont_idx as i32);
                    dp[p] = max(dp[p], candidate);
                    best_dp = max(best_dp, (dp[p].0, p as i32));
                }
            }

            max_col_dp.set(ev.1 as usize, (dp[p].0, p as u32));
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
    //dp_vector: dp,
}

pub fn find_kmer_matches(seq1: &[u8], seq2: &[u8], k: usize) -> Vec<(u32, u32)> {
    if seq1.len() < seq2.len() {
        let set = hash_kmers(seq1, k);
        find_kmer_matches_seq1_hashed(&set, seq2, k)
    } else {
        let set = hash_kmers(seq2, k);
        find_kmer_matches_seq2_hashed(seq1, &set, k)
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

pub fn find_kmer_matches_seq2_hashed(
    seq1: &[u8],
    seq2_set: &HashMapFx<&[u8], Vec<u32>>,
    k: usize,
) -> Vec<(u32, u32)> {
    let mut matches = Vec::new();

    for i in 0..(seq1.len() + 1).saturating_sub(k) {
        let slc = &seq1[i..i + k];

        if let Some(matches1) = seq2_set.get(slc) {
            for pos1 in matches1 {
                matches.push((i as u32, *pos1));
            }
        }
    }

    matches.sort_unstable();
    matches
}