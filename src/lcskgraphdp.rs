
//obtained from rustbio library and modified//
use std::cmp::{max, Ordering};
use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;

use petgraph::{Directed, Graph, Incoming};

pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs
pub type POAGraph = Graph<u8, i32, Directed, usize>;

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
#[derive(Debug, Clone)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
}

pub struct Alignment {
    pub score: i32,
    //    xstart: Edge,
    pub operations: Vec<AlignmentOperation>,
}

#[derive(Debug, Clone)]
pub struct TracebackCell {
    dp_score: i32,
    temp_match: i32,
    start_dp: i32,
    op: AlignmentOperation,
}

impl Ord for TracebackCell {
    fn cmp(&self, other: &TracebackCell) -> Ordering {
        self.dp_score.cmp(&other.dp_score)
    }
}

impl PartialOrd for TracebackCell {
    fn partial_cmp(&self, other: &TracebackCell) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TracebackCell {
    fn eq(&self, other: &TracebackCell) -> bool {
        self.dp_score == other.dp_score
    }
}

//impl Default for TracebackCell { }

impl Eq for TracebackCell {}

pub struct Traceback {
    rows: usize,
    cols: usize,

    // store the last visited node in topological order so that
    // we can index into the end of the alignment when we backtrack
    last: NodeIndex<usize>,
    matrix: Vec<(Vec<TracebackCell>, usize, usize)>,
}

impl Traceback {
    /// Create a Traceback matrix with given maximum sizes
    ///
    /// # Arguments
    ///
    /// * `m` - the number of nodes in the DAG
    /// * `n` - the length of the query sequence
    fn initialize_scores(&mut self, gap_open: i32) {
        for j in 0..=self.cols {
            self.matrix[0].0.push(
                TracebackCell {
                    dp_score: (j as i32) * gap_open,
                    temp_match: 0,
                    start_dp: 0,
                    op: AlignmentOperation::Ins(None),
                });
        }
        self.matrix[0].0[0] = TracebackCell {
            dp_score: 0,
            temp_match: 0,
            start_dp: 0,
            op: AlignmentOperation::Match(None),
        };
    }
    fn with_capacity(m: usize, n: usize, gap_open: i32) -> Self {
        // each row of matrix contain start end position and vec of traceback cells
        let mut matrix: Vec<(Vec<TracebackCell>, usize, usize)> = vec![(vec![], 0, n); m + 1];
        for j in 0..=n {
            matrix[0].0.push(TracebackCell {
                dp_score: (j as i32) * gap_open,
                temp_match: 0,
                start_dp: 0,
                op: AlignmentOperation::Ins(None),
            });
        }
        matrix[0].0[0] = TracebackCell {
            dp_score: 0,
            temp_match: 0,
            start_dp: 0,
            op: AlignmentOperation::Match(None),
        };
        Traceback {
            rows: m,
            cols: n,
            last: NodeIndex::new(0),
            matrix,
        }
    }
    fn new_row(&mut self, row: usize, size: usize, gap_open: i32, start: usize, end: usize) {
        //println!("row {} start {} end {}", row, start, end);
        self.matrix[row].1 = start;
        self.matrix[row].2 = end;
        if start == 0 {
            self.matrix[row].0.push(TracebackCell {
                dp_score: 0,
                temp_match: 0,
                start_dp: 0,
                op: AlignmentOperation::Del(None),
            });
        }
        else {
            self.matrix[row].0.push(TracebackCell {
                dp_score: 0,
                temp_match: 0,
                start_dp: 0,
                op: AlignmentOperation::Match(None),
            });
        }
        for _ in 1..=size {
            self.matrix[row].0.push(TracebackCell {
                dp_score: 0,
                temp_match: 0,
                start_dp: 0,
                op: AlignmentOperation::Match(None),
            });
        }
    }
    fn new() -> Self {
        Traceback {
            rows: 0,
            cols: 0,
            last: NodeIndex::new(0),
            matrix: Vec::new(),
        }
    }

    fn set(&mut self, i: usize, j: usize, cell: TracebackCell) {
        // set the matrix cell if in band range
        let real_position = j - self.matrix[i].1;
        self.matrix[i].0[real_position] = cell;
    }

    fn get(&self, i: usize, j: usize) -> &TracebackCell {
        let real_position = j - self.matrix[i].1;
        return &self.matrix[i].0[real_position];
        
    }

    pub fn get_score(&self) -> i32 {
        let i = self.last.index() + 1;
        let (_m, n) = (self.cols, self.cols);
        //println!("dp score {} {} {}", i, n, self.get(i, n).dp_score);
        self.get(i, n).dp_score
    }

    pub fn alignment(&self) -> Alignment {
        // optimal AlignmentOperation path
        let mut ops: Vec<AlignmentOperation> = vec![];
        
        // Now backtrack through the matrix to construct an optimal path
        let mut i = self.last.index() + 1;
        let mut j = self.cols;
        let mut test_vector: Vec<usize> = vec![0; 6];
        while i > 0 || j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(self.get(i,j).op.clone());
            match self.get(i,j).op {
                AlignmentOperation::Match(Some((p, _))) => {
                    i = p + 1;
                    j -= 1;
                    test_vector[0] += 1;
                    //println!("MATCHSOME");
                }
                AlignmentOperation::Del(Some((p, _))) => {
                    i = p + 1;
                    test_vector[1] += 1;
                    //println!("DELSOM");
                }
                AlignmentOperation::Ins(Some(p)) => {
                    i = p + 1;
                    j -= 1;
                    test_vector[2] += 1;
                    //println!("INSSOM");
                }
                AlignmentOperation::Match(None) => {
                    i = 0;
                    j -= 1;
                    test_vector[3] += 1;
                    //println!("MATNON");
                }
                AlignmentOperation::Del(None) => {
                    i -= 1; // j -= 1;
                    test_vector[4] += 1;
                    //println!("DELNON");
                }
                AlignmentOperation::Ins(None) => {
                    j -= 1; // i -= 1;
                    test_vector[5] += 1;
                    //println!("INSNON");
                }
            }
        }
        //println!("MS DS IS MN DN IN {:?}", test_vector);
        ops.reverse();
        Alignment {
            score: self.get(self.last.index() + 1, self.cols).dp_score,
            operations: ops,
        }
    }
}

/// A partially ordered aligner builder
///
/// Uses consuming builder pattern for constructing partial order alignments with method chaining
pub struct Aligner {
    pub traceback: Traceback,
    query: Vec<u8>,
    pub poa: Poa,
}

impl Aligner {
    pub fn new(match_score: i32, mismatch_score: i32, gap_open_score: i32, reference: &Vec<u8>) -> Self {
        Aligner {
            traceback: Traceback::new(),
            query: reference.to_vec(),
            poa: Poa::from_string(match_score, mismatch_score, gap_open_score, reference),
        }
    }

    /// Add the alignment of the last query to the graph.
    pub fn add_to_graph(&mut self) -> &mut Self { 
        let alignment = self.traceback.alignment();
        self.poa.add_alignment(&alignment, &self.query);
        self
    }
    /// Globally align a given query against the graph.
    pub fn global(&mut self, query: &Vec<u8>, k: usize) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.global(query, k);
        let mut topo =  Topo::new(&self.poa.graph);
        let mut topo_indices = Vec::new();
        while let Some(node) = topo.next(&self.poa.graph) {
            topo_indices.push(node);
        }
        self
    }
    pub fn graph(&self) -> &POAGraph {
        &self.poa.graph
    } 
}

/// A partially ordered alignment graph
///
/// A directed acyclic graph datastructure that represents the topology of a
/// traceback matrix.
pub struct Poa {
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    pub graph: POAGraph,
}

impl Poa {
    pub fn from_string(match_score: i32, mismatch_score: i32, gap_open_score: i32, seq: &Vec<u8>) -> Self {
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        let mut node: NodeIndex<usize>;
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);
            graph.add_edge(prev, node, 1);
            prev = node;
        }
        Poa { match_score: match_score, mismatch_score: mismatch_score, gap_open_score: gap_open_score, graph}
    }

    /// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    pub fn global(&self, query: &Vec<u8>, k: usize) -> Traceback {
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        let mut traceback = Traceback::with_capacity(m, n, self.gap_open_score);
        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);
        let mut last_node_match = false;
        // the band required nodes
        while let Some(node) = topo.next(&self.graph) {
            let start = 0;
            let end = query.len();
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            //println!("query len {}", query.len());
            traceback.new_row(i, (end - start) + 1, self.gap_open_score, start, end);
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> = self.graph.neighbors_directed(node, Incoming).collect();
            for (j_p, q) in query.iter().enumerate().skip(start) {
                let j = j_p + 1;
                // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    let temp_score;
                    let dp_score;
                    if r == *q {
                        //println!("MAtch");
                        temp_score = traceback.get(0, j - 1).temp_match + 1;
                    }
                    else {
                        temp_score = 0;
                    }
                    if temp_score >= k as i32 {
                        dp_score = max(temp_score, traceback.get(0, j - 1).dp_score);
                    }
                    else {
                        dp_score = traceback.get(0, j - 1).dp_score;
                    }
                    TracebackCell {
                        dp_score: dp_score,
                        temp_match: temp_score,
                        start_dp: 0,
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = TracebackCell {
                        dp_score: 0,
                        temp_match: 0,
                        start_dp: 0,
                        op: AlignmentOperation::Match(None),
                    };
                    for prev_node in &prevs {
                        let i_p: usize = prev_node.index() + 1; // index of previous node
                        let temp_score;
                        let start_dp;
                        if r == *q {
                            last_node_match = true;
                            //println!("MAtch {}", traceback.get(i_p, j - 1).temp_match);
                            if traceback.get(i_p, j - 1).temp_match == 0 {
                                //temp_score = max(traceback.get(i_p, j - 1).dp_score + 1, 1);
                                temp_score = 1;
                                start_dp = traceback.get(i_p, j - 1).dp_score;
                            }
                            else if (traceback.get(i_p, j - 1).temp_match != 0) && (traceback.get(i_p, j - 1).temp_match + traceback.get(i_p, j - 1).start_dp + k as i32 - 2 < traceback.get(i_p, j - 1).dp_score) {
                                //println!("New start");
                                //temp_score = traceback.get(i_p, j - 1).temp_match + 1;
                                // make this the new start if start_dp calculation is higher than tmp_score + previous start dp
                                temp_score = 1;
                                start_dp = traceback.get(i_p, j - 1).dp_score;
                            }
                            else {
                                temp_score = traceback.get(i_p, j - 1).temp_match + 1;
                                start_dp = traceback.get(i_p, j - 1).start_dp;
                            }
                            //println!("{}", temp_score);
                        }
                        else {
                            last_node_match = false;
                            temp_score = 0;
                            start_dp = 0;
                        }
                        //println!("{}", temp_score);
                        if (temp_score >= k as i32) {
                            //println!("{} {} {} {}", i, j, temp_score, traceback.get(i_p, j - 1).dp_score);
                            //println!("temp {}", temp_score);
                            //println!("{} {} {}", traceback.get(i, j - 1).dp_score, traceback.get(i_p, j).dp_score, temp_score);
                            max_cell = max(
                                max_cell,
                                max(
                                    
                                    TracebackCell {
                                        dp_score: traceback.get(i, j - 1).dp_score,
                                        temp_match: 0,
                                        start_dp: 0,
                                        op: AlignmentOperation::Ins(Some(i - 1)),
                                    },
                                    TracebackCell {
                                        dp_score: temp_score + start_dp,
                                        temp_match: temp_score,
                                        start_dp: start_dp,
                                        op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                    }
                                ),
                            );
                            //println!("max dp {} temp {}", max_cell.dp_score, max_cell.temp_match);
                            max_cell = max(
                                max(
                                    TracebackCell {
                                        dp_score: traceback.get(i_p, j).dp_score,
                                        temp_match: 0,
                                        start_dp: traceback.get(i_p, j).start_dp,
                                        op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                    },
                                    TracebackCell {
                                        dp_score: traceback.get(i, j - 1).dp_score,
                                        temp_match: 0,
                                        start_dp: 0,
                                        op: AlignmentOperation::Ins(Some(i - 1)),
                                    },
                                ),
                                max_cell,
                            );
                            
                        }
                        else {
                            //println!("{} {} ", i, j);
                            max_cell = max(
                                max_cell,
                                max(
                                    TracebackCell {
                                        dp_score: traceback.get(i_p, j).dp_score,
                                        temp_match: temp_score,
                                        start_dp: start_dp,
                                        op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                    },
                                    TracebackCell {
                                        dp_score: traceback.get(i, j - 1).dp_score,
                                        temp_match: temp_score,
                                        start_dp: start_dp,
                                        op: AlignmentOperation::Ins(Some(i - 1)),
                                    },
                                ),
                            );
                            max_cell = max(
                                max_cell,
                                max(
                                    TracebackCell {
                                        dp_score: traceback.get(i_p, j).dp_score,
                                        temp_match: temp_score,
                                        start_dp: start_dp,
                                        op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                    },
                                    TracebackCell {
                                        dp_score: traceback.get(i, j - 1).dp_score,
                                        temp_match: temp_score,
                                        start_dp: start_dp,
                                        op: AlignmentOperation::Ins(Some(i - 1)),
                                    },
                                ),
                            );
                        }
                    }
                    max_cell
                };
                let score = max_cell;
                //println!("{} {} {} {}", i, j , score.dp_score, score.temp_match);
                traceback.set(i, j, score);
            }
        }
        //println!("{}", last_node_match);
        /*if last_node_match {
            let maxcell = TracebackCell {
                dp_score: traceback.get( traceback.last.index() + 1, query.len()).dp_score + 1,
                temp_match: traceback.get( traceback.last.index() + 1, query.len()).temp_match,
                start_dp: traceback.get( traceback.last.index() + 1, query.len()).start_dp,
                op: AlignmentOperation::Ins(Some(1)),
            };
            traceback.set(traceback.last.index() + 1, query.len(), maxcell);
        }*/
        let last_row = traceback.last.index() + 2;
        let last_coloumn = query.len() + 1;
        // print the whole matrix for debugging stuff
        for i in 0.. last_coloumn {
            for j in 0.. last_row {
                print!("*{} {} {}* ", traceback.get(j, i).temp_match, traceback.get(j, i).dp_score, traceback.get(j, i).start_dp);
            }
            println!();
        }
        
        traceback
    }
    pub fn custom(&mut self, query: &Vec<u8>) -> Traceback {
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        // save score location of the max scoring node for the query for suffix clipping
        let mut max_in_column = vec![(0, 0); n + 1];
        let mut traceback = Traceback::with_capacity(m, n, -2);
        traceback.initialize_scores(self.gap_open_score);
        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);
        while let Some(node) = topo.next(&self.graph) {
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1; // 0 index is for initialization so we start at 1
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            traceback.new_row(
                i,
                n + 1,
                self.gap_open_score,
                0,
                n + 1,
            );
            // query base and its index in the DAG (traceback matrix rows)
            for (query_index, query_base) in query.iter().enumerate() {
                let j = query_index + 1; // 0 index is initialized so we start at 1
                                         // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    let temp_score;
                    if r == *query_base {
                        temp_score = self.match_score;
                    }
                    else {
                        temp_score = self.mismatch_score;
                    }
                    TracebackCell {
                        dp_score: traceback.get(0, j - 1).dp_score + temp_score,
                        start_dp: 0,
                        temp_match: 0,
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = 
                        TracebackCell {
                            dp_score: MIN_SCORE,
                            start_dp: 0,
                            temp_match: 0,
                            op: AlignmentOperation::Match(None)};
                    for prev_node in &prevs {
                        let i_p: usize = prev_node.index() + 1; // index of previous node
                        let temp_score;
                        if r == *query_base {
                            temp_score = self.match_score;
                        }
                        else {
                            temp_score = self.mismatch_score;
                        }
                        max_cell = max(
                            max_cell,
                            max(
                                TracebackCell {
                                    dp_score: traceback.get(i_p, j - 1).dp_score
                                        + temp_score,
                                    start_dp: 0,
                                    temp_match: 0,
                                    op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    dp_score: traceback.get(i_p, j).dp_score + self.gap_open_score,
                                    start_dp: 0,
                                    temp_match: 0,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                    }
                    max_cell
                };
                let score = max(
                    max_cell,
                    TracebackCell {
                        dp_score: traceback.get(i, j - 1).dp_score + self.gap_open_score,
                        start_dp: 0,
                        temp_match: 0,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                traceback.set(i, j, score);
            }
        }

        traceback
    }
    /// Incorporate a new sequence into a graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment of the new sequence to the graph
    /// * `seq` - The sequence being incorporated
    pub fn add_alignment(&mut self, aln: &Alignment, seq: &Vec<u8>) {
        let head = Topo::new(&self.graph).next(&self.graph).unwrap();
        let mut prev: NodeIndex<usize> = NodeIndex::new(head.index());
        let mut i: usize = 0;
        let mut edge_not_connected: bool = false;
        for op in aln.operations.iter() {
            match op {
                AlignmentOperation::Match(None) => {
                    let node: NodeIndex<usize> = NodeIndex::new(head.index());
                    if (seq[i] != self.graph.raw_nodes()[head.index()].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        if edge_not_connected {
                            self.graph.add_edge(prev, node, 1);
                        }
                        edge_not_connected = false;
                        prev = node;
                    }
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                        edge_not_connected = false;
                    }
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(*p);
                    if (seq[i] != self.graph.raw_nodes()[*p].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                    } else {
                        // increment node weight
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                *self.graph.edge_weight_mut(edge).unwrap() += 1;
                            }
                            None => {
                                if prev.index() != head.index() && prev.index() != node.index() {
                                    self.graph.add_edge(prev, node, 1);
                                }
                            }
                        }
                        prev = NodeIndex::new(*p);
                    }
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {
                    let node = self.graph.add_node(seq[i]);
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                    }
                    prev = node;
                    edge_not_connected = true;
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes and xclip
            }
        }
    }
}
