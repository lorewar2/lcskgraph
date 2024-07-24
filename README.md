# LCSKPOA

LCSKPOA: Enabling banded semi-global partial order alignments via efficient and accurate backbone generation through extended lcsk++.

## Setting up / Requirements

Rust 1.65.0+ should be installed.

## Usage

Download the repository.

```bash
git clone https://github.com/lorewar2/lcskgraph.git
cd LCSKGRAPH
```

Replicate the results in paper:

```bash
bash replicate_results.sh
```

Run bam file:

```bash
cargo run --release -- -p "data/sample_pacbio.bam" "output.fa" -b 200 -k 12
```

Run fasta file:

```bash
cargo run --release -- -p "data/pacbio.fa" "output.fa" -b 200 -k 12
```

## How to cite

If you are using LCSKPOA in your work, please cite:

[LCSKPOA: Enabling banded semi-global partial order alignments via efficient and accurate backbone generation through extended lcsk++](https://www.biorxiv.org/content/10.1101/2024.07.18.604181v1)