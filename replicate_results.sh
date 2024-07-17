cargo run --release -- -s -t 20 -l 100 -b 10 -k 2
cargo run --release -- -s -t 20 -l 100 -b 10 -k 4
cargo run --release -- -s -t 20 -l 1000 -b 50 -k 4
cargo run --release -- -s -t 20 -l 1000 -b 50 -k 8
cargo run --release -- -s -t 20 -l 10000 -b 200 -k 8
cargo run --release -- -s -t 20 -l 10000 -b 200 -k 12
cargo run --release -- -s -t 20 -l 30000 -b 300 -k 8
cargo run --release -- -s -t 20 -l 30000 -b 300 -k 12
cargo run --release -- -p "data/sample_pacbio.bam" "" -t 20 -b 200 -k 8
cargo run --release -- -p "data/sample_pacbio.bam" "" -t 20 -b 200 -k 12