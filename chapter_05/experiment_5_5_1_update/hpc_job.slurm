#!/bin/bash
#SBATCH --job-name=variant_caller
#SBATCH --output=variant_caller_%j.out
#SBATCH --error=variant_caller_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00

PROJECT_DIR=/home/zack/BGVR/chapter_05/t1/naive_variant_caller
CHUNK_SIZE=2
export PATH="$HOME/.cargo/bin:$PATH"

if ! command -v cargo >/dev/null 2>&1; then
    echo "Error: Cargo not found"
    exit 1
fi

cd $PROJECT_DIR || exit 1

# Run for each chunk
for CHUNK_INDEX in 0 1; do
    cargo run --release -- \
      --pileup-input pileup.json \
      --hypotheses-input hypotheses.json \
      --chunk-size $CHUNK_SIZE \
      --chunk-index $CHUNK_INDEX \
      --output-dir partial_variants \
      --merged-output merged_variants.json
done
