# Dockerfile: Containerize naive_variant_caller for HPC nodes
FROM rust:1.82-slim
WORKDIR /app
COPY Cargo.toml .
COPY src ./src
COPY pileup.json .
COPY hypotheses.json .
RUN cargo build --release
CMD ["cargo", "run", "--release", "--", \
     "--pileup-input", "pileup.json", \
     "--hypotheses-input", "hypotheses.json", \
     "--chunk-size", "2", \
     "--output-dir", "partial_variants", \
     "--merged-output", "merged_variants.json"]
