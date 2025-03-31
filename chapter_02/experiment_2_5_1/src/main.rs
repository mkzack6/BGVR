use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use tch::{nn, nn::Module, nn::OptimizerConfig, Device, Tensor, Kind, Reduction};

fn load_csv(csv_path: &str, n_features: usize) -> Result<(Tensor, Tensor), Box<dyn Error>> {
    let file = File::open(csv_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Abaikan header
    let _header = lines.next();

    let mut features: Vec<Vec<f32>> = Vec::new();
    let mut labels: Vec<f32> = Vec::new();

    for line in lines {
        let line = line?;
        let values: Vec<f32> = line
            .split(',')
            .map(|s| s.parse::<f32>().unwrap_or(0.0)) // Default 0.0 jika gagal parse
            .collect();

        if values.len() != n_features + 1 {
            eprintln!("Warning: Skipping row with incorrect feature count: {:?}", values);
            continue;
        }

        let (feature_vec, label) = values.split_at(n_features);
        features.push(feature_vec.to_vec());
        labels.push(label[0]);
    }

    let feature_tensor = Tensor::from_slice2(&features).to_kind(Kind::Float).to(Device::Cpu);
    let label_tensor = Tensor::from_slice(&labels).to_kind(Kind::Float).to(Device::Cpu).unsqueeze(-1); // Tambahkan dimensi batch

    // Normalisasi Min-Max Scaling
    let min_vals = feature_tensor.min_dim(0, true).0;
    let max_vals = feature_tensor.max_dim(0, true).0;
    let feature_tensor = (feature_tensor - &min_vals) / (&max_vals - &min_vals + 1e-8);

    Ok((feature_tensor, label_tensor))
}

#[derive(Debug)]
struct Model {
    fc1: nn::Linear,
    fc2: nn::Linear,
    fc3: nn::Linear,
}

impl nn::Module for Model {
    fn forward(&self, xs: &Tensor) -> Tensor {
        let mut logits = xs.apply(&self.fc1).relu();
        logits = logits.apply(&self.fc2).relu();
        logits = logits.apply(&self.fc3);
        logits // Tanpa sigmoid karena BCEWithLogitsLoss sudah menangani ini
    }
}

fn train_model(train_features: Tensor, train_labels: Tensor, output_path: &str) -> Result<(), Box<dyn Error>> {
    let vs = nn::VarStore::new(Device::Cpu);
    let model = Model {
        fc1: nn::linear(vs.root(), train_features.size()[1], 128, Default::default()),
        fc2: nn::linear(vs.root(), 128, 64, Default::default()),
        fc3: nn::linear(vs.root(), 64, 1, Default::default()),
    };

    let mut opt = nn::Adam::default().build(&vs, 1e-3)?;
    let batch_size = 64;
    let num_epochs = 10;

    let mut output_file = File::create(output_path)?;

    for epoch in 1..=num_epochs {
        let mut total_loss = 0.0;

        for batch_start in (0..train_features.size()[0]).step_by(batch_size) {
            let batch_end = (batch_start + batch_size as i64).min(train_features.size()[0]);
            let batch_xs = train_features.narrow(0, batch_start, batch_end - batch_start);
            let batch_ys = train_labels.narrow(0, batch_start, batch_end - batch_start);

            let logits = model.forward(&batch_xs).squeeze(); // Hapus argumen -1
            let loss = logits.binary_cross_entropy_with_logits::<Tensor>(&batch_ys.squeeze(), None, None, Reduction::Mean); // Hapus argumen -1

            println!(
                "Epoch {}, Batch {}: Logits Mean {:?}, Loss {:?}",
                epoch,
                batch_start,
                logits.mean(Kind::Float),
                loss.mean(Kind::Float)
            );

            opt.backward_step(&loss);
            total_loss += loss.double_value(&[]);
        }

        writeln!(output_file, "Epoch {}, Validation loss: {:.4}", epoch, total_loss)?;
        println!("Epoch {}, Validation loss: {:.4}", epoch, total_loss);
    }

    println!("Training complete. Check output file at: {}", output_path);
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let n_features = 10000_usize;
    let csv_path = "C:\\Users\\trian\\BGVR\\chapter_01\\experiment_15\\src\\gene_expression.csv";
    let output_path = "C:\\Users\\trian\\BGVR\\chapter_01\\experiment_15\\output.txt";

    let (train_features, train_labels) = load_csv(csv_path, n_features)?;

    println!("First 10 feature values: {:?}", train_features.narrow(0, 0, 1));

    train_model(train_features, train_labels, output_path)?;

    Ok(())
}
