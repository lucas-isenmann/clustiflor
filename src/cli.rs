use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Cli {
    /// cluster|bicluster|onesided-bicluster
    #[arg(value_name = "CLUSTER_TYPE")]
    pub cluster_type: String,

    /// Either a file or a directory
    #[arg(value_name = "DATA_PATH")]
    pub data_path: String,

    /// Split threshold:
    #[arg(short, long, default_value_t = 1.0)]
    pub split_threshold: f64,

    /// Size sensitivity:
    #[arg(long, default_value_t = 1.0)]
    pub size_sensitivity: f64,

    /// Samples size
    #[arg(long, default_value_t = 10)]
    pub samples_size: usize,

    /// Ignore weights (all weights are set to 1)
    #[arg(short, long)]
    pub ignore_weights: bool,

    /// Simple file format
    #[arg(long)]
    pub simple_file_format: bool,

    #[arg(long, default_value_t = 16)]
    pub matrix_power: usize,

    #[arg(long)]
    pub split_rows: bool,

    /// If true, the file format should be a matrix with numerical values
    #[arg(long)]
    pub matrix_format: bool,

    /// CSV separator
    #[arg(long)]
    pub sep: String,

    #[arg(long)]
    pub matrix_labels: bool,

    /// Verbose: 0 means print nothing,
    /// 1 prints details for each iteration,
    /// 2 prints a lot of details
    #[arg(short, long, default_value_t = 0)]
    pub verbose: usize,
}
