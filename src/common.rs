use std::io::{stdout, Write};

use ndarray::Array2;
use std::time::{Duration, Instant};


pub fn print_matrix(matrix: &Array2<f64>) {
    let rows = matrix.shape()[0];
    let cols = matrix.shape()[1];

    for i in 0..rows {
        for j in 0..cols {
            let value = matrix[[i, j]];
            let formatted: String = format!("{:.2}", value)
                .replace("NaN", "");
            
            print!("{}", formatted);
            if j < cols - 1 {
                print!(" ");
            }
        }
        println!();
    }
}



pub fn progress_bar(current: usize, total: usize, start_instant: Instant) {
    let elapsed_time = start_instant.elapsed();
    
    let percentage = (current as f64) / (total as f64);
    let bar_length = 100;
    
    let bar = format!(
        "\r[{:>bar_length$}] {}% ({}/{}) | Elapsed: {:>8} | ETA: {:>8}",
        ">".repeat((percentage * (bar_length as f64)).floor() as usize),
        (percentage * 100.).round(),
        current,
        total,
        human_duration(elapsed_time),
        estimate_eta(elapsed_time, percentage)
    );
    
    stdout().write_all(bar.as_bytes()).unwrap();
}

fn estimate_eta(elapsed: Duration, progress: f64) -> String {
    if progress == 0. {
        return String::new();
    }
    let total_duration_estimate = (elapsed.as_secs() as f64) / progress;
    let remaining_estimate = total_duration_estimate - (elapsed.as_secs() as f64);
    let estimated_eta = Duration::from_secs(remaining_estimate as u64);
    human_duration(estimated_eta)
}



/// - Durations greater than 60 seconds are formatted as "Xm Ys" (e.g., "1m 5s")
/// - Durations 60 seconds or less are formatted as "Xs" (e.g., "15s")
///
/// # Examples
/// 
///     use std::time::Duration; 
///     assert_eq!(human_duration(Duration::from_secs(65)), "1m 5s"); 
///     assert_eq!(human_duration(Duration::from_secs(15)), "15s");
/// 
fn human_duration(d: Duration) -> String {
    let total_seconds = d.as_secs();
    
    if total_seconds > 60 {
        let minutes = total_seconds / 60;
        let seconds = total_seconds % 60;
        
        format!("{}m {}s", minutes, seconds)
    } else {
        format!("{}s", total_seconds)
    }
}


