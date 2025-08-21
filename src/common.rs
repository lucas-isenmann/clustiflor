use std::{io::{stdout, Write}};
use clap::Parser;
use colored::Colorize;
use ndarray::Array2;
use std::time::{Duration, Instant};

extern crate terminal_size;
use terminal_size::{terminal_size, Width};

// extern crate terminal_size::{Width, Height, terminal_size};

fn swap_columns(matrix: &mut Array2<f64>, col1: usize, col2: usize) {
   for i in 0..matrix.nrows(){
        let t = matrix[[i,col2]];
        matrix[[i,col2]] = matrix[[i,col1]];
        matrix[[i,col1]] = t;
   }
}

fn swap_rows(matrix: &mut Array2<f64>, row1: usize, row2: usize) {
   for j in 0..matrix.ncols(){
        let t = matrix[[row1,j]];
        matrix[[row1,j]] = matrix[[row2,j]];
        matrix[[row2,j]] = t;
   }
}

fn eliminate_below_pivot(matrix: &mut Array2<f64>, pivot_row: usize, row: usize) {
    let factor = matrix[[row, pivot_row]];
    for j in 0..matrix.ncols() {
        matrix[[row, j]] = matrix[[row, j]] - matrix[[pivot_row, j]] * factor;
    }
}



pub fn compute_statio_distrib_by_pivot(matrix: &Array2<f64>) -> Option<Array2<f64>> {
    let n = matrix.nrows() as usize;
    let mut u = matrix.clone();
    for i in 0..n {
        u[[i,i]] -= 1.;
    }

    
    // Pivot algorithm implementation
    for i in 0..n {
        // Find pivot element
        let max_row = (i..n).max_by(|&a, &b| u[[a, i]].abs().partial_cmp(&u[[b, i]].abs()).unwrap()).unwrap();

        // Check if column is all zeros
        if u.column(i).iter().all(|&x| x.abs() < 1e-10) {
            // Try next columns
            for j in i+1..n {
                if !u.column(j).iter().all(|&x| x.abs() < 1e-10) {
                    // Swap columns
                    swap_columns(&mut u, i, j);
                    break;
                }
            }
            // If we couldn't find a non-zero column, matrix might be full rank
            if u.column(i).iter().all(|&x| x.abs() < 1e-10) {
                continue;
            }
            
            // update max_row
            // max_row = (i..n).max_by(|&a, &b| u[[a, i]].abs().partial_cmp(&u[[b, i]].abs()).unwrap()).unwrap();
        }
        
        // Swap rows if needed
        if max_row != i {
            swap_rows(&mut u, i, max_row);
        }
        
        // Normalize pivot row
        let pivot = u[[i, i]];
        if pivot.abs() < 1e-10 {
            continue;
        }
        u.row_mut(i).mapv_inplace(|x| x / pivot);
        
        // Eliminate column entries above and below pivot
        for j in 0..n {
            if j != i {
                eliminate_below_pivot(&mut u, i, j);
            }
        }
    }

    // Find free variables
    let free_vars: Vec<usize> = (0..n).filter(|&i| u[[i,i]].abs() < 1e-10).collect();
    
    if free_vars.is_empty() {
        return None; // Matrix is full rank
    }
    
    // Construct kernel vector
    let mut kernel = Array2::zeros((n,1));
    kernel[[free_vars[0],0]] = 1.0;
    
    // Back-substitution
    for i in (0..n).rev() {
        if !free_vars.contains(&i) {
            let sum: f64 = free_vars.iter().map(|j| u[[i, *j]] * kernel[[*j,0]]).sum();
            kernel[[i,0]] = -sum;
        }
    }

    let mut sum = 0.;
    for i in 0..n {
        sum += kernel[[i,0]];
    }

    for i in 0..n {
        kernel[[i,0]] /= sum;
    }
    
    Some(kernel)
}




pub fn compute_statio_distrib_by_exp(tm: &Array2<f64>, exp: usize, verbose: usize) -> Array2<f64> {
    let mut tm_powered = tm.clone();
    let d = tm.nrows();

    for _ in 0..exp {
        tm_powered = tm_powered.dot(&tm_powered);
    }
    
    let mut v = Array2::zeros((d, 1));
    let p = 1./(d as f64);
    for i in 0..d {
        v[[i, 0]] = p;
    }
    // v[[vertex_id, 0]] = 1.0;

    
    
    let result = tm_powered.dot(&v);
    if verbose >= 3 {
        println!("{result:?}");
    }
    result
}


/// Fastest way to compute the stationnary distribution of a transition matrix.
/// 2x faster than by exponentiating the matrix
/// 4x faster than by pivot
/// 
pub fn compute_statio_distrib_by_iter(tm: &Array2<f64>, exp: usize, verbose: usize) -> Array2<f64>  {  
    let d = tm.nrows();
    
    let mut v = Array2::zeros((d, 1));
    let p = 1./(d as f64);
    for i in 0..d {
        v[[i, 0]] = p;
    }
    
    for _ in 0..exp {
        v = tm.dot(&v);
    }
    if verbose >= 3 {
        println!("{v:?}");
    }
    v
}

pub fn approx_statio_distrib_by_indegree(tm: &Array2<f64>, verbose: usize) -> Array2<f64>  {  
    let d = tm.nrows();
    let mut v = Array2::zeros((d, 1));
    
    let mut total = 0.;
    for i in 0..d {
        let mut indegree = 0.;
        for j in 0..d {
            indegree += tm[[i,j]];
        }
        v[[i, 0]] = indegree;
        total += indegree;
    }

    for i in 0..d {
        v[[i,0]] /= total;
    }
   
    if verbose >= 3 {
        println!("{v:?}");
    }
    v
}

// pub fn compute_1_eigenvector(tm: &Array2<f64>, verbose: usize)  {  

//     // Compute centered transition matrix
//     let tm = transition_matrix_centered(vertex_id, tm_common, subset);
//     let mut tm = tm.t().into_owned();

    
//     let result = TruncatedSvd::new(tm, TruncatedOrder::Smallest)
//         .decompose(1)
//         .unwrap();

//     let (u, sigma, v_t) = result.values_vectors();

//     println!("{u} {sigma} {v_t}");
    


// }



pub fn dist(v: &Array2<f64>, w: &Array2<f64>) -> f64{
    let mut sum = 0.;
    for i in 0..v.nrows() {
        sum +=  (v[[i,0]] - w[[i,0]])* (v[[i,0]] - w[[i,0]])
    }
    return f64::sqrt(sum)
}


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
    let mut bar_length: usize = 50;

    if let Some((Width(w), _)) = terminal_size() {
        if w >= 55 {
            bar_length = (w-55) as usize;
        } else {
            bar_length = 0;
        }
    } else {
        println!("Unable to get term size :(")
    }
    
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



pub fn print_error(message: &str){
    eprintln!("{} {message}", "error:".red().bold());
    std::process::exit(1);
}
  



  /// Clustering and biclustering algorithms with overlaps
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

    #[arg(long)]
    pub matrix_format: bool,


    /// Verbose: 0 means print nothing,
    /// 1 prints details for each iteration,
    /// 2 prints a lot of details 
    #[arg(short, long, default_value_t = 0)]
    pub verbose: usize,

   
}
