use std::collections::HashMap;
use std::hash::Hash;
use std::io::{stdout, Write};

use ndarray::Array2;



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



pub fn progress_bar(current: usize, total: usize) {
    let percentage = (current as f64) / (total as f64);
    let bar_length = 100;
    
    let bar = format!(
        "\r[{:>bar_length$}] {}% ({}/{})",
        ">".repeat((percentage * (bar_length as f64)).floor() as usize),
        percentage.round(),
        current,
        total
    );
    
    stdout().write_all(bar.as_bytes()).unwrap();
}



use std::fs::{self, File};
use std::io::{self, BufRead, BufReader};


pub fn load_r_biclusters(file_path: &str, node_map_a: &HashMap<String, usize>,  node_map_b: &HashMap<String, usize>) -> Vec<Vec<usize>> {

    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);
    let mut biclusters = vec![];
    let mut bicluster = vec![];

    for (i, line) in reader.lines().enumerate() {
        
        if let Ok(line) = line {
            if i % 3 == 0 {
                continue
            }
            else if i % 3 == 1 {
                bicluster.clear();
                let values: Vec<&str> = line.split(" ").collect();
                for  x in values {
                    let nx = node_map_a.get(x).unwrap();
                    bicluster.push(*nx);
                }
            } else {
                let values: Vec<&str> = line.split(" ").collect();
                for  x in values {
                    let nx = node_map_b.get(x).unwrap();
                    bicluster.push(node_map_a.len()+ *nx);
                }
                biclusters.push(bicluster.clone());
            }

        }
    }
    biclusters
}