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
