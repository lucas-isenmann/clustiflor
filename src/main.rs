mod lib;
use std::{collections::HashMap, fs::File, io::BufWriter, io::Write};

use lib::{load_adj_list_file, load_edges_file, solve};




fn main() {
    let (matrix, node_map) = load_adj_list_file("real/kegg.txt", ' ');
    // let matrix = load_edges_file("real/biodata_Shaw/output_file.txt", ' ');
    // let matrix = generate_random_matrix(10, 0.6);    
    // println!("{matrix:?}");

    let n = matrix.shape()[0];
    let mut m = 0;
    for i in 0..n {
        for j in i+1..n{
            if matrix[[i,j]] > 0. {
                m += 1;
            }
        }
    }
    println!("n={n} m={m}");
    let clusters = solve(matrix, true, 2, 1, 10.);


    let mut reversed_map = HashMap::new();
    for (key, value) in node_map.iter() {
        reversed_map.insert(*value, key.clone());
    }

    let output_file = File::create("clusters.txt").expect("Failed to open file");
    let mut writer = BufWriter::new(output_file);

    for cluster in &clusters {
        for v in cluster {
            write!(writer, "{} ", reversed_map.get(v).unwrap()).expect("Failed to write to file");
        }
        writeln!(writer).expect("Failed to write end-of-line to file");
    }
}
