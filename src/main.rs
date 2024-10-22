mod lib;
mod bicluster;
use std::{collections::HashMap, fs::File, io::BufWriter, io::Write};

use bicluster::load_wadj_from_csv;
use lib::{cluster_graph, load_adj_list_file, load_edges_file};

use std::env;
use std::path::Path;


fn main() {

    let args: Vec<String> = env::args().collect();
    
    if args.len() < 3 {
        eprintln!("Usage: {} <file_path> <delimiter>", args[0]);
        std::process::exit(1);
    }

    let file_path = &args[1];
    let delimiter = &args[2];

    println!("File path: {}", file_path);
    println!("Delimiter: {}", delimiter);

    match load_wadj_from_csv(file_path, delimiter) {
        (wadj, n, m, metaA, metaB) => {
            println!("Graph loaded successfully!");
            println!("Number of rows: {}", n);
            println!("Number of cols: {}", m);
            println!("{metaA:?}");
            println!("{metaB:?}");
        }
    }



    // let (matrix, node_map) = load_adj_list_file("real/kegg.txt", ' ');
    // // let matrix = load_edges_file("real/biodata_Shaw/output_file.txt", ' ');
    // // let matrix = generate_random_matrix(10, 0.6);    
    // // println!("{matrix:?}");

    // let n = matrix.shape()[0];
    // let mut m = 0;
    // for i in 0..n {
    //     for j in i+1..n{
    //         if matrix[[i,j]] > 0. {
    //             m += 1;
    //         }
    //     }
    // }
    // println!("n={n} m={m}");
    // let clusters = cluster_graph(matrix, true, 2, 1, 10.);


    // let mut reversed_map = HashMap::new();
    // for (key, value) in node_map.iter() {
    //     reversed_map.insert(*value, key.clone());
    // }

    // let output_file = File::create("clusters.txt").expect("Failed to open file");
    // let mut writer = BufWriter::new(output_file);

    // for cluster in &clusters {
    //     for v in cluster {
    //         write!(writer, "{} ", reversed_map.get(v).unwrap()).expect("Failed to write to file");
    //     }
    //     writeln!(writer).expect("Failed to write end-of-line to file");
    // }
}
