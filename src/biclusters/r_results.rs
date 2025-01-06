
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use super::biclust::Biclust;
use super::biclustering::Biclustering;



/// First line of Bimax is the title
pub fn load_r_biclusters(
    file_path: &str, 
    node_map_a: &HashMap<String, usize>,  
    node_map_b: &HashMap<String, usize>) -> Biclust {

    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);

    let n = node_map_a.len();
    let m = node_map_b.len();
    let mut clustered_rows = vec![false; n];
    let mut clustered_cols = vec![false; m];
    let mut biclusters = Biclust::new(n, m);

    let mut bicluster = vec![];

    for (i, line) in reader.lines().enumerate() {
        if i == 0 {
            continue
        }
        
        if let Ok(line) = line {
            if i % 3 == 1 {
                continue
            }
            else if i % 3 == 2 {
                // Rows
                bicluster.clear();
                let values: Vec<&str> = line.split(" ").collect();
                for  x in values {
                    let &nx = node_map_a.get(x).unwrap();
                    clustered_rows[nx] = true;
                    bicluster.push(nx);
                }
            } else {
                // Cols

                let values: Vec<&str> = line.split(" ").collect();
                for  x in values {
                    

                    if let Some(&nx) = node_map_b.get(x) {
                        clustered_cols[nx] = true;
                        bicluster.push(n+ nx);
                    } else {
                        println!("col {line}");

                        println!("{x}");
                        println!("{:?}", node_map_b);
                        panic!("aha");
                    }
                    
                }
                biclusters.add_bicluster(bicluster.clone());
            }

        }
    }
    let mut isolated_rows = vec![];
    for row in 0..n {
        if clustered_rows[row] == false {
            isolated_rows.push(row);
        }
    }
    if isolated_rows.len() > 0 {
        biclusters.add_bicluster(isolated_rows);
    }
    let mut isolated_cols = vec![];
    for col in 0..m {
        if clustered_cols[col] == false {
            isolated_cols.push(n+col);
        }
    }
    if isolated_cols.len() > 0 {
        biclusters.add_bicluster(isolated_cols);
    }
    biclusters
}