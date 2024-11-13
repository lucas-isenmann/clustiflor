
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use super::biclust::Biclust;
use super::biclustering::Biclustering;



/// Dont forget to delete the first line of results.txt obtained by Bimax
pub fn load_r_biclusters(
    file_path: &str, 
    node_map_a: &HashMap<String, usize>,  
    node_map_b: &HashMap<String, usize>) -> Biclust {

    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);

    let n = node_map_a.len();
    let m = node_map_b.len();
    let mut biclusters = Biclust::new(n, m);

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
                biclusters.add_bicluster(bicluster.clone());
            }

        }
    }
    biclusters
}