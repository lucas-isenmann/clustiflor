use std::collections::HashMap;
use ndarray::Array2;

use std::fs::File;
use std::io::{BufRead, BufReader};


use std::collections::HashSet;
use rand::{seq::SliceRandom, thread_rng};
use std::time::Instant;




pub fn load_adj_list_file(file_name: &str, delimiter: char) -> (Array2<f64>, HashMap<String, usize>) {
    let reader = BufReader::new(File::open(file_name).expect("Failed to open file"));
    let lines = reader.lines();

    let mut data: Vec<(usize, usize)> = Vec::new();
    let mut node_map = std::collections::HashMap::<String, usize>::new();
    let mut n = 0;

    // Read nodes and edges
    for line in lines {
        if let Ok(line) = line {
            if line.starts_with("#"){
                continue;
            }
            let values: Vec<&str> = line.split(delimiter).collect();
            let v = String::from(values[0]);
            if !node_map.contains_key(&v) {
                node_map.insert(v.clone(), n);
                n += 1;
            }
            let nv = node_map.get(&v).unwrap().clone();

            for i in 1..values.len(){
                let neighbor = String::from(values[i]);
                if !node_map.contains_key(&neighbor) {
                    node_map.insert(neighbor.clone(), n);
                    n += 1;
                }
                let n2 = node_map.get(&neighbor).unwrap();
                data.push((nv, *n2))
                
            }
        }
    }

    // Create adjacency matrix
    let mut adj_matrix = Array2::zeros((n, n));

    for (u,v) in data {
        adj_matrix[[u,v]] = 1.;
        adj_matrix[[v,u]] = 1.;
    }

    (adj_matrix, node_map)
}






pub fn load_edges_file(file_name: &str, delimiter: char) -> Array2<f64> {
    let reader = BufReader::new(File::open(file_name).expect("Failed to open file"));
    let lines = reader.lines();

    let mut data: Vec<(usize, usize)> = Vec::new();
    let mut node_map = std::collections::HashMap::<String, usize>::new();
    let mut n = 0;

    // Read nodes and edges
    for line in lines {
        if let Ok(line) = line {
            let values: Vec<&str> = line.split(delimiter).collect();
            if values.len() == 2 {
                let node1 = String::from(values[0]);
                let node2 = String::from(values[1]);

                // Add both nodes to the map if not already present
                if !node_map.contains_key(&node1) {
                    node_map.insert(node1.clone(), n);
                    n += 1;
                }
                if !node_map.contains_key(&node2) {
                    node_map.insert(node2.clone(), n);
                    n += 1;
                }

                let n1 = node_map.get(&node1).unwrap();
                let n2 = node_map.get(&node2).unwrap();
                data.push((*n1, *n2))
            }
        }
    }

    // Create adjacency matrix
    let mut adj_matrix = Array2::zeros((n, n));

    for (u,v) in data {
        adj_matrix[[u,v]] = 1.;
        adj_matrix[[v,u]] = 1.;
    }

    adj_matrix
}









fn remove_edge(matrix: &mut Array2<f64>, u: usize, v: usize){
    matrix[[u,v]] = 0.;
    matrix[[v,u]] = 0.
}

fn add_edge(matrix: &mut Array2<f64>, u: usize, v: usize, weight: f64){
    matrix[[u,v]] = weight;
    matrix[[v,u]] = weight
}





/**
    Given the adjacency matrix, compute the neighbors of vertex `v` at distance at most `d`.
    It is the closed neighborhood.
 */
fn compute_neighbors(matrix: &Array2<f64>, v: usize, d: usize) -> Vec<usize>{
    let mut neighbors = vec![v];
    let mut last = vec![v];

    for _ in 0..d {
        let mut new_last = Vec::new();
        // let l = last.clone();
        // last.clear();
        for &w in &last {
            for u in 0..matrix.shape()[0] {
                if matrix[[w,u]] > 0. && !neighbors.contains(&u) {
                    neighbors.push(u);
                    new_last.push(u)
                }
            }
        }
        last = new_last;
    }
    
    neighbors 
}


fn compute_transition_matrix(matrix: &Array2<f64>, n: usize) -> Array2<f64> {
    let mut tm: Array2<f64> = Array2::zeros((n, n));

    let mut degree = vec![0.0; n];

    for i in 0..n {
        for j in 0..n {
            if i != j {
                degree[i] += matrix[[i,j]];
            }
        }
    }

    for i in 0..n {
        for j in 0..n {
            if degree[i] > 0.0 {
                tm[[i,j]] = matrix[[i,j]] / degree[i];
            }
        }
    }

    tm
}



fn best(matrix: &Array2<f64>, order: &Vec<(usize, f64)>, split_threshold: f64) -> (f64, Vec<usize>) {
    let n = matrix.shape()[0];
    let mut best_weight = std::f64::INFINITY;
    let mut best_cluster = Vec::new();
    let mut cluster = Vec::new();
    let mut problematic_pairs: HashSet<(usize, usize)> = HashSet::new();
    let mut neighbors: Vec<Vec<usize>> = vec![vec![]; n];
    let mut in_cluster = vec![false; n];
    // let mut problematic_pairs: Vec<(usize, usize)> = Vec::new();

    let mut outdegree = vec![0.0; n];
    let mut c: f64 = 0.;

    for &(i,_) in order.iter() {
        in_cluster[i] = true;
        // Remove problematic pairs
        problematic_pairs.retain(|&pair| !(matrix[[pair.0,i]] > 0.0 && matrix[[pair.1,i]] > 0.0));
    
        // // Add new problematic pairs
        for &v in cluster.iter() {
            if matrix[[i,v]] == 0.0 {
                let mut common = false;
                for &w in neighbors[v].iter() {
                    if matrix[[i,w]] > 0.0 {
                        common = true;
                        break;
                    }
                }
                // for &w in cluster.iter(){
                //     if matrix[[w,v]] > 0.0 && matrix[[i,w]] > 0.0 {
                //         common = true;
                //         break;
                //     }
                // }
                if !common {
                    // problematic_pairs.push((i,v));
                    problematic_pairs.insert((i, v));
                }
            }
        }
        cluster.push(i);
    
        // Handle single-node clusters
        if cluster.len() == 1 {
            outdegree[i] = 0.0;
            for j in 0..n {
                outdegree[i] += matrix[[i,j]];
            }
            if outdegree[i] > split_threshold {
                c += split_threshold;
            } else {
                c += outdegree[i];
            }
            best_weight = outdegree[i];
            best_cluster = vec![i];
            // println!("outdegree {i} {} cost: {c}", outdegree[i]);
            continue;
        }
    
        // Update outdegrees
        outdegree[i] = 0.0;
        for j in 0..n {
            if  matrix[[i,j]] > 0.0 {
                if in_cluster[j] == false {
                // if !cluster.contains(&j){
                    outdegree[i] += matrix[[i,j]];
                } else{
                    neighbors[i].push(j);
                    neighbors[j].push(i);
                    if outdegree[j] > split_threshold && outdegree[j] - matrix[[i,j]] <= split_threshold {
                        c -= split_threshold;
                        c += outdegree[j] - matrix[[i,j]];
                    } else if outdegree[j] <= split_threshold {
                        c -= matrix[[i,j]];
                    }
                    outdegree[j] -= matrix[[i,j]];
                }
            }
        }

        //
        if outdegree[i] > split_threshold {
            c += split_threshold;
        } else {
            c += outdegree[i]
        }

    
        // Check if there are any problematic pairs
        if problematic_pairs.is_empty() {
            // Calculate cost
            let cost = c* (cluster.len() as f64).powf(-1.0);
            if cost == 0.0 {
                return (cost, cluster);
            }
            if cost <= best_weight {
                // println!("prefix: {cluster:?} {c} {cost:.2}");

                best_weight = cost;
                best_cluster = cluster.clone();
            }
        } else {
            // println!("Non 2club: {}", cluster.len())
        }
    }
    
    return (best_weight, best_cluster);
    
}




fn transition_matrix_centered(vertex_neighbors_index: usize, tm_common: &Array2<f64>, neighbors: &Vec<usize>) -> Array2<f64> {
    
    let d = neighbors.len();
    let mut tmc: Array2<f64> = Array2::zeros((d, d));

    for ni in 0..d{
        let i = neighbors[ni];
        let i_index = i;
        let mut s = 0.;

        for nj in 0..d{
            let j = neighbors[nj];
            let j_index = j;
            tmc[[ni,nj]] = tm_common[[i_index,j_index]];
            s += tmc[[ni,nj]];
        }
        tmc[[ni,vertex_neighbors_index]] += 1. - s;
    }
    return tmc
}



fn compute_order(subset: &Vec<usize>, vertex: usize, tm_common: &Array2<f64>) -> Vec<(usize, f64)> {  
    let d = subset.len();
    let vertex_id = subset.iter().position(|&v| v == vertex).unwrap();

    // Compute centered transition matrix
    let tm = transition_matrix_centered(vertex_id, tm_common, subset);

    // Compute tm^8*v_0
    let mut tm_powered = tm.t().into_owned();
    tm_powered = tm_powered.dot(&tm_powered);
    tm_powered = tm_powered.dot(&tm_powered);
    // tm_powered = tm_powered.dot(&tm_powered);
    
    let mut v = Array2::zeros((d, 1));
    v[[vertex_id, 0]] = 1.0;
    
    let v_result = tm_powered.dot(&v);

    // Order subset by decreasing probability
    let mut order: Vec<(usize, f64)> = subset.iter().enumerate()
        .map(|(ni, &i)| (i, v_result[[ni, 0]]))
        .collect();
    
    order.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    order

}



fn choose_unique_random_numbers(n: usize, k: usize, assigned: &HashSet<usize>) -> Vec<usize> {
    let mut all_numbers: Vec<usize> = (0..n).collect();
    
    all_numbers.retain(|&x| !assigned.contains(&x));
    
    if all_numbers.len() < k {
        return all_numbers
    }

    let mut rng = thread_rng();
    all_numbers.shuffle(&mut rng);

    all_numbers.drain(..k).collect()
}



fn clusters_size_stats(clusters: &Vec<Vec<usize>>) {
    let mut sorted_vectors: Vec<Vec<usize>> = clusters.clone();

    sorted_vectors.sort_by_key(|v| -((*v).len() as i32));

    // println!("Lengths of inner vectors in descending order:");
    // for v in &sorted_vectors {
    //     println!("{}", v.len());
    // }

    // Calculate the average length
    let total_length: usize = sorted_vectors.iter().map(|v| v.len()).sum();
    let average_length = total_length as f64 / clusters.len() as f64;

    println!("Clusters size distribution");
    println!("- Average size: {:.2}", average_length);
    println!("- Max size: {}", sorted_vectors[0].len());
    // println!("- Second max size: {}", sorted_vectors[1].len());
    // println!("- Third max size: {}", sorted_vectors[2].len());
    println!("- Min size: {}", sorted_vectors[sorted_vectors.len()-1].len());
}



/// Solves a clustering problem.
///
/// # Arguments
///
/// - `matrix`: A 2D array representing a graph with weights between 0 and 1.
/// - `verbose`: A boolean flag to control verbosity of output.
/// - `dist`: An unsigned integer specifying the maximal inner distance of clusters.
/// - `nb_vertices`: An unsigned integer specifying the number of random vertices considered at each step.
/// - `split_threshold`: A float value representing the split threshold.
///
/// # Returns
///
/// A vector of vectors, where each inner vector contains the indices of the biclusters found.
///
/// # Examples
///
/// ```

pub fn cluster_graph(mut matrix: Array2<f64>, verbose: bool, dist: usize, samples_size: usize, split_threshold: f64) -> Vec<Vec<usize>> {

    let n = matrix.shape()[0];
    // println!("n = {n} {}", original_vertices.len());
    let mut assigned = HashSet::new();
    let mut c = 0.;

    let mut nb_splits = 0;
    let mut assignation: Vec<Vec<usize>> = vec![vec![]; n];
    let mut nb_clusters = 0;
    let mut clusters = vec![];

    loop {
        println!("-------\n");
        // println!("{} {}", matrix[[2381, 47]], matrix[[47, 2381]]);
        let mut best_cost = std::f64::INFINITY;
        let mut best_cluster = vec![];

        let start = Instant::now();
        let tm = compute_transition_matrix(&matrix, n);
        let duration = Instant::now().duration_since(start);
        println!("Transition matrix {:.6} seconds", duration.as_millis() as f64 / 1000.0);


        let mut mindeg = 100000;
        let mut minv = None;
        let mut maxdeg = 0;
        // let mut maxv = None;

        for v in 0..n {
            if !assigned.contains(&v) {
                let mut degree = 0;
                for j in 0..n {
                    if matrix[[v,j]] > 0.0 {
                        degree += 1
                    }
                }
                if degree > maxdeg {
                    maxdeg = degree;
                    // maxv = Some(v);
                }
                if degree < mindeg {
                    mindeg = degree;
                    minv = Some(v);
                }
            }
        }

        if let Some(minv) = minv {
            println!("Center: {minv} degree: {mindeg}");
            // let x = compute_neighbors(&matrix, minv, 1);

            // println!("vertex: {minv} N[v]: {x:?} {}", x.len());


            let start = Instant::now();
            let x1 = compute_neighbors(&matrix, minv, dist);
            let duration = Instant::now().duration_since(start);
            println!("Compute neighbors {:.6} seconds", duration.as_millis() as f64 / 1000.0);

            println!("Try size: {}", x1.len());
            // println!("vertex: {minv} N[v]: {x1:?} {}", x1.len());

            let start = Instant::now();
            let order1 = compute_order(&x1, minv, &tm);
            let duration = Instant::now().duration_since(start);
            println!("order {:.6} seconds", duration.as_millis() as f64 / 1000.0);

            // println!("Order: {order1:?}");
            let start = Instant::now();
            let (cost1, cluster1) = best(&matrix, &order1, split_threshold);
            let duration = Instant::now().duration_since(start);
            println!("best {:.6} seconds", duration.as_millis() as f64 / 1000.0);

            if cost1 < best_cost {
                best_cost = cost1;
                best_cluster = cluster1;
            }
        }



        // Samples
        let sample = choose_unique_random_numbers(n , samples_size-1, &assigned);

        for v in sample {
            let start = Instant::now();
            let subset: Vec<usize> = compute_neighbors(&matrix, v, dist);
            let duration = Instant::now().duration_since(start);
            println!("Compute neighbors {:.6} seconds", duration.as_millis() as f64 / 1000.0);

            println!("Subset: {}", subset.len());
            // println!("vertex: {minv} N[v]: {x1:?} {}", x1.len());

            let start = Instant::now();
            let order1: Vec<(usize, f64)> = compute_order(&subset, v, &tm);
            let duration = Instant::now().duration_since(start);
            println!("order {:.6} seconds", duration.as_millis() as f64 / 1000.0);

            // println!("Order: {order1:?}");
            let start = Instant::now();
            let (cost1, cluster1) = best(&matrix, &order1, split_threshold);
            let duration = Instant::now().duration_since(start);
            println!("best {:.6} seconds", duration.as_millis() as f64 / 1000.0);

            if cost1 < best_cost {
                best_cost = cost1;
                best_cluster = cluster1;
            }
        }

        


        // Commented out part for max vertex (as in original Python code)
        /*
        if let Some(maxv) = maxv {
            let x1 = compute_2neighbors(&graph, maxv);
            let order1 = compute_order(&x1, maxv, &tm, &original_vertices);
            let (cost1, cluster1) = best(&graph, order1.clone());
            if cost1 < best_cost {
                best_cost = cost1;
                best_cluster = cluster1;
            }
        }
        */

        // Commented out part for checking every non-assigned vertex (as in original Python code)
        /*
        for &v in original_vertices.iter() {
            if !assigned.contains(&v) {
                let x = compute_2neighbors(&graph, v);
                let order = compute_order(&x, v, &tm, &original_vertices);
                let (cost, cluster) = best(&graph, order.clone());
                if cost == 0.0 {
                    best_cost = 0.0;
                    best_cluster = cluster;
                    break;
                }
                if cost < best_cost {
                    best_cost = cost;
                    best_cluster = cluster;
                }
            }
        }
        */

        if best_cluster.is_empty() {
            break;
        } else {
            if verbose {
                println!("Cluster #{nb_clusters}: {:?} (size: {})", best_cluster, best_cluster.len());
                println!("Cost: {best_cost}");
            }

            let mut nb_splits_cluster = 0;

            if best_cluster.len() == 1 {
                for &v in &best_cluster {
                    for j in 0..n {
                        if matrix[[v,j]] > 0. {
                            c += matrix[[v,j]];
                            remove_edge(&mut matrix, v, j);
                        }
                    }
                }
            }

            for &v in &best_cluster {
                assignation[v].push(nb_clusters);
                // println!("check {v}");

                let mut d = 0.0;
                for j in 0..n {
                    if matrix[[v,j]] > 0. && !best_cluster.contains(&j){
                        // println!("{v} {j}");
                        d += matrix[[v,j]];
                    }
                }
                // println!("outdegree[{v}] = {d}");

                if d > split_threshold {
                    // Split v
                    nb_splits_cluster += 1;
                    nb_splits += 1;
                    c += split_threshold; 
                    // println!("spl {v}");
                    for j in best_cluster.iter() {
                        remove_edge(&mut matrix, v, *j)
                    }
                } else {
                    // Delete out edges
                    assigned.insert(v);
                    for j in 0..n {
                        if !best_cluster.contains(&j) {
                            c += matrix[[v,j]];
                            // println!("del {v} {j}");
                        }
                        remove_edge(&mut matrix, v, j);
                    }
                }
            }

            println!("Nb splits: {nb_splits_cluster}/{}", best_cluster.len());
            println!("Assigned: {}/{n}", assigned.len());
            

            nb_clusters += 1;
            clusters.push(best_cluster);
        }
    }

    let nb_deletions =c - split_threshold*(nb_splits as f64);

    println!("# Parameters");
    println!("Split threshold: {split_threshold}");
    println!("Markov power: {}", 4);
    println!("Samples size: {samples_size}");
    println!("Cluster size coef: {}", 1);

    println!("# Results");
    println!("Nb clusters: {}", clusters.len());
    // println!("Nb_operations: {c}");
    println!("Nb splits: {nb_splits}" );
    println!("Nb deletions: {nb_deletions:.0}",  );
    println!("Overlapping: {:.3}", 1.+(nb_splits as f64)/(n as f64));

    clusters_size_stats(&clusters);

    clusters
}



