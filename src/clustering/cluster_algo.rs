use std::collections::HashMap;
use std::env;
use ndarray::Array2;
use rand::rngs::ThreadRng;

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use std::collections::HashSet;
use rand::{seq::SliceRandom, thread_rng};
use std::time::{Instant};

use crate::common::{approx_statio_distrib_by_indegree, compute_statio_distrib_by_exp, compute_statio_distrib_by_iter, compute_statio_distrib_by_pivot, dist, print_matrix, progress_bar};




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





///
/// format
/// ```txt
/// # comments
/// 0 1
/// 1 2
/// 2 0
/// ```
pub fn load_edges_file(file_name: &str, delimiter: char, ignore_weights: bool) -> (Array2<f64>, HashMap<String, usize>) {
    let reader = BufReader::new(File::open(file_name).expect("Failed to open file"));
    let lines = reader.lines();

    let mut data: Vec<(usize, usize, f64)> = Vec::new();
    let mut node_map = std::collections::HashMap::<String, usize>::new();
    let mut n = 0;

    // Read nodes and edges
    for line in lines {
        if let Ok(line) = line {
            if line.starts_with("#"){
                continue;
            }
            let values: Vec<&str> = line.split(delimiter).collect();
            if values.len() >= 2 {
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

                let mut weight = 1.;
                if values.len() >= 3 && ignore_weights == false{
                    weight = values[2].parse().unwrap();

                    if weight < 0.0 || weight > 1.0 {
                        panic!("Weight should be in [0,1]");
                    }
                }
                if weight > 0. {
                    data.push((*n1, *n2, weight));

                }
                // println!("{:?} {ignore_weights}", data.last());
            }
        }
    }

    // Create adjacency matrix
    let mut adj_matrix = Array2::zeros((n, n));

    for (u,v, weight) in data {
        adj_matrix[[u,v]] = weight;
        adj_matrix[[v,u]] = weight;
    }

    (adj_matrix, node_map)
}









fn remove_edge(matrix: &mut Array2<f64>, u: usize, v: usize){
    matrix[[u,v]] = 0.;
    matrix[[v,u]] = 0.
}

// fn add_edge(matrix: &mut Array2<f64>, u: usize, v: usize, weight: f64){
//     matrix[[u,v]] = weight;
//     matrix[[v,u]] = weight
// }





/**
    Given the adjacency matrix, compute the neighbors of vertex `v` at distance at most `d`.
    It is the closed neighborhood.
 */
fn compute_neighbors(matrix: &Array2<f64>, v: usize, d: usize) -> Vec<usize>{
    let n = matrix.nrows();
    let mut neighbors = vec![v];
    let mut last = vec![v];

    for _ in 0..d {
        let mut new_last = Vec::new();
        // let l = last.clone();
        // last.clear();
        for &w in &last {
            for u in 0..n {
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



fn best(matrix: &Array2<f64>, order: &Vec<(usize, f64)>, split_threshold: f64, verbose: usize) -> (f64, Vec<usize>) {
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
            if verbose >= 2{
                println!("cluster: {cluster:?} cost: {c}");
            }
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
            tmc[[nj,ni]] = tm_common[[i_index,j_index]];
            s += tmc[[nj,ni]];
        }
        tmc[[vertex_neighbors_index, ni]] += 1. - s;
    }
    return tmc
}





// static mut DIFF: i128 = 0;




fn compute_order(subset: &Vec<usize>, vertex: usize, tm_common: &Array2<f64>, verbose: usize) -> Vec<(usize, f64)> {  
    let vertex_id = subset.iter().position(|&v| v == vertex).unwrap();

    // Compute centered transition matrix
    let tm = transition_matrix_centered(vertex_id, tm_common, subset);

    let v_result = compute_statio_distrib_by_iter(&tm, 16, verbose);

    // if let Some(v_perfect) = compute_statio_distrib_by_pivot(&tm) {
    //     // println!("{}", subset.len());
    //     // println!("{:.5}", dist(&v_perfect, &v_result));
    //     let t_pivot = Instant::now();
    //     compute_statio_distrib_by_pivot(&tm);
    //     let d_pivot = t_pivot.elapsed().as_millis() as i128;

    //     let mut t_exp = Instant::now();
    //     let v_exp = compute_statio_distrib_by_exp(&tm, 5, 0);
    //     let d_exp = t_exp.elapsed().as_millis() as i128;
    //     let diff_exp = dist(&v_perfect, &v_exp);
    //     // println!("{:.5}", diff_exp);

    //     let mut t_exp = Instant::now();
    //     let v_iter = compute_statio_distrib_by_iter(&tm, 16, 0);
    //     let d_iter = t_exp.elapsed().as_millis() as i128;
    //     let diff_iter = dist(&v_perfect, &v_iter);
    //     // println!("{:.5}", diff_iter);

    //     unsafe { 
    //         // println!("{}", diff_exp - diff_iter);
    //         DIFF += d_approx - d_iter;
    //         println!("{:.5}", DIFF);
    //      };
    // }


    // Order subset by decreasing probability
    let mut order: Vec<(usize, f64)> = subset.iter().enumerate()
        .map(|(ni, &i)| (i, v_result[[ni, 0]]))
        .collect();
    
    order.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    if verbose >= 3 {
        println!("vertex: {vertex}");
        print_matrix(&tm);
        println!("{order:?}");
    }

    order

}



/// Pick a random sample of size k in (0..n-1) such that no indices are in assigned.
/// If k > the number of unassigned vertices, then retain all the unassigned vertices 
fn pick_unassigned_sample(rng: &mut ThreadRng, k: usize, unassigned: &HashSet<usize>) -> Vec<usize> {
    let mut all_numbers: Vec<usize> =vec![];
    for v in unassigned {
        all_numbers.push(*v);
    }
    // all_numbers.retain(|&x| !assigned.contains(&x));
    
    if all_numbers.len() < k {
        return all_numbers
    }


    all_numbers.shuffle(rng);
    let mut result = vec![];
    for i in 0..k {
        result.push(all_numbers[i]);
    }
    result
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

    println!("# Clusters size distribution");
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

pub fn cluster_graph(mut matrix: Array2<f64>, verbose: usize, dist: usize, samples_size: usize, split_threshold: f64) -> Vec<Vec<usize>> {

    let mut rng = thread_rng();

    let n = matrix.nrows();
    let mut unassigned: HashSet<usize> = (0..n).collect();
    let mut c = 0.;

    let mut nb_splits = 0;
    let mut nb_deletions: usize = 0;
    let mut deletions_cost = 0.;
    let mut assignation: Vec<Vec<usize>> = vec![vec![]; n];
    let mut nb_clusters = 0;
    let mut clusters = vec![];

    let mut isolated_vertices = vec!();

    let start_instant: Instant = Instant::now();


    loop {
        if verbose >= 1 {
            println!("-------");
        } else {
            progress_bar(n-unassigned.len(), n, start_instant);
        }

        let mut mindeg = 100000;
        let mut minv = None;
        let mut maxdeg = 0;
        let mut maxv = None;
        isolated_vertices.clear();

        // Search min/max degree and isolated vertices
        for &v in unassigned.iter() {
                let mut degree = 0;
                for j in 0..n {
                    if matrix[[v,j]] > 0.0 {
                        degree += 1
                    }
                }
                if degree == 0 {
                    isolated_vertices.push(v);
                }
                if degree > maxdeg {
                    maxdeg = degree;
                    maxv = Some(v);
                }
                if degree < mindeg {
                    mindeg = degree;
                    minv = Some(v);
                }
        }

        // Cluster isolated vertices
        for &v in isolated_vertices.iter() {
            unassigned.remove(&v);
            assignation[v].push(nb_clusters);
            nb_clusters += 1;
            clusters.push(vec![v]);
        }
        if isolated_vertices.is_empty() == false {
            if verbose >= 1 {
                println!("Nb isolated vertices: {}", isolated_vertices.len());
                println!("Unassigned: {}/{n}", unassigned.len());
            }
            continue;
        }



        let mut best_cost = std::f64::INFINITY;
        let mut best_cluster = vec![];
        let tm = compute_transition_matrix(&matrix, n);

        // Pick a sample
        let mut sample = pick_unassigned_sample(&mut rng,  samples_size, &unassigned);

        // Add minv and maxv if they are not in the sample
        if let Some(minv) = minv {
            if sample.contains(&minv) == false{
                sample.push(minv);
            }
        }
        if let Some(maxv) = maxv {
            if sample.contains(&maxv) == false {
                sample.push(maxv);
            }
        }

        for v in sample {
            let neighbors: Vec<usize> = compute_neighbors(&matrix, v, dist);

            let order: Vec<(usize, f64)> = compute_order(&neighbors, v, &tm, verbose);

            let (cost, cluster) = best(&matrix, &order, split_threshold, verbose);

            if cost < best_cost {
                best_cost = cost;
                best_cluster = cluster;
            }
        }

        


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
            if verbose >= 1 {
                println!("Cluster #{nb_clusters}: {:?}", best_cluster);
                println!("Size: {}", best_cluster.len());
                println!("Reduced Cost: {best_cost:.2}");
            }

            let mut nb_splits_cluster = 0;
            let mut nb_deletions_cluster = 0;

            if best_cluster.len() == 1 {
                for &v in &best_cluster {
                    for j in 0..n {
                        if matrix[[v,j]] > 0. {
                            c += matrix[[v,j]];
                            nb_deletions += 1;
                            deletions_cost += matrix[[v,j]];
                            if verbose >= 1 {
                                println!("del {v} {j} {:.2}", matrix[[v,j]]);
                            }
                            nb_deletions_cluster += 1;
                            remove_edge(&mut matrix, v, j);
                        }
                    }
                }
            }

            for &v in &best_cluster {
                assignation[v].push(nb_clusters);

                // Compute the out degree
                let mut d = 0.0;
                for j in 0..n {
                    if matrix[[v,j]] > 0. && !best_cluster.contains(&j){
                        d += matrix[[v,j]];
                    }
                }

                if d > split_threshold {
                    // Split v
                    nb_splits_cluster += 1;
                    nb_splits += 1;
                    c += split_threshold; 
                    if verbose >= 1 {
                        println!("spl {v}");
                    }
                    for j in best_cluster.iter() {
                        remove_edge(&mut matrix, v, *j)
                    }
                } else {
                    // Delete out edges
                    unassigned.remove(&v);
                    for j in 0..n {
                        if !best_cluster.contains(&j) {
                            c += matrix[[v,j]];
                            if matrix[[v,j]] > 0. {
                                nb_deletions += 1;
                                nb_deletions_cluster += 1;
                                deletions_cost += matrix[[v,j]];
                                if verbose >= 1 {
                                    println!("del {v} {j} {:.2}", matrix[[v,j]]);
                                }
                            }
                        }
                        remove_edge(&mut matrix, v, j);
                    }
                }
            }

            if verbose >= 1 {
                println!("Nb splits: {nb_splits_cluster}");
                println!("Nb deletions: {nb_deletions_cluster}");
                println!("Unassigned vertices remaining: {}/{n}", unassigned.len());
            }
            
            

            nb_clusters += 1;
            clusters.push(best_cluster);
        }
    }

    println!("Clustering ended successfully");


    println!("# Parameters");
    println!("- Split threshold: {split_threshold}");
    println!("- Markov power: {}", 16);
    println!("- Samples size: {samples_size}");
    println!("- Cluster size coef: {}", 1);

    println!("# Results");
    println!("- Nb clusters: {}", clusters.len());
    // println!("Nb_operations: {c}");
    println!("- Nb splits: {nb_splits}" );
    println!("- Nb deletions: {nb_deletions}",  );
    println!("- Deletions Cost: {deletions_cost:.3}",  );
    println!("- Overlapping: {:.3}", 1.+(nb_splits as f64)/(n as f64));

    clusters_size_stats(&clusters);

    clusters
}






pub fn run_cluster_solver() {

    let program_args: Vec<String> = env::args().collect();

    if program_args.len() < 2 {
        eprintln!("Usage: {} cluster <data_path>", program_args[0]);
        std::process::exit(1);
    }

    let mut split_threshold = 1.;
    let mut verbose_level = 0;
    let data_path = &program_args[2];
    let mut samples_size = 10;
    let mut ignore_weights = false;

    println!("Data path: {}", data_path);

    for arg in program_args.iter() {
        if arg.starts_with("--split-threshold=") {
            split_threshold = arg.split_at(18).1.parse::<f64>().unwrap_or(split_threshold);
        }
        if arg.starts_with("--verbose=") {
            verbose_level = arg.split_at(10).1.parse::<usize>().unwrap_or(verbose_level);
        }
        if arg.starts_with("--samples-size=") {
            samples_size = arg.split_at(15).1.parse::<usize>().unwrap_or(samples_size);
        }

        if arg.starts_with("--ignore-weights") {
            ignore_weights = true;
        }
    }



    
    let (matrix, node_indices) = load_edges_file(&data_path, ' ', ignore_weights);


    // Compute the reverse node map
    let mut node_labels = HashMap::new();
    for (key, value) in node_indices.iter() {
        node_labels.insert(*value, key.clone());
    }

    // Compute the number of edges
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
    println!("Start clustering...");
    let clusters = cluster_graph(matrix, verbose_level, 2, samples_size, split_threshold);


    
    let results_path = data_path.to_string() + ".clusters";
    let output_file = File::create(results_path).expect("Failed to create file");
    let mut writer = BufWriter::new(output_file);

    for cluster in &clusters {
        for v in cluster {
            write!(writer, "{} ", node_labels.get(v).unwrap()).expect("Failed to write to file");
        }
        writeln!(writer).expect("Failed to write end-of-line to file");
    }
}

