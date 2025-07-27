use core::f64;
use std::collections::HashMap;
use std::time::Instant;

use ndarray::Array2;


use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::biclusters::biclust::Biclust;
use crate::common::{print_matrix, progress_bar};

use super::common::transition_matrix_b;
use super::weighted_biadj::WeightedBiAdjacency;





/// The center vertex is supposed to be at index 0 in neighbors (so it is the closed neighborhood)
/// 
fn centered_transition_matrix(  tm: &Array2<f64>, neighbors: &Vec<usize>) -> Array2<f64> {
    let d = neighbors.len();
    let mut ctm = Array2::zeros((d,d));

    for ni in 0..d {
        let i = neighbors[ni];
        let mut s = 0.;
        for nj in 0..d {
            let j = neighbors[nj];
            ctm[[ni,nj]] = tm[[i,j]];
            s += ctm[[ni,nj]];
        }
        ctm[[ni,0]] += 1.-s;
    }

    ctm
}

fn compute_order(wadj: &WeightedBiAdjacency, m: usize, vertex: usize, tm_common: &Array2<f64>, markov_power: usize, verbose: usize)
 -> Vec<(usize,f64)> {
    
    // Search the B vertices which have a common neighbor with vertex
    // neighbors[0] = vertex
    let mut neighbors: Vec<usize> = vec![vertex];
    for b in 0..m {
        if b == vertex {
            continue;
        }
        for (x,_) in wadj.iter(b){
            let mut found = false;
            for (y,_) in wadj.iter(vertex) {
                if y == x {
                    neighbors.push(b);
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
    }
    let d = neighbors.len();
    let mut ctm = centered_transition_matrix(tm_common, &neighbors).t().into_owned();
    if verbose >= 2 {
        println!("Centered to {vertex} transition matrix:");
        print_matrix(&ctm);
    }


    for _ in 0..markov_power {
        ctm = ctm.dot(&ctm);
    }
    

    let mut v = Array2::zeros((d, 1));
    v[[0, 0]] = 1.0;
    
    let v_result = ctm.dot(&v);

    // Order subset by decreasing probability
    let mut order: Vec<(usize, f64)> = neighbors.iter().enumerate()
        .map(|(ni, &i)| (i, v_result[[ni, 0]]))
        .collect();
    
    order.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    order

}

#[derive(PartialEq, Eq)]
enum State {
    Negative,
    Positive
}



fn degree(wadj: &WeightedBiAdjacency, a: usize, m: usize) -> f64 {
    let mut d = 0.;
    for (_,w) in wadj.iter(a+m) {
        d += w;
    }
    d
}

///
/// cost_coef in [0,1]
fn best(n:usize, m: usize, order: Vec<(usize, f64)>, wadj: &WeightedBiAdjacency, cost_coef: f64, split_threshold: f64) -> (f64, Vec<usize>)  {
    let mut best_cost = f64::NAN;
    let mut best_cluster = Vec::new();
    let mut b_cluster = Vec::new();

    let mut s = 0.;
    let mut c = 0.;
    let mut state = HashMap::new();
    let mut indegree = vec![0.;n];
    let mut outdegree = vec![0.;n];

    for (i, proba_i) in order {
        if proba_i == 0. {
            break
        }
        b_cluster.push(i);
        s += 1.;

        for (x, st) in state.iter_mut(){
            if wadj.has_edgee(*x, i) == false {
                if *st == State::Positive {
                    if indegree[*x] <= s*0.5 {
                        *st = State::Negative;
                        // state.insert(*x, State::Negative);
                        c += indegree[*x];
                        c -= (s-1.) - indegree[*x];
                        if outdegree[*x] <= split_threshold {
                            c -= outdegree[*x];
                        } else {
                            c -= 1.;
                        }
                    } else {
                        c += 1.
                    }
                }
            }
        }
        for (&x,&w) in wadj.iter(i) {
            if state.contains_key(&x) == false {
                indegree[x] = w;
                outdegree[x] = degree(wadj, x, m) - w;
                if indegree[x] > s*0.5 {
                    state.insert(x, State::Positive);
                    c += s - indegree[x];
                    if outdegree[x] > split_threshold {
                        c += 1.
                    } else {
                        c += outdegree[x];
                    }
                } else {
                    state.insert(x, State::Negative);
                    c += indegree[x];
                }
            } else if  state[&x] == State::Positive {
                indegree[x] += w;
                outdegree[x] -= w;

                if indegree[x] > s*0.5 {
                    // Keep connected
                    c += 1.-w;
                    if outdegree[x]+w > 1. && outdegree[x] <= split_threshold {
                        c += -1. + outdegree[x];
                    } else if outdegree[x] <= 1. {
                        c -= w;
                    }
                } else {
                    // Unconnect x
                    state.insert(x, State::Negative);
                    c -= (s-1.) - (indegree[x]-w);
                    if outdegree[x] + w > split_threshold {
                        c -= 1.;
                    } else {
                        c -= outdegree[x] + w;
                    }
                    c += indegree[x];
                }
            } else {
                indegree[x] += w;
                outdegree[x] -= w;
                if indegree[x] > s*0.5 {
                    state.insert(x, State::Positive);
                    c -= indegree[x] - w;
                    c += s - indegree[x];
                    if outdegree[x] > split_threshold {
                        c += 1.;
                    } else {
                        c += outdegree[x];
                    }
                } else {
                    c += w;
                }
            }
        }

        let cost = c*f64::powf(s, -cost_coef);
        if best_cost.is_nan() || cost < best_cost {
            best_cost = cost;
            best_cluster = b_cluster.clone();
        }
    }

    (best_cost, best_cluster)
}






fn compute_unclustered_a(n: usize, a_clusters: &Vec<Vec<usize>>) -> Vec<usize> {
    let mut a_unclustered = vec![];
    for a in 0..n {
        let mut found = false;
        for a_cluster in a_clusters.iter() {
            if a_cluster.contains(&a) {
                found = true;
                break;
            }
        }
        if found == false {
            a_unclustered.push(a);
        }
    }

    a_unclustered
}


pub struct AlgoStats {
    pub adjusted_error: f64,
    pub nb_operations: f64,
    pub nb_splits: f64,
    pub nb_deletions: f64,
    pub nb_additions: f64
}



pub fn bicluster( wadj: &mut WeightedBiAdjacency, cost_coef: f64, split_threshold: f64, markov_power: usize, verbose: usize) -> (Biclust, AlgoStats) {

    let min_error = wadj.compute_min_error();


    let n = wadj.get_n();
    let m = wadj.get_m();
    let mut nb_operations = 0.;
    let mut nb_deletions = 0.;
    let mut nb_splits = 0.;
    let mut nb_additions = 0.;

    let mut a_clusters = vec![];
    let mut b_clusters = vec![];
    let mut isolated_b_vertices = vec![];

    let mut assigned = vec![false; m];
    let mut nb_assigned = 0;

    let start_instant = Instant::now();

    // While there exists some B vertices to cluster
    loop {
        
        if verbose >= 1 {
            println!("{nb_assigned} B / {m} vertices are assigned");
        }
        if nb_assigned == m {
            break;
        }
        
        if verbose == 0 {
            progress_bar(nb_assigned, m, start_instant);
        }

        // Check if there is a B vertex of degree 0
        let mut has_isolated_b_vertices = false;
        for b in 0..m {
            if assigned[b] == false && wadj.col_degree(b) == 0 {
                assigned[b] = true;
                nb_assigned += 1;
                isolated_b_vertices.push(b);
                has_isolated_b_vertices = true;
            }
        }
        if has_isolated_b_vertices {
            continue;
            
        }

        // Compute the transition matrix between B vertices
        let tm = transition_matrix_b(wadj, n, m);

        if verbose >= 2 {
            println!("Transition matrix:");
            print_matrix(&tm);
        }

        let mut best_cluster = Vec::new();
        let mut best_cost = f64::NAN;


        // Find the B_cluster with minimal cost
        for b in 0..m {
            if assigned[b] {
                continue;
            }

            // Compute order of b
            let mut d = 0.;
            for (_,w) in wadj.iter(b) {
                d += w;
            }
            let order = 
                if d == 0. {
                    vec![(b,1.)]
                } else {
                    compute_order(&wadj, m, b, &tm, markov_power, verbose)
                };
            if verbose >= 2 {
                println!("Step 1: compute order of {b}");
                println!("{order:?}");
            }
            
            // Compute best cost
            let (cost, b_cluster) = best(n, m, order, &wadj, cost_coef, split_threshold);
            if best_cost.is_nan() || cost < best_cost {
                best_cost = cost;
                best_cluster = b_cluster.clone();
            }
            if verbose >= 2 {
                println!("Step 2: best_cost: {best_cost}, cluster: {best_cluster:?}");
            }
        }


        if best_cluster.len() == 0 { // This case should not happen
            break
        }
        
        for &b in best_cluster.iter() {
            if assigned[b] {
                println!("bug {b} already assigned");
            }
            assigned[b] = true;
            nb_assigned += 1;
        }
        b_clusters.push(best_cluster.clone());

        if verbose >= 2 {
            println!("Step 3: apply {best_cluster:?}");
        }
        let (a_cluster, del, add, spl )= apply_operations(n, m, best_cluster, wadj, split_threshold, verbose); 
        nb_deletions += del;
        nb_additions += add;
        nb_splits += spl;
        nb_operations += del + add + spl;

        a_clusters.push(a_cluster);
    }

    // Add isolated B vertices
    if isolated_b_vertices.len() > 0 {
        b_clusters.push(isolated_b_vertices.clone());
        a_clusters.push(vec![]);
    }

    // After having clustered every B vertices, it is possible that there remains unclustered A vertices
    let unclustered_a = compute_unclustered_a(n,  &a_clusters);

    if unclustered_a.len() > 0 {
        a_clusters.push(unclustered_a.clone());
        b_clusters.push(vec![]);
    }

    let eta = (nb_additions + nb_deletions ) / ((n *m) as f64);
    let eta = eta-min_error;
    

    let bicluster_stats  = AlgoStats {
        adjusted_error: eta,
        nb_operations: nb_operations,
        nb_splits: nb_splits,
        nb_additions: nb_additions,
        nb_deletions: nb_deletions
    };

    // Check integrity
    check_integrity(&a_clusters, &b_clusters, n, m);

    
    // compute_clusters(&a_clusters, &b_clusters, n, m)
    (Biclust::from_separate_biclusters(n,m, &a_clusters, &b_clusters), bicluster_stats)

}      










fn apply_operations(n: usize, m: usize, b_cluster: Vec<usize>, wadj: &mut WeightedBiAdjacency, split_threshold: f64, verbose: usize) -> (Vec<usize>, f64, f64, f64) {
    let mut a_cluster = vec![];
    let mut nb_deletions = 0.;
    let mut nb_splits = 0.;
    let mut nb_additions = 0.;
    let blen = b_cluster.len() as f64;

    for a in 0..n {
        let mut indegree = 0.;
        for (i,w) in wadj.iter(a+m) {
            if b_cluster.contains(&i) {
                indegree += w;
            }
        } 
        if indegree > blen*0.5 {
            a_cluster.push(a);
            nb_additions += blen - indegree;

            for &b in b_cluster.iter() {
                if wadj.has_edgee(a, b) {
                    wadj.delete_edge(a, b);
                }   
            }

            // Compute the number of neighbors of vertex 'a' outside of X
            let mut out_degree = 0.;
            for (i,w) in wadj.iter(a+m) {
                if b_cluster.contains(&i) == false {
                    out_degree += w;
                }
            }

            // Delete case
            if out_degree <= split_threshold {
                nb_deletions += out_degree;
                for b in 0..m {
                    if b_cluster.contains(&b) == false && wadj.has_edgee(a, b) {
                        if verbose >= 2 {
                            println!("delete {a} {b}")
                        }
                        wadj.delete_edge(a, b);
                    }
                }
            } 
            // Split case
            else { 
                nb_splits += 1.;
                if verbose >= 2 {
                    println!("split {a}")
                }
                for &b in b_cluster.iter() {
                    if wadj.has_edgee(a, b) {
                        wadj.delete_edge(a, b);
                    }
                }
            }



        }
        else {
            nb_deletions += indegree;
            for &b in b_cluster.iter() {
                if wadj.has_edgee(a, b) {
                    if verbose >= 2 {
                        println!("delete {a} {b}")
                    }
                    wadj.delete_edge(a, b);
                }
            }
        }
    }

    (a_cluster, nb_deletions, nb_additions, nb_splits)
}



fn compute_overlapping(a_clusters: &Vec<Vec<usize>>, n: usize) -> f64 {
    let mut sum = 0.;
    for a_cluster in a_clusters {
        sum += a_cluster.len() as f64;
    }
    sum / (n as f64)
}



fn compute_clusters(a_clusters: &Vec<Vec<usize>>, b_clusters: &Vec<Vec<usize>>, n: usize, m: usize) -> Vec<Vec<usize>> {
    let mut clusters = vec![];
    for i in 0..b_clusters.len() {
        let mut cluster = vec![];
        for &b in b_clusters[i].iter() {
            cluster.push(b+n);
        }
        for &a in a_clusters[i].iter() {
            cluster.push(a);
        }
        cluster.sort();
        clusters.push(cluster);
    }
    clusters
}


/// Check if the b_clusters form a partition of {0, ..., m-1}
/// and that a_clusters are covering {0, ..., n-1}.
fn check_integrity(a_clusters: &Vec<Vec<usize>>, b_clusters: &Vec<Vec<usize>>, n: usize, m: usize) -> bool {

    // Check B_clusters is a partition of {0, ..., m-1}  (every integer is covered exactly once)
    let mut b_hit = vec![0; m];
    for b_cluster in b_clusters {
        for &b in b_cluster {
            b_hit[b] += 1;
        }
    }
    for b in 0..m {
        if b_hit[b] != 1 {
            println!("ERROR: B vertex {b} is hit {} times", b_hit[b]);
            return false;
        }
    }

    // Check that all A vertices are clustered
    let mut a_hit = vec![0; n];
    for a_cluster in a_clusters {
        for &a in a_cluster {
            a_hit[a] += 1;
        }
    }
    for a in 0..n {
        if a_hit[a] == 0 {
            println!("ERROR: A vertex: {a} is not hit");
            return false;
        }
    }

    true
}

    









pub fn load_wadj_from_csv(file_path: &str, del: &str, split_rows: bool) -> (WeightedBiAdjacency, usize, usize, Vec<String>, Vec<String>, HashMap<String, usize>, HashMap<String, usize>) {

    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);

    let mut n = 0;
    let mut m = 0;
    let mut node_map_a = HashMap::<String, usize>::new();
    let mut node_map_b = HashMap::<String, usize>::new();
    let mut labels_a = Vec::new();
    let mut labels_b = Vec::new();
    let mut edges: Vec<(usize, usize, f64)> = vec![];

    for line in reader.lines() {
        if let Ok(line) = line {
            if line.starts_with('#') {
                continue;
            }


            let values: Vec<&str> = line.split(del).collect();
            
            if values.len() <= 1 {
                eprintln!("Error: Line has less than 2 columns. Skipping.");
                continue; 
            }


            let node1 = if split_rows {
                String::from(values[0])
            } else {
                String::from(values[1])
            };
            let node2 = if split_rows {
                String::from(values[1])
            } else {
                String::from(values[0])
            };

            // Add both nodes to the maps if not already present
            if !node_map_a.contains_key(&node1) {
                node_map_a.insert(node1.clone(), n);
                labels_a.push(node1.clone());
                n += 1;
            }
            if !node_map_b.contains_key(&node2) {
                node_map_b.insert(node2.clone(), m);
                labels_b.push(node2.clone());
                m += 1;
            }

            let n1 = node_map_a.get(&node1).unwrap();
            let n2 = node_map_b.get(&node2).unwrap();

            let mut weight = 1.;
            if values.len() >= 3 {
                weight = values[2].parse().unwrap();
            }
            edges.push((*n1, *n2, weight))
        }
    }

    let mut wadj = WeightedBiAdjacency::new(n, m);
    for &(a,b,w) in edges.iter() {
        wadj.add_edge(a, b, w);
    }



    (wadj, n, m, labels_a, labels_b, node_map_a, node_map_b)
}


/// 
/// 
pub fn load_wadj_from_matrix(file_path: &str) -> WeightedBiAdjacency {
    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);

    let mut n = 0;
    let mut m = 0;
    let mut edges: Vec<(usize, usize, f64)> = vec![];
    let mut biclusters: Vec<Vec<usize>> = vec![];

    for line in reader.lines() {
        if let Ok(line) = line {

             // Skip empty lines and comments
             if line.is_empty() || line.starts_with('#') {
                continue;
            }
            
            // Parse bicluster definition
            if line.starts_with("BC") {
                let members_str = &line[3..];
                if members_str.len() == 0{
                    continue;
                }
                let members: Vec<usize> = members_str
                    .split(',')
                    .map(|x| x.parse().expect("Invalid bicluster member"))
                    .collect();
                biclusters.push(members);
                continue;
            }

            let values: Vec<&str> = line.split(" ").collect();
            m = values.len();

            for j in 0..m {
                let weight: f64 = values[j].parse().unwrap();
                if weight != 0. {
                    edges.push((n, j, weight));
                }
            }
            n += 1;

        }
    }

    let mut wadj = WeightedBiAdjacency::new(n, m);
    for &(a,b,w) in edges.iter() {
        wadj.add_edge(a, b, w);
    }

    wadj.set_ground_biclusters(Biclust::from_biclusters(n, m, &biclusters));

    wadj
}




pub fn print_wadj_stats(wadj: &WeightedBiAdjacency, n: usize, m: usize){



    println!("# Data statistics");
    println!("- n: {n}");
    println!("- m: {m}");

    let mut minw = f64::NAN;
    let mut maxw = f64::NAN;
    let mut ne = 0;
    let mut density = 0.;
    for b in 0..m {
        ne += wadj.col_degree(b);
        for (_, &w) in wadj.iter(b) {
            if maxw.is_nan() || w > maxw {
                maxw = w;
            }
            if minw.is_nan() || w < minw {
                minw = w;
            }
            density += w;
        }
    }
    density /= (n*m) as f64;
    println!("- Number of edges: {ne}");
    println!("- Density of edges: {:.3}", (ne as f64)/((n*m) as f64));

    println!("- Min weight: {minw}");
    println!("- Max weight: {maxw}");
    println!("- Avg weight: {density:.3}");
    println!("- Min error: {:.3}", wadj.compute_min_error());


}


fn biclusters_contains_subset(biclusters: &Vec<Vec<usize>>, subset: &Vec<usize>) -> bool{
    for bicluster in biclusters {
        let mut ok = true;
        for x in subset {
            if bicluster.contains(x) == false {
                ok = false;
                break;
            }
        }
        if ok {
            return true;
        }
    }
    false
}


pub fn compute_edition_diff(biclusters: &Vec<Vec<usize>>, wadj: &Vec<HashMap<usize, f64>>, n: usize, m: usize) -> f64 {
    let mut r = 0;

    for a in 0..n {
        for b in 0..m {
            if wadj[b].contains_key(&a) != biclusters_contains_subset(biclusters, &vec![a,b+n]) {
                println!("{a} {} ",b+n);
                r += 1;
            }
        }
    }
    
    
    r as f64
}



pub fn compute_nb_unclustered(biclusters: &Vec<Vec<usize>>, n: usize, m: usize) -> (usize, usize) {
    let mut ra = 0;
    for a in 0..n {
        if biclusters_contains_subset(biclusters, &vec![a]) == false {
            ra += 1;
        }
    }    
    let mut rb = 0;
    for b in 0..m {
        if biclusters_contains_subset(biclusters, &vec![b+n]) == false {
            rb += 1;
        }
    }
    (ra, rb)
}