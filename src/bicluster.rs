use core::f64;
use std::collections::HashMap;

use ndarray::Array2;


/// Return the transition matrix between B vertices
/// 
fn transition_matrix_b(wadj: &Vec<HashMap<usize, f64>>, n: usize, m: usize) -> Array2<f64> {
    let mut tm = Array2::zeros((m,m));

    let d = compute_degrees(wadj, n, m);

    // Compute tm[i][j] = probability i -> j
    // tm[i][j] = (sum wxi wxj) / (d[i] d[j])
    for i in 0..m{
        for j in 0..m {
            for (x,wxi) in wadj[i].iter() {
                for (y, w) in wadj[j].iter(){
                    if x == y {
                        tm[[i,j]] += (wxi*w)/(d[i]*d[x+m]);
                        break;
                    }
                }
            }
        }
    }
    tm

}


/// The center vertex is supposed to be at index 0
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

fn compute_order(wadj: &Vec<HashMap<usize, f64>>, m: usize, vertex: usize, tm_common: &Array2<f64>, markov_power: usize)
 -> Vec<(usize,f64)> {
    
    // Search the B vertices which have a common neighbor with vertex
    // neighbors[0] = vertex
    let mut neighbors: Vec<usize> = vec![vertex];
    for b in 0..m {
        if b == vertex {
            continue;
        }
        for (x,_) in wadj[b].iter() {
            let mut found = false;
            for (y,_) in wadj[vertex].iter() {
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
    let vertex_id = 0; // index in neighbors
    let ctm = centered_transition_matrix( tm_common, &neighbors);

    let mut tm_powered = ctm.t().into_owned();
    for _ in 0..markov_power {
        tm_powered = tm_powered.dot(&tm_powered);
    }
    

    let mut v = Array2::zeros((d, 1));
    v[[vertex_id, 0]] = 1.0;
    
    let v_result = tm_powered.dot(&v);

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



fn degree(wadj: &Vec<HashMap<usize, f64>>, a: usize, m: usize) -> f64 {
    let mut d = 0.;
    for (_,w) in wadj[a+m].iter() {
        d += w;
    }
    d
}

///
/// cost_coef in [0,1]
fn best(m: usize, order: Vec<(usize, f64)>, wadj: &Vec<HashMap<usize, f64>>, cost_coef: f64) -> (f64, Vec<usize>)  {
    let mut best_weight = f64::NAN;
    let mut best_cluster = Vec::new();
    let mut b_cluster = Vec::new();

    let mut s = 0.;
    let mut c = 0.;
    let mut state = HashMap::new();
    let mut indegree = vec![0.;m];
    let mut outdegree = vec![0.;m];

    for (i, proba_i) in order {
        if proba_i == 0. {
            break
        }
        b_cluster.push(i);
        s += 1.;

        for (x, st) in state.iter_mut(){
            if wadj[i].contains_key(x) == false {
                if *st == State::Positive {
                    if indegree[*x] <= s*0.5 {
                        *st = State::Negative;
                        // state.insert(*x, State::Negative);
                        c += indegree[*x];
                        c -= (s-1.) - indegree[*x];
                        if outdegree[*x] <= 1. {
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
        for (&x,&w) in wadj[i].iter() {
            if state.contains_key(&x) == false {
                indegree[x] = w;
                outdegree[x] = degree(wadj, x, m) - w;
                if indegree[x] > s*0.5 {
                    state.insert(x, State::Positive);
                    c += s - indegree[x];
                    if outdegree[x] > 1. {
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
                    if outdegree[x]+w > 1. && outdegree[x] <= 1. {
                        c += -1. + outdegree[x];
                    } else if outdegree[x] <= 1. {
                        c -= w;
                    }
                } else {
                    // Unconnect x
                    state.insert(x, State::Negative);
                    c -= (s-1.) - (indegree[x]-w);
                    if outdegree[x] + w > 1. {
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
                    if outdegree[x] > 1. {
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
        if best_weight.is_nan() || cost < best_weight {
            best_weight = cost;
            best_cluster = b_cluster.clone();
        }
    }

    (best_weight, best_cluster)
}



fn compute_degrees(wadj: &Vec<HashMap<usize, f64>>, n: usize, m: usize) -> Vec<f64> {
    let mut degrees = vec![0.; n+m];
    for i in 0..m {
        for (_,w) in wadj[i].iter() {
            degrees[i] += w;
        }
    }
    for x in 0..n {
        for (_,w) in wadj[x+m].iter() {
            degrees[x+m] += w;
        }
    }
    degrees
}



fn compute_unclustered_a(wadj: &Vec<HashMap<usize, f64>>, n: usize, m: usize, a_clusters: &Vec<Vec<usize>>) -> Vec<usize> {
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



pub fn bicluster( wadj: &mut Vec<HashMap<usize, f64>>, n:usize, m: usize, cost_coef: f64, split_coef: f64, markov_power: usize) -> Vec<Vec<usize>> {
    let mut nb_operations = 0.;
    let mut nb_deletions = 0.;
    let mut nb_splits = 0.;
    let mut nb_additions = 0.;

    let mut a_clusters = vec![];
    let mut b_clusters = vec![];

    let mut assigned = vec![false; m];

    // While there exists B vertices to cluster
    loop {
        
        // Compute the transition matrix between B vertices
        let tm = transition_matrix_b(wadj, n, m);

        print_matrix(&tm);

        let mut best_cluster = Vec::new();
        let mut best_weight = f64::NAN;


        // Find the B_cluster with minimal cost
        for b in 0..m {
            if assigned[b] {
                continue;
            }

            println!("Step 1: compute order of {b}");
            let mut d = 0.;
            for (_,w) in wadj[b].iter() {
                d += w;
            }
            let order = 
                if d == 0. {
                    vec![(b,1.)]
                } else {
                    compute_order(&wadj, m, b, &tm, markov_power)
                };
            println!("Step 2: compute best cost");
            let (w, b_cluster) = best(m, order, &wadj, cost_coef);
            if best_weight.is_nan() || w < best_weight {
                best_weight = w;
                best_cluster = b_cluster.clone();
            }
        }


        if best_cluster.len() == 0 { // This case should not happen
            break
        }
        
        for &b in best_cluster.iter() {
            assigned[b] = true;
        }
        b_clusters.push(best_cluster.clone());

        println!("Step 3: apply {best_cluster:?}");
        let a_cluster = apply_operations(n, m, best_cluster, wadj); 
        // todo!("modify nb_...");

        a_clusters.push(a_cluster);
    }

    // After having clustered every B vertices, it is possible that there remains unclustered A vertices
    let unclustered_a = compute_unclustered_a(wadj, n, m, &a_clusters);

    if unclustered_a.len() > 0 {
        println!("Unclustered A: {unclustered_a:?}");
        a_clusters.push(unclustered_a);
        b_clusters.push(vec![]);
    }

    // Print results
    println!("# Hyperparameters");
    println!("- cost coef: {cost_coef}");
    println!("- split coef: {split_coef}");
    println!("- markov power: {markov_power}");
    println!("
# Results
- eta: 
- nb biclusters: 
- nb operations: {nb_operations}
");

    // Check integrity
    check_integrity(&a_clusters, &b_clusters, n, m);

    compute_clusters(&a_clusters, &b_clusters, n, m)

}      




//     printv(verbose, "# Biclustering results")
//     printv(verbose, f"- eta: {(nb_additions+nb_deletions)/(init_n*init_m):.3}")
//     printv(verbose, f"- nb biclusters: {len(Bclusters)}")

//     printv(verbose, f"- nb operations: {nb_operations} (del: {nb_deletions}, add: {nb_additions}, spl: {nb_splits})")




fn delete_edge(wadj: &mut Vec<HashMap<usize, f64>>, a: usize, b: usize, m: usize) {
    wadj[a+m].remove(&b);
    wadj[b].remove(&a);
}

fn add_edge(wadj: &mut Vec<HashMap<usize, f64>>, a: usize, b: usize, weight: f64, m: usize) {
    wadj[a+m].insert(b, weight);
    wadj[b].insert(a, weight);
}




fn apply_operations(n: usize, m: usize, b_cluster: Vec<usize>, wadj: &mut Vec<HashMap<usize, f64>>) -> Vec<usize> {
    let mut a_cluster = vec![];
    let mut nb_operations = 0.;
    let mut nb_deletions = 0.;
    let mut nb_splits = 0.;
    let mut nb_additions = 0.;
    let blen = b_cluster.len() as f64;

    for a in 0..n {
        let mut indegree = 0.;
        for (i,w) in wadj[a+m].iter() {
            if b_cluster.contains(&i) {
                indegree += w;
            }
        } 
        if indegree > blen*0.5 {
            a_cluster.push(a);
            nb_operations += blen - indegree;
            nb_additions += blen - indegree;

            for &b in b_cluster.iter() {
                if wadj[b].contains_key(&a) {
                    delete_edge(wadj, a, b, m);
                }   
            }

            // Compute the number of neighbors of a outside of X
            let mut out_degree = 0.;
            for (i,w) in wadj[a+m].iter() {
                if b_cluster.contains(&i) == false {
                    out_degree += w;
                }
            }

            // Delete case:
            if out_degree <= 1. {
                nb_operations += out_degree;
                nb_deletions += out_degree;
                for i in 0..m {
                    if b_cluster.contains(&i) == false && wadj[i].contains_key(&a) {
                        delete_edge(wadj, a, i, m);
                    }
                }
            } else { // Split case
                nb_operations += 1.;
                nb_splits += 1.;
                for b in b_cluster.iter() {
                    if wadj[*b].contains_key(&a) {
                        delete_edge(wadj, a, *b, m)
                    }
                }
            }



        }
        else {
            nb_operations += indegree;
            nb_deletions += indegree;
            for i in b_cluster.iter() {
                if wadj[*i].contains_key(&a) {
                    delete_edge(wadj, a,*i, m);
                }
            }
        }
    }

    a_cluster
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
        clusters.push(cluster);
    }
    clusters
}



fn check_integrity(a_clusters: &Vec<Vec<usize>>, b_clusters: &Vec<Vec<usize>>, n: usize, m: usize) -> bool {
    let mut sum = 0;
    for b_cluster in b_clusters {
        sum += b_cluster.len();
    }
    if sum != m {
        println!("ERROR: sum={sum} of B clusters size is different than m={m}");
        return false;
    }


    // Check B clusters are disjoint

    // Check that all A vertices are clustered

    true
}

    

//     # Check disjointness
//     for k in range(len(Bclusters)):
//         for k2 in range(len(Bclusters)):
//             if k < k2:
//                 for i in Bclusters[k]:
//                     if i in Bclusters[k2]:
//                         print(f"bug {i} in Bclusters {k} and {k2}")
//     if m != Bsum:
//         raise f"m: {m} Bsum: {Bsum}"


//     memberships = memberships_from_clusters_list(Aclusters, Bclusters, init_n, init_m)
    


//     return 0




use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::common::print_matrix;

pub fn load_wadj_from_csv(file_path: &str, del: &str) -> (Vec<HashMap<usize, f64>>, usize, usize, HashMap<String, usize>, HashMap<String, usize>) {

    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);

    let mut n = 0;
    let mut m = 0;
    let mut node_map_a = HashMap::<String, usize>::new();
    let mut node_map_b = HashMap::<String, usize>::new();
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

            let node1 = String::from(values[0]);
            let node2 = String::from(values[1]);

            // Add both nodes to the maps if not already present
            if !node_map_a.contains_key(&node1) {
                node_map_a.insert(node1.clone(), n);
                n += 1;
            }
            if !node_map_b.contains_key(&node2) {
                node_map_b.insert(node2.clone(), m);
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

    let mut wadj = vec![HashMap::new(); n+m];
    for &(a,b,w) in edges.iter() {
        wadj[a+m].insert(b, w);
        wadj[b].insert(a,w);
    }



    (wadj, n, m, node_map_a, node_map_b)
}


pub fn print_biclusters_stats(biclusters: &Vec<Vec<usize>>, n: usize, m: usize) {
    println!("Number of biclusters: {}", biclusters.len())
}



pub fn print_wadj_stats(wadj: &Vec<HashMap<usize, f64>>, n: usize, m: usize) {
    println!("# Data statistics:");
    println!("n={n}");
    println!("m={m}");

    let mut minw = f64::NAN;
    let mut maxw = f64::NAN;
    let mut ne = 0;
    let mut density = 0.;
    for b in 0..m {
        ne += wadj[b].len();
        for (_, &w) in wadj[b].iter() {
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
    println!("Number of edges: {ne}");
    println!("Edge density: {}", (ne as f64)/((n*m) as f64));

    println!("Min weight: {minw}");
    println!("Max weight: {maxw}");
    println!("Density: {density:.2}")

}