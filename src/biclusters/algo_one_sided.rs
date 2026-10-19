use core::f64;
use ndarray::Array2;
use rand::rng;
use rand::rngs::ThreadRng;
use rand::seq::SliceRandom;
use std::collections::{HashMap, HashSet};
use std::time::Instant;

use super::common::transition_matrix_b;
use super::weighted_biadj::WeightedBiAdjacency;
use crate::biclusters::biclust::Biclust;
use crate::biclusters::common::AlgoStats;
use crate::common::{compute_statio_distrib_by_iter, print_matrix, progress_bar};


fn compute_order_local(
    wadj: &WeightedBiAdjacency, 
    center_col: usize,
    iters: usize, 
    tol: f64,
) -> Vec<(usize, f64)> {
    // Boule S, b en position 0
    let mut local: HashMap<usize, usize> = HashMap::from([(center_col, 0)]);
    let mut cols = vec![center_col];
    for (&row, _) in wadj.col_adj[center_col].iter() {
        for (&c, _) in wadj.row_adj[row].iter() {
            local.entry(c).or_insert_with(|| { cols.push(c); cols.len() - 1 });
        }
    }
    let d = cols.len();

    // Demi-pas précalculés : fwd[j] = [(ligne, W/d_B(j))], back[ligne] = [(j dans S, W/d_A(a))]
    let mut row_local: HashMap<usize, usize> = HashMap::new();
    let mut back: Vec<Vec<(usize, f64)>> = vec![];
    let mut fwd: Vec<Vec<(usize, f64)>> = vec![vec![]; d];
    for (j, &c) in cols.iter().enumerate() {
        for (&row, &w) in wadj.col_adj[c].iter() {
            let r = *row_local.entry(row).or_insert_with(|| {
                back.push(wadj.row_adj[row].iter()
                    .filter_map(|(c2, &w2)| local.get(c2).map(|&j2| (j2, w2 / wadj.row_wdeg(row))))
                    .collect());
                back.len() - 1
            });
            fwd[j].push((r, w / wadj.col_wdeg(c)));
        }
    }

    // Itération de puissance creuse, la masse sortante retourne en b
    let (mut pi, mut next, mut h) = (vec![0.; d], vec![0.; d], vec![0.; back.len()]);
    pi[0] = 1.;
    for _ in 0..iters {
        h.fill(0.); next.fill(0.);
        for j in 0..d { if pi[j] > 0. { for &(r, p) in &fwd[j] { h[r] += pi[j] * p; } } }
        for (r, &hr) in h.iter().enumerate() { if hr > 0. { for &(j2, p) in &back[r] { next[j2] += hr * p; } } }
        let stay: f64 = next.iter().sum();
        next[0] += 1. - stay;
        let diff: f64 = pi.iter().zip(&next).map(|(x, y)| (x - y).abs()).sum();
        std::mem::swap(&mut pi, &mut next);
        if diff < tol { break; }
    }
    let mut order: Vec<_> = cols.into_iter().zip(pi).collect();
    order.sort_by(|x, y| y.1.partial_cmp(&x.1).unwrap());
    order
}




/// The center vertex is supposed to be at index 0 in neighbors (so it is the closed neighborhood)
///
fn centered_transition_matrix(tm: &Array2<f64>, neighbors: &Vec<usize>) -> Array2<f64> {
    let d = neighbors.len();
    let mut ctm = Array2::zeros((d, d));

    for ni in 0..d {
        let i = neighbors[ni];
        let mut s = 0.;
        for nj in 0..d {
            let j = neighbors[nj];
            ctm[[ni, nj]] = tm[[i, j]];
            s += ctm[[ni, nj]];
        }
        ctm[[ni, 0]] += 1. - s;
    }

    ctm
}

fn compute_order(
    wadj: &WeightedBiAdjacency,
    m: usize,
    vertex: usize,
    tm_common: &Array2<f64>,
    markov_power: usize,
    verbose: usize,
) -> Vec<(usize, f64)> {
    // Search the B vertices which have a common neighbor with vertex
    // neighbors[0] = vertex
    let mut neighbors: Vec<usize> = vec![vertex];
    for b in 0..m {
        if b == vertex {
            continue;
        }
        for (x, _) in wadj.col_adj[b].iter() {
            let mut found = false;
            for (y, _) in wadj.col_adj[vertex].iter() {
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

    let ctm = centered_transition_matrix(tm_common, &neighbors)
        .t()
        .into_owned();
    if verbose >= 2 {
        println!("Centered to {vertex} transition matrix:");
        print_matrix(&ctm);
    }

    let v_result = compute_statio_distrib_by_iter(&ctm, markov_power, verbose);

    // By exponentiation
    // let d = neighbors.len();
    // for _ in 0..markov_power {
    //     ctm = ctm.dot(&ctm);
    // }
    // let mut v = Array2::zeros((d, 1));
    // v[[0, 0]] = 1.0;
    // let v_result = ctm.dot(&v);

    // Order subset by decreasing probability
    let mut order: Vec<(usize, f64)> = neighbors
        .iter()
        .enumerate()
        .map(|(ni, &i)| (i, v_result[[ni, 0]]))
        .collect();

    order.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    order
}

#[derive(PartialEq, Eq)]
enum State {
    Negative,
    Positive,
}


///
/// cost_coef in [0,1]
fn best(
    order: Vec<(usize, f64)>,
    wadj: &WeightedBiAdjacency,
    cost_coef: f64,
    split_threshold: f64,
) -> (f64, Vec<usize>) {
    let n = wadj.get_n();
    let mut best_cost = f64::NAN;
    let mut best_cluster = Vec::new();
    let mut cols_cluster = Vec::new();

    let mut s = 0.;
    let mut c = 0.;
    let mut state = HashMap::new();
    let mut indegree = vec![0.; n];
    let mut outdegree = vec![0.; n];

    for (i, proba_i) in order {
        if proba_i == 0. {
            break;
        }
        cols_cluster.push(i);
        s += 1.;

        for (x, st) in state.iter_mut() {
            if wadj.has_edge(*x, i) == false {
                if *st == State::Positive {
                    if indegree[*x] <= s * 0.5 {
                        *st = State::Negative;
                        // state.insert(*x, State::Negative);
                        c += indegree[*x];
                        c -= (s - 1.) - indegree[*x];
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
        for (&x, &w) in wadj.col_adj[i].iter() {
            if state.contains_key(&x) == false {
                indegree[x] = w;
                outdegree[x] = wadj.row_wdeg(x) - w;
                if indegree[x] > s * 0.5 {
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
            } else if state[&x] == State::Positive {
                indegree[x] += w;
                outdegree[x] -= w;

                if indegree[x] > s * 0.5 {
                    // Keep connected
                    c += 1. - w;
                    if outdegree[x] + w <= split_threshold {
                        c -= w;
                    } else if outdegree[x] <= split_threshold {
                        c += outdegree[x] - 1.;
                    }
                } else {
                    // Unconnect x
                    state.insert(x, State::Negative);
                    c -= (s - 1.) - (indegree[x] - w);
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
                if indegree[x] > s * 0.5 {
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

        let cc = compute_nb_operations(&cols_cluster, wadj, split_threshold, 0);

        println!("{s}, {cols_cluster:?} cost: {c}, ope {cc:?}");
        if c != (cc.1 + cc.2 + cc.3) {
            panic!("{c} {cc:?}");
        }
        

        let cost = c * f64::powf(s, -cost_coef);
        if best_cost.is_nan() || cost < best_cost {
            best_cost = cost;
            best_cluster = cols_cluster.clone();
        }
    }

    (best_cost, best_cluster)
}


/// return (nb_rows, nb_deletions, nb_additions, nb_splits)
fn compute_nb_operations (
    cols_cluster: &Vec<usize>,
    wadj: & WeightedBiAdjacency,
    split_threshold: f64,
    verbose: usize,
) -> (f64, f64, f64, f64) {
    let mut nb_rows = 0.;
    let mut nb_deletions = 0.;
    let mut nb_splits = 0.;
    let mut nb_additions = 0.;
    let cluster_size = cols_cluster.len() as f64;

    let mut in_cols = vec![false; wadj.get_m()];
    for &col in cols_cluster {
        in_cols[col] = true;
    }

    for row in 0..wadj.get_n() {
        // the indegree is the weights of the edges with cols_cluster
        // the out_degree is the sum of the weight of the edges outside cols_cluster
        let (indegree, out_degree) = wadj.get_row_in_out_degrees(row, &in_cols);

        // row is accepted in the bicluster
        if indegree > cluster_size * 0.5 {
            nb_rows += 1.;
            nb_additions += cluster_size - indegree;

            // Delete case
            if out_degree <= split_threshold {
                nb_deletions += out_degree;
            }
            // Split case
            else {
                nb_splits += 1.;
            }
        }
        // row is not accepted
        else {
            nb_deletions += indegree;
        }
    }

    (nb_rows, nb_deletions, nb_additions, nb_splits)
}



fn compute_cost(
    cols_cluster: &Vec<usize>,
    wadj: & WeightedBiAdjacency,
    split_threshold: f64,
    size_sensitivity: f64,
    verbose: usize,
) -> f64 {
    let (nb_rows, add, del, spl) = compute_nb_operations(cols_cluster, wadj, split_threshold, verbose);
    let c = cols_cluster.len() as f64;
    let cost = (add)/c + del+ spl/ (nb_rows+1.);
    // println!("{cost:.2} {} {nb_rows:.0} {add:.0} {del:.0} {spl:.0}", cols_cluster.len());
    cost
}

fn compute_old_cost(
    cols_cluster: &Vec<usize>,
    wadj: & WeightedBiAdjacency,
    split_threshold: f64,
    size_sensitivity: f64,
    verbose: usize,
) -> f64 {

    let nb_cols_cluster = cols_cluster.len() as f64;
    let (nb_rows, add, del, spl) = compute_nb_operations(cols_cluster, wadj, split_threshold, verbose);

     (add + del + spl) * f64::powf(nb_cols_cluster, -size_sensitivity)
}



fn best_sure(
    order: &Vec<(usize, f64)>,
    wadj: &WeightedBiAdjacency,
    cost_coef: f64,
    split_threshold: f64,
    verbose: usize
) -> (f64, Vec<usize>) {
    let mut best_cost = f64::NAN;
    let mut best_cluster = Vec::new();
    let mut cols_cluster = Vec::new();

    for &(i,proba_i) in order.iter() {
        if proba_i == 0. {
            break;
        }
        cols_cluster.push(i);
        let cost = compute_cost(&cols_cluster, wadj, split_threshold, cost_coef, verbose);
        
        if  best_cost.is_nan() || cost < best_cost {
            best_cluster = cols_cluster.clone();
            best_cost = cost;
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



/// Pick a random sample of size k in (0..n-1) such that no indices are in assigned.
/// If k > the number of unassigned vertices, then retain all the unassigned vertices
fn pick_unassigned_sample(
    rng: &mut ThreadRng,
    k: usize,
    unassigned: &HashSet<usize>,
) -> Vec<usize> {
    let mut all_numbers: Vec<usize> = vec![];
    for v in unassigned {
        all_numbers.push(*v);
    }
    // all_numbers.retain(|&x| !assigned.contains(&x));

    if all_numbers.len() < k {
        return all_numbers;
    }

    all_numbers.shuffle(rng);
    let mut result = vec![];
    for i in 0..k {
        result.push(all_numbers[i]);
    }
    result
}

/// Main one sided biclustering function
///
/// Partition columns
///
/// Rows may overlap
pub fn bicluster_one_sided(
    wadj: &mut WeightedBiAdjacency,
    cost_coef: f64,
    split_threshold: f64,
    markov_power: usize,
    sample_size: usize,
    verbose: usize,
) -> (Biclust, AlgoStats) {
    let min_error = wadj.compute_min_error();
    let mut rng = rng();

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

    let start_instant: Instant = Instant::now();


    // While there exists some B vertices to cluster
    loop {
        if verbose >= 1 {
            println!("---");
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
            if assigned[b] == false && wadj.col_nb_neighbors(b) == 0 {
                assigned[b] = true;
                nb_assigned += 1;
                isolated_b_vertices.push(b);
                has_isolated_b_vertices = true;
            }
        }
        if has_isolated_b_vertices {
            continue;
        }

        // Search min/max degree
        let mut mindeg = 100000;
        let mut minv = None;
        let mut maxdeg = 0;
        let mut maxv = None;

        for v in 0..m {
            if assigned[v] {
                continue;
            }
            let degree = wadj.col_nb_neighbors(v);
            if degree > maxdeg {
                maxdeg = degree;
                maxv = Some(v);
            }
            if degree < mindeg {
                mindeg = degree;
                minv = Some(v);
            }
        }
        if verbose >= 1 {
            println!("vertex {minv:?} of degree {mindeg}, vertex {maxv:?} of degree {maxdeg}");
        }

        // Compute the transition matrix between B vertices
        let tm = transition_matrix_b(wadj, n, m);

        if verbose >= 2 {
            println!("Transition matrix:");
            print_matrix(&tm);
        }

        let mut best_cluster = Vec::new();
        let mut best_cost = f64::NAN;

        let mut unassigned: HashSet<usize> = (0..0).collect();
        for b in 0..m {
            if assigned[b] == false {
                unassigned.insert(b);
            }
        }

        // Pick a sample
        let mut sample = pick_unassigned_sample(&mut rng, sample_size, &unassigned);


        // Add minv and maxv if they are not in the sample
        // if let Some(minv) = minv {
        //     if sample.contains(&minv) == false {
        //         sample.push(minv);
        //     }
        // }
        if let Some(maxv) = maxv {
            if sample.contains(&maxv) == false {
                sample.push(maxv);
            }
        }

        let mut best_col = 0;

        // Find the B_cluster with minimal cost
        for col in sample {

            println!("center: {col}");

            // Compute local order of col
            let d = wadj.col_wdeg(col);
            let order = if d == 0. {
                vec![(col, 1.)]
            } else {
                compute_order_local(&wadj, col,  markov_power, 0.01)
            };

            println!("order computed");
            
            if verbose >= 2 {
                println!("Step 1: compute order of {col}");
                println!("{order:?}");
            }



            // Compute best cost
            // let (cost, b_cluster) = best( order, &wadj, cost_coef, split_threshold);
            let (cost, b_cluster) = best_sure(&order, &wadj, cost_coef, split_threshold, verbose);

            println!(
                "Col: {col} deg: {} B-neighbors: {} best_cost: {cost:.2}",
                wadj.col_nb_neighbors(col),
                order.len()
            );


            if best_cost.is_nan() || cost < best_cost {
                best_cost = cost;
                best_cluster = b_cluster.clone();
                best_col = col;
            }
            if verbose >= 2 {
                println!("Step 2: best_cost: {best_cost}, cluster: {best_cluster:?}");
            }
        }

        if best_cluster.len() == 0 {
            // This case should not happen
            panic!("best_cluster is empty");
        }

        for &b in best_cluster.iter() {
            if assigned[b] {
                panic!("col {b} already assigned");
            }
            assigned[b] = true;
            nb_assigned += 1;
        }
        b_clusters.push(best_cluster.clone());

        if verbose >= 0 {
            best_cluster.sort();
            println!("--- New Bicluster n° {}", b_clusters.len());
            println!("Discoverd from {best_col}");
            // println!("Cols indices:  {best_cluster:?}");
            // let mut cols_labels = vec![];
            // for i in 0..best_cluster.len() {
            //     cols_labels.push(wadj.get_label(n + best_cluster[i]));
            // }
            // println!("Cols labels: [{}]", cols_labels.join(" "));
            println!("Nb cols:  {}", best_cluster.len());
        }

        let (a_cluster, del, add, spl) =
            apply_operations(n, m, best_cluster, wadj, split_threshold, verbose);
        nb_deletions += del;
        nb_additions += add;
        nb_splits += spl;
        nb_operations += del + add + spl;

        if verbose >= 0 {
            // println!("Rows indices: {a_cluster:?}");
            // let mut rows_labels = vec![];
            // for i in 0..a_cluster.len() {
            //     rows_labels.push(wadj.get_label(a_cluster[i]));
            // }
            // println!("Rows labels: [{}]", rows_labels.join(" "));
            println!("Nb rows: {}", a_cluster.len());
            println!("Deletions: {del}");
            println!("Additions: {add}");
            println!("Splits: {spl}");
            println!("Cost: {best_cost:.3}");
        }

        a_clusters.push(a_cluster);
    }

    // Add isolated B vertices
    if isolated_b_vertices.len() > 0 {
        b_clusters.push(isolated_b_vertices.clone());
        a_clusters.push(vec![]);
    }

    // After having clustered every B vertices, it is possible that there remains unclustered A vertices
    let unclustered_a = compute_unclustered_a(n, &a_clusters);

    if unclustered_a.len() > 0 {
        a_clusters.push(unclustered_a.clone());
        b_clusters.push(vec![]);
    }

    let eta = (nb_additions + nb_deletions) / ((n * m) as f64);
    let eta = eta - min_error;

    let bicluster_stats = AlgoStats {
        adjusted_error: eta,
        nb_operations: nb_operations,
        nb_splits: nb_splits,
        nb_additions: nb_additions,
        nb_deletions: nb_deletions,
    };

    // Check integrity
    check_integrity(&a_clusters, &b_clusters, n, m);

    // compute_clusters(&a_clusters, &b_clusters, n, m)
    (
        Biclust::from_separate_biclusters(n, m, &a_clusters, &b_clusters),
        bicluster_stats,
    )
}


fn apply_operations(
    n: usize,
    m: usize,
    b_cluster: Vec<usize>,
    wadj: &mut WeightedBiAdjacency,
    split_threshold: f64,
    verbose: usize,
) -> (Vec<usize>, f64, f64, f64) {
    let mut a_cluster = vec![];
    let mut nb_deletions = 0.;
    let mut nb_splits = 0.;
    let mut nb_additions = 0.;
    let blen = b_cluster.len() as f64;

    for a in 0..n {
        let mut indegree = 0.;
        for (&col, w) in wadj.row_adj[a].iter() {
            if b_cluster.contains(&col) {
                indegree += w;
                // println!("{col}");
                // println!("{}", wadj.has_edgee(a, col));
            }
        }

        if indegree > blen * 0.5 {
            // println!("{a} indegree: {indegree}");
            a_cluster.push(a);
            nb_additions += blen - indegree;

            if verbose >= 2 {
                for &col in b_cluster.iter() {
                    if wadj.has_edge(a, col) == false {
                        println!("add {a} {col}");
                    }
                }
            }

            // Compute the number of neighbors of row a outside b_cluster
            let mut out_degree = 0.;
            for (col, w) in wadj.row_adj[a].iter() {
                if b_cluster.contains(&col) == false {
                    out_degree += w;
                }
            }

            // Delete case
            if out_degree <= split_threshold {
                nb_deletions += out_degree;
                for b in 0..m {
                    if b_cluster.contains(&b) == false && wadj.has_edge(a, b) {
                        if verbose >= 2 {
                            println!("delete {a} {b}")
                        }
                        wadj.delete_edge(a, b);
                    }
                }
                for &b in b_cluster.iter() {
                    // I think it is useless
                    if wadj.has_edge(a, b) {
                        wadj.delete_edge(a, b);
                    }
                }
            }
            // Split case
            else {
                nb_splits += 1.;
                if verbose >= 1 {
                    println!("split {a} {}", wadj.row_label(a))
                }
                for &b in b_cluster.iter() {
                    if wadj.has_edge(a, b) {
                        // println!("spl del {a} {b}");
                        wadj.delete_edge(a, b);
                    }
                }
            }
        } else {
            nb_deletions += indegree;
            for &b in b_cluster.iter() {
                if wadj.has_edge(a, b) {
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

// fn compute_clusters(a_clusters: &Vec<Vec<usize>>, b_clusters: &Vec<Vec<usize>>, n: usize, m: usize) -> Vec<Vec<usize>> {
//     let mut clusters = vec![];
//     for i in 0..b_clusters.len() {
//         let mut cluster = vec![];
//         for &b in b_clusters[i].iter() {
//             cluster.push(b+n);
//         }
//         for &a in a_clusters[i].iter() {
//             cluster.push(a);
//         }
//         cluster.sort();
//         clusters.push(cluster);
//     }
//     clusters
// }

/// Check if the b_clusters form a partition of {0, ..., m-1}
/// and that a_clusters are covering {0, ..., n-1}.
fn check_integrity(
    a_clusters: &Vec<Vec<usize>>,
    b_clusters: &Vec<Vec<usize>>,
    n: usize,
    m: usize,
) -> bool {
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

fn biclusters_contains_subset(biclusters: &Vec<Vec<usize>>, subset: &Vec<usize>) -> bool {
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

pub fn compute_edition_diff(
    biclusters: &Vec<Vec<usize>>,
    wadj: &Vec<HashMap<usize, f64>>,
    n: usize,
    m: usize,
) -> f64 {
    let mut r = 0;

    for a in 0..n {
        for b in 0..m {
            if wadj[b].contains_key(&a) != biclusters_contains_subset(biclusters, &vec![a, b + n]) {
                println!("{a} {} ", b + n);
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
        if biclusters_contains_subset(biclusters, &vec![b + n]) == false {
            rb += 1;
        }
    }
    (ra, rb)
}





#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_compute_nb_operations() {
        let input = "\
+ a b c
X 1 0 0
Y 1 0 1
Z 1 1 1
";

        let reader = Cursor::new(input);

        let wadj = WeightedBiAdjacency::parse_wadj_from_reader(
            reader,
            false,
            " ",
            true,
        )
        .unwrap();

        let cols_cluster = vec![0];

        // Threshold 0.5
        let (nb_rows, nb_deletions, nb_additions, nb_splits) =
            compute_nb_operations(
                &cols_cluster,
                &wadj,
                0.5,
                0,
            );

        assert_eq!(nb_rows, 3.);

        assert!((nb_deletions - 0.0).abs() < 1e-9);
        assert!((nb_additions - 0.0).abs() < 1e-9);
        assert!((nb_splits - 2.0).abs() < 1e-9);

        // Threshold 1.0
        let (nb_rows, nb_deletions, nb_additions, nb_splits) =
            compute_nb_operations(
                &cols_cluster,
                &wadj,
                1.0,
                0,
            );

        assert_eq!(nb_rows, 3.);

        assert!((nb_deletions - 1.0).abs() < 1e-9);
        assert!((nb_additions - 0.0).abs() < 1e-9);
        assert!((nb_splits - 1.0).abs() < 1e-9);


        // Threshold 2.0
        let (nb_rows, nb_deletions, nb_additions, nb_splits) =
            compute_nb_operations(
                &cols_cluster,
                &wadj,
                2.5,
                0,
            );

        assert_eq!(nb_rows, 3.);

        assert!((nb_deletions - 3.0).abs() < 1e-9);
        assert!((nb_additions - 0.0).abs() < 1e-9);
        assert!((nb_splits - 0.0).abs() < 1e-9);
    }

   
}