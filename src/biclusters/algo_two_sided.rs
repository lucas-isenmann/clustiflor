use core::f64;
use std::collections::HashSet;
use std::time::Instant;

use ndarray::Array2;



use crate::biclusters::biclust::Biclust;
use crate::common::{print_matrix, progress_bar};

use super::algo::AlgoStats;
use super::biclustering::Biclustering;
use super::common::{rows_transition_matrix, transition_matrix_b};
use super::weighted_biadj::WeightedBiAdjacency;







// /// centered on cols[0]
// fn centered_cols_transition_matrix(  cols_tm: &Array2<f64>, cols: &Vec<usize>) -> Array2<f64> {
//     let d = cols.len();
//     let mut ctm = Array2::zeros((d,d));

//     for ni in 0..d {
//         let i = cols[ni];
//         let mut s = 0.;
//         for nj in 0..d {
//             let j = cols[nj];
//             ctm[[ni,nj]] = cols_tm[[i,j]];
//             s += ctm[[ni,nj]];
//         }
//         ctm[[ni,0]] += 1.-s;
//     }

//     ctm
// }

/// Centered on cols1
fn centered_transition_matrix2(  tm: &Array2<f64>, set1: &Vec<usize>, set2: &Vec<usize>) -> Array2<f64> {
    let d1 = set1.len();
    let d2 = set2.len();
    let mut ctm = Array2::zeros((d1+d2,d1+d2));

    for ni in 0..d1 {
        let i = set1[ni];
        for nj in 0..d1 {
            let j = set1[nj];
            ctm[[ni,nj]] = tm[[i,j]];
        }
        for nj in 0..d2 {
            let j = set2[nj];
            ctm[[ni,d1+nj]] = tm[[i,j]];
        }
    }
    for ni in 0..d2 {
        let i = set2[ni];
        let mut s =   0.;
        for nj in 0..d1 {
            let j = set1[nj];
            ctm[[d1+ni, nj]] = tm[[i,j]];
            s += ctm[[d1+ni, nj]];
        }
        for nj in 0..d2 {
            let j = set2[nj];
            ctm[[d1+ni, d1+nj]] = tm[[i,j]];
            s += ctm[[d1+ni, d1+nj]];
        }

        let p = (1.-s)/ ( d1 as f64);
        for nj in 0..d1 {
            ctm[[ni,nj]] += p;
        }
    }
    ctm
}


fn centered_ordering(tm: &Array2<f64>, set1: Vec<usize>, set2: Vec<usize>, markov_power: usize) -> Vec<(usize, f64)>{
    let d1 = set1.len();
    let d3 = set2.len();
    let mut ctm = centered_transition_matrix2(tm, &set1, &set2).t().into_owned();
    for _ in 0..markov_power {
        ctm = ctm.dot(&ctm);
    }

    let mut v = Array2::zeros((d1+d3, 1));
    for i in 0..d1 {
        v[[i, 0]] = 1.0 / (d1 as f64);
    }
    
    let v_result = ctm.dot(&v);
    let set = [set1, set2].concat();

    // Order subset by decreasing probability
    let mut ordered_set: Vec<(usize, f64)> = set.iter().enumerate()
        .map(|(ni, &i)| (i, v_result[[ni, 0]]))
        .collect();
    ordered_set.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    ordered_set
}



fn neighbors_3_col(wadj: &WeightedBiAdjacency, m: usize, col: usize) -> (HashSet<usize>, HashSet<usize>, HashSet<usize>){
    let mut row_1neighbors = HashSet::new();
    let mut col_2neighbors = HashSet::new();
    let mut row_3neighbors = HashSet::new();

    for (row,_) in wadj.iter(col) {
        row_1neighbors.insert(*row);
        for (col2, _) in wadj.iter(row+m){
            if *col2 != col && col_2neighbors.contains(col2) == false {
                col_2neighbors.insert(*col2);
            }
        }
    }

    for col2 in col_2neighbors.iter() {
        for (row3, _) in wadj.iter(*col2) {
            if row_1neighbors.contains(row3) == false {
                row_3neighbors.insert(*row3);
            }
        }
    }

    (row_1neighbors, col_2neighbors, row_3neighbors)
}


fn neighbors_3_row(wadj: &WeightedBiAdjacency, m: usize, row: usize) -> (HashSet<usize>, HashSet<usize>, HashSet<usize>){
    let mut col_1neighbors = HashSet::new();
    let mut row_2neighbors = HashSet::new();
    let mut col_3neighbors = HashSet::new();

    for (col,_) in wadj.iter(row+m) {
        col_1neighbors.insert(*col);
        for (row2, _) in wadj.iter(*col){
            if *row2 != row && row_2neighbors.contains(row2) == false {
                row_2neighbors.insert(*row2);
            }
        }
    }

    for row2 in row_2neighbors.iter() {
        for (col3, _) in wadj.iter(*row2+m) {
            if col_1neighbors.contains(col3) == false {
                col_3neighbors.insert(*col3);
            }
        }
    }

    (col_1neighbors, row_2neighbors, col_3neighbors)
}


fn compute_bicluster_col(wadj: &WeightedBiAdjacency, m: usize, col: usize, rows_tm: &Array2<f64>, cols_tm: &Array2<f64>, markov_power: usize, verbose: usize)
 -> (Vec<(usize,f64)>, Vec<(usize, f64)>) {
    
    let (row_1neighbors, col_2neighbors, row_3neighbors) = neighbors_3_col(wadj, m, col);

    if verbose >= 3 {
        println!("   col {col}");
        println!("   adjacent rows1: {row_1neighbors:?}");
        println!("   adjacent cols2: {col_2neighbors:?}");
        println!("   adjacent rows3: {row_3neighbors:?}");
    }

    // Columns
    let c0: Vec<usize> = vec![col];
    let c2: Vec<usize> = col_2neighbors.iter().map(|x| *x).collect();
    let cols_order = centered_ordering(cols_tm, c0, c2, markov_power);

    // Rows
    let r1: Vec<usize> = row_1neighbors.iter().map(|x| *x).collect();
    let r3: Vec<usize> = row_3neighbors.iter().map(|x| *x).collect();
    let rows_order = centered_ordering(rows_tm, r1, r3, markov_power);

    (rows_order, cols_order)
}




fn bicluster_size(n: usize, m: usize) -> f64{
    let n = n as f64;
    let m = m as f64;
    n + m 
}


fn best_cluster_col(wadj: &WeightedBiAdjacency, n: usize, m: usize, col: usize, rows_tm: &Array2<f64>, cols_tm: &Array2<f64>, markov_power: usize, split_threshold: f64, cost_coef: f64, verbose: usize)
 -> (Vec<usize>, Vec<usize>, f64) {
    
    let (row_1neighbors, col_2neighbors, row_3neighbors) = neighbors_3_col(wadj, m, col);

    if verbose >= 3 {
        println!("   col {col}");
        println!("   adjacent rows1: {row_1neighbors:?}");
        println!("   adjacent cols2: {col_2neighbors:?}");
        println!("   adjacent rows3: {row_3neighbors:?}");
    }

    // Columns
    let c0: Vec<usize> = vec![col];
    let c2: Vec<usize> = col_2neighbors.iter().map(|x| *x).collect();
    let cols_order = centered_ordering(cols_tm, c0, c2, markov_power);
    // println!("col: {col}, cols_order: {cols_order:?}");


    // Rows is the concatenation of row_1 and row_3
    let mut rows = vec![];
    for &row in row_1neighbors.iter() {
        rows.push(row);
    }
    for &row in row_3neighbors.iter() {
        rows.push(row);
    }
    // println!("rows1: {row_1neighbors:?}");
    // println!("rows3: {row_3neighbors:?}");

    let mut best_cost = f64::NAN;
    let mut best_cols_cluster = vec![];
    let mut best_rows_cluster = vec![];
    let mut cols_cluster = vec![];
    for i in 0..cols_order.len() {
        cols_cluster.push(cols_order[i].0);

        // println!("cols_cluster: {cols_cluster:?}");

        let mut cost = vec![0;n];
        for (i, &row) in rows.iter().enumerate() {
            let mut d = 0;
            for &col in cols_cluster.iter() {
                if wadj.has_edgee(row, col) {
                    d += 1;
                }
            }
            cost[row] = cols_cluster.len()-d; // edge additions
            if (wadj.row_degree(row)-d) as f64 > split_threshold {
                cost[row] += 1; // split
            } else {
                cost[row] += wadj.row_degree(row)-d;
            }
        }

        rows.sort_by(
            | row1, row2| { cost[*row1].cmp(&cost[*row2]) }
        );

        // println!("rows_cost: {cost:?}");
        // println!("sorted rows: {rows:?}");

        let mut rows_cost = 0;

        let mut rows_cluster = vec![];
        for (i, &row) in rows.iter().enumerate() {
            rows_cost += cost[row];
            rows_cluster.push(row);
            let mut cols_cost = 0;
            
            for &col in cols_cluster.iter() {
                let mut d = 0;
                for &row in rows_cluster.iter(){
                    if wadj.has_edgee(row, col){
                        d += 1;
                    }
                }
                if (wadj.col_degree(col) - d) as f64 > split_threshold {
                    cols_cost += 1;
                } else {
                    cols_cost += wadj.col_degree(col) - d;
                }
            }

            let size = bicluster_size(rows_cluster.len(), cols_cluster.len());
            // let size = (rows_cluster.len() + cols_cluster.len()) as f64;
            let total_cost = (rows_cost + cols_cost) as f64;
            let final_cost = total_cost * size.powf(-cost_coef);

            if best_cost.is_nan() || final_cost < best_cost {
                best_cost = final_cost;
                best_cols_cluster = cols_cluster.clone();
                best_rows_cluster = rows_cluster.clone();
            }
        }
    }

    // println!("best_cost: {best_cost:.2}");
    // println!("best_cols_cluster: {best_cols_cluster:?}");
    // println!("best_rows_cluster: {best_rows_cluster:?}");


    // panic!("salt");

    (best_rows_cluster, best_cols_cluster, best_cost)
}



fn best_cluster_row(wadj: &WeightedBiAdjacency, n: usize, m: usize, row: usize, rows_tm: &Array2<f64>, cols_tm: &Array2<f64>, markov_power: usize, split_threshold: f64, cost_coef: f64, verbose: usize)
 -> (Vec<usize>, Vec<usize>, f64) {
    
    let (col_1neighbors, row_2neighbors, col_3neighbors) = neighbors_3_row(wadj, m, row);


    // Rows
    let r0: Vec<usize> = vec![row];
    let r2: Vec<usize> = row_2neighbors.iter().map(|x| *x).collect();
    let rows_order = centered_ordering(rows_tm, r0, r2, markov_power);


    // Cols is the concatenation of row_1 and col_3
    let mut cols = vec![];
    for &col in col_1neighbors.iter() {
        cols.push(col);
    }
    for &col in col_3neighbors.iter() {
        cols.push(col);
    }
    // println!("cols1: {col_1neighbors:?}");
    // println!("cols3: {col_3neighbors:?}");

    let mut best_cost = f64::NAN;
    let mut best_cols_cluster = vec![];
    let mut best_rows_cluster = vec![];
    let mut rows_cluster = vec![];
    for i in 0..rows_order.len() {
        rows_cluster.push(rows_order[i].0);

        // println!("cols_cluster: {cols_cluster:?}");

        let mut cost = vec![0;m];
        for (i, &col) in cols.iter().enumerate() {
            let mut d = 0;
            for &row in rows_cluster.iter() {
                if wadj.has_edgee(row, col) {
                    d += 1;
                }
            }
            cost[col] = rows_cluster.len()-d; // edge additions
            if (wadj.col_degree(col)-d) as f64 > split_threshold {
                cost[col] += 1; // split
            } else {
                cost[col] += wadj.col_degree(col)-d;
            }
        }

        cols.sort_by(
            | x, y| { cost[*x].cmp(&cost[*y]) }
        );

        // println!("rows_cost: {cost:?}");
        // println!("sorted rows: {rows:?}");

        let mut cols_cost = 0;

        let mut cols_cluster = vec![];
        for (i, &col) in cols.iter().enumerate() {
            cols_cost += cost[col];
            cols_cluster.push(col);
            let mut rows_cost = 0;
            
            for &row in rows_cluster.iter() {
                let mut d = 0;
                for &col in cols_cluster.iter(){
                    if wadj.has_edgee(row, col){
                        d += 1;
                    }
                }
                if (wadj.row_degree(row) - d) as f64 > split_threshold {
                    rows_cost += 1;
                } else {
                    rows_cost += wadj.row_degree(row) - d;
                }
            }

            // let size = (cols_cluster.len() + rows_cluster.len()) as f64;
            let size = bicluster_size(rows_cluster.len(), cols_cluster.len());
            let total_cost = (rows_cost + cols_cost) as f64;
            let final_cost = total_cost * size.powf(-cost_coef);

            if best_cost.is_nan() || final_cost < best_cost {
                best_cost = final_cost;
                best_cols_cluster = cols_cluster.clone();
                best_rows_cluster = rows_cluster.clone();
            }
        }
    }

    // println!("best_cost: {best_cost:.2}");
    // println!("best_cols_cluster: {best_cols_cluster:?}");
    // println!("best_rows_cluster: {best_rows_cluster:?}");


    // panic!("salt");

    (best_rows_cluster, best_cols_cluster, best_cost)
}




fn compute_bicluster_row(wadj: &WeightedBiAdjacency, m: usize, row: usize, rows_tm: &Array2<f64>, cols_tm: &Array2<f64>, markov_power: usize, verbose: usize)
 -> (Vec<(usize,f64)>, Vec<(usize, f64)>) {
    
    let (col_1neighbors, row_2neighbors, col_3neighbors) = neighbors_3_row(wadj, m, row);

    if verbose >= 3 {
        println!("   row {row}");
        println!("   adjacent rows1: {col_1neighbors:?}");
        println!("   adjacent cols2: {row_2neighbors:?}");
        println!("   adjacent rows3: {col_3neighbors:?}");
    }

    // Rows
    let c0: Vec<usize> = vec![row];
    let c2: Vec<usize> = row_2neighbors.iter().map(|x| *x).collect();
    let rows_order = centered_ordering(rows_tm, c0, c2, markov_power);

    // Columns
    let r1: Vec<usize> = col_1neighbors.iter().map(|x| *x).collect();
    let r3: Vec<usize> = col_3neighbors.iter().map(|x| *x).collect();
    let cols_order = centered_ordering(cols_tm, r1, r3, markov_power);

    (rows_order, cols_order)
}





fn compute_cost(m: usize, split_threshold: f64, rows_cluster: &Vec<usize>, cols_cluster: &Vec<usize>, wadj: &WeightedBiAdjacency) -> f64{
    let mut cost = 0.0;
    
    for &row in rows_cluster {

        for col in cols_cluster {
            let mut found = false;
            for (col2, &w) in wadj.iter(row+m) {
                if col2 == col {
                    cost += 1. - w; // addition
                    found = true;
                    break;
                }
            }
            if !found {
                cost += 1.; // addition
            }
        }
    }

    for &row in rows_cluster {
        let mut deletion_cost = 0.;
        for (col, &w) in wadj.iter(row+m) {
            if cols_cluster.contains(col) == false {
                deletion_cost += w;
            }
        }
        if deletion_cost <= split_threshold {
            cost += deletion_cost // delete edges
        } else {
            cost += 1.; // split row
        }
    }

    for &col in cols_cluster {
        let mut deletion_cost = 0.;
        for (row, &w) in wadj.iter(col) {
            if rows_cluster.contains(row) == false {
                deletion_cost += w;
            }
        }
        if deletion_cost <= split_threshold {
            cost += deletion_cost // delete edges
        } else {
            cost += 1.; // split col
        }
    }

    cost
}

///
/// cost_coef in [0,1]
fn best(n: usize, m: usize, rows_order: Vec<(usize, f64)>, cols_order: Vec<(usize, f64)>, wadj: &WeightedBiAdjacency, cost_coef: f64, split_threshold: f64) 
-> (f64, Vec<usize>, Vec<usize>)  {
    let mut min_cost = f64::NAN;
    let mut min_cost_rows_cluster = Vec::new();
    let mut min_cost_cols_cluster = Vec::new();

    let c = cols_order.len();

    let mut rows_cluster = Vec::new();

    for i in 0..rows_order.len() {
        rows_cluster.push(rows_order[i].0);

        let mut cols_cluster = Vec::new();
        for j in 0..c {
            cols_cluster.push(cols_order[j].0);


            let cost = compute_cost(m, split_threshold, &rows_cluster, &cols_cluster, wadj);
            let s = i+2+j; // = rows_cluster.len() + cols_cluster.len()
            let cost = cost*f64::powf(s as f64, -cost_coef);
            
            // if cost > min_cost*4. {
            //     break;
            // }

            if min_cost.is_nan() || cost < min_cost {
                min_cost = cost;
                min_cost_rows_cluster = rows_cluster.clone();
                min_cost_cols_cluster = cols_cluster.clone();
            }
        }
    }
    
    (min_cost, min_cost_rows_cluster, min_cost_cols_cluster)
}









/// 
/// 
pub fn bicluster_two_sided( wadj: &mut WeightedBiAdjacency,
    cost_coef: f64,
    split_threshold: f64, 
    markov_power: usize,
    verbose: usize) -> (Biclust, AlgoStats) {

    let min_error = wadj.compute_min_error();


    let n = wadj.get_n();
    let m = wadj.get_m();
    let mut nb_operations = 0.;
    let mut nb_deletions = 0.;
    let mut nb_splits = 0.;
    let mut nb_additions = 0.;

    let mut a_clusters = vec![];
    let mut b_clusters = vec![];

    let mut assigned = vec![false; n+m];
    let mut nb_assigned = 0;

    let start_instant = Instant::now();

    // While there exists some B vertices to cluster
    while nb_assigned < n+m {
        
        
        
        if verbose == 0 {
            progress_bar(nb_assigned, n+m, start_instant);
        }

        // Check if there is a B vertex of degree 0
        for b in 0..m {
            if assigned[b+n] == false && wadj.col_degree(b) == 0 {
                assigned[b+n] = true;
                nb_assigned += 1;
            }
        }

        // Check if there is a row vertex of degree 0
        for a in 0..n {
            if assigned[a] == false && wadj.row_degree(a) == 0 {
                assigned[a] = true;
                nb_assigned += 1;
            }
        }

        if verbose >= 1 {
            println!("####################");
            println!("{nb_assigned} / {} vertices are assigned", n+m);
            let mut nb_edges = 0;
            for col in 0..m {
                for _ in  wadj.iter(col){
                    nb_edges += 1;
                }
            }
            println!("nb edges: {nb_edges}");
        }

        // Compute the transition matrix between B vertices
        let cols_tm = transition_matrix_b(wadj, n, m);
        let rows_tm = rows_transition_matrix(wadj, n, m);


        if verbose >= 3 {
            println!("Rows transition matrix:");
            print_matrix(&rows_tm);
            println!("Cols transition matrix:");
            print_matrix(&cols_tm);
        }

        let mut best_rows_cluster = Vec::new();
        let mut best_cols_cluster = Vec::new();
        let mut best_cost = f64::NAN;


        // V2
        for col in 0..m {
            if assigned[col+n] { continue; }
            let (rows_cluster,cols_cluster, cost) = best_cluster_col(wadj, n, m, col, &rows_tm, &cols_tm, markov_power, split_threshold, cost_coef, verbose);
            if best_cost.is_nan() || cost < best_cost {
                best_cost = cost;
                best_rows_cluster = rows_cluster.clone();
                best_cols_cluster = cols_cluster.clone();
            }
        }

        for row in 0..n {
            if assigned[row] { continue; }
            let (rows_cluster,cols_cluster, cost) = best_cluster_row(wadj, n, m, row, &rows_tm, &cols_tm, markov_power, split_threshold, cost_coef, verbose);
            if best_cost.is_nan() || cost < best_cost {
                best_cost = cost;
                best_rows_cluster = rows_cluster.clone();
                best_cols_cluster = cols_cluster.clone();
            }
        }

        // V1
        // Find the B_cluster with minimal cost
        for col in 0..m {
            break; // ###########################
            if assigned[col+n] {
                continue;
            }
            break;

            // Compute order of N^3[b]
            let mut d = 0.;
            for (_,w) in wadj.iter(col) {
                d += w;
            }
            let (rows_order, cols_order) = 
                if d == 0. {
                    (vec![], vec![(col,1.)])
                } else {
                    compute_bicluster_col(&wadj, m, col, &rows_tm, &cols_tm, markov_power, verbose)
                };
            if verbose >= 2 {
                println!("--- Step 1: compute rows_order of {col}");
                let mut d = 0.;
                for (_,w) in wadj.iter(col) {
                    d += w;
                }
                println!("degree: {}", d);
                println!("adjacent_rows: {rows_order:.2?}");
                println!("adjacent cols: {cols_order:.2?}");
            }
            
            // Compute best cost
            let (cost, rows_cluster, cols_cluster) = best(n, m, rows_order, cols_order, &wadj, cost_coef, split_threshold);
            if best_cost.is_nan() || cost < best_cost {
                best_cost = cost;
                best_rows_cluster = rows_cluster.clone();
                best_cols_cluster = cols_cluster.clone();
            }
            if verbose >= 2 {
                println!("cost: {cost:.3}\n rows_cluster: {rows_cluster:?}, cols_cluster: {cols_cluster:?}");
                println!("--- Step 2: best_cost: {best_cost:.3}, rows_cluster: {best_rows_cluster:?}, cols_cluster: {best_cols_cluster:?}");
            }

            if best_cost == 0. {
                break;
            }
        }

        // Find the Row centered biclusters with minimal cost
        for row in 0..n {
            break; // #######################
            if assigned[row] {
                continue;
            }

            // Compute order of N^3[b]
            let mut d = 0.;
            for (_,w) in wadj.iter(row+m) {
                d += w;
            }
            let (rows_order, cols_order) = 
                if d == 0. {
                    (vec![(row,1.)], vec![])
                } else {
                    compute_bicluster_row(&wadj, m, row, &rows_tm, &cols_tm, markov_power, verbose)
                };
            if verbose >= 2 {
                println!("--- Step 1: compute centered near element of Row {row}");
                println!("near rows: {rows_order:.2?}");
                println!("near cols: {cols_order:.2?}");
            }
            
            // Compute best cost
            let (cost, rows_cluster, cols_cluster) = best(n, m, rows_order, cols_order, &wadj, cost_coef, split_threshold);
            if best_cost.is_nan() || cost < best_cost {
                best_cost = cost;
                best_rows_cluster = rows_cluster.clone();
                best_cols_cluster = cols_cluster.clone();
            }
            if verbose >= 2 {
                println!("cost: {cost:.3}\n rows_cluster: {rows_cluster:?}, cols_cluster: {cols_cluster:?}");
                println!("--- Step 2: best_cost: {best_cost:.3}, rows_cluster: {best_rows_cluster:?}, cols_cluster: {best_cols_cluster:?}");
            }

            if best_cost == 0. {
                break;
            }
        }


        if best_cols_cluster.len() == 0 { // This case should not happen
            println!("break");
            break
        }
        
        if verbose >= 1 {
            println!("----------------");
            let s = best_rows_cluster.len() + best_cols_cluster.len();
            println!("Best bicluster: nrows={} ncols={} cost={:.3} rawcost={:.3}", best_rows_cluster.len(), best_cols_cluster.len(), best_cost, best_cost/(s as f64).powf(-cost_coef));
            best_rows_cluster.sort();
            best_cols_cluster.sort();
            println!("rows: {best_rows_cluster:?}");
            println!("cols: {best_cols_cluster:?}");
        }
        
        b_clusters.push(best_cols_cluster.clone());

        if verbose >= 2 {
            println!("--- Step 3: apply");
        }
        let ( del, add, spl )= apply_operations(n, m, &best_rows_cluster, &best_cols_cluster, wadj, split_threshold, verbose); 
        nb_deletions += del;
        nb_additions += add;
        nb_splits += spl;
        nb_operations += del + add + spl;
        if verbose >= 1 {
            println!("nb_del: {del:.1} nb_add: {add:.1} nb_splits: {spl:.0}");
        }

        a_clusters.push(best_rows_cluster);
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


    
    // compute_clusters(&a_clusters, &b_clusters, n, m)
    (Biclust::from_separate_biclusters(n,m, &a_clusters, &b_clusters), bicluster_stats)

}      










fn apply_operations(n: usize, m: usize, rows_cluster: &Vec<usize>, cols_cluster: &Vec<usize>, wadj: &mut WeightedBiAdjacency, split_threshold: f64, verbose: usize) -> (f64, f64, f64) {
    let mut nb_deletions = 0.;
    let mut nb_splits = 0.;
    let mut nb_additions = 0.;

    let mut edges_to_delete = Vec::new();
    for &row in rows_cluster.iter() {
        let mut indegree = 0;
        for (col, _) in wadj.iter(row+m) {
            if cols_cluster.contains(col) {
                indegree += 1;
                // println!("DEL {row} {col}");
                edges_to_delete.push((row, *col))
            } 
        }
        nb_additions += (cols_cluster.len() - indegree) as f64;
        // println!("add {}", cols_cluster.len() - indegree);
    }

    for &row in rows_cluster {
        let mut deletion_cost = 0.;
        for (col, &w) in wadj.iter(row+m) {
            if cols_cluster.contains(col) == false {
                deletion_cost += w;
            }
        }
        if deletion_cost <= split_threshold {
            for (col, &w) in wadj.iter(row+m) {
                if cols_cluster.contains(col) == false {
                    nb_deletions += w;
                    // println!("del {row} {col}");
                    edges_to_delete.push((row, *col));
                }
            }
        } else {
            // println!("spl row {row}");
            nb_splits += 1.;
        }
    }

    for &col in cols_cluster {
        let mut deletion_cost = 0.;
        for (row, &w) in wadj.iter(col) {
            if rows_cluster.contains(row) == false {
                deletion_cost += w;
            }
        }
        if deletion_cost <= split_threshold {
            for (row, &w) in wadj.iter(col) {
                if rows_cluster.contains(row) == false {
                    nb_deletions += w;
                    // println!("del {row} {col}");
                    edges_to_delete.push((*row, col))
                }
            }
        } else {
            // println!("spl col {col}");
            nb_splits += 1.;
        }
    }

    for (row, col) in edges_to_delete {
        wadj.delete_edge(row, col);
    }
    

    (nb_deletions, nb_additions, nb_splits)
}



    


pub fn analyze_ground_biclusters(wadj: &WeightedBiAdjacency){
    let n = wadj.get_n();
    let m = wadj.get_m();
    if let Some(biclusters) = wadj.get_ground_truth() {
        for bicluster in biclusters.biclusters(){
            println!("{bicluster:?}");
        }
        let mut sum = 0;
        for row in 0..n{
            let mut overlap = vec![];
            for (i, bicluster) in biclusters.biclusters().iter().enumerate(){
                if bicluster.contains(&row){
                    overlap.push(i);
                }
            }
            sum += overlap.len();
            if overlap.len() > 1 {
                println!("row {row}: {overlap:?}");
            }
        }

        for col in 0..wadj.get_m(){
            let mut overlap = vec![];
            for (i, bicluster) in biclusters.biclusters().iter().enumerate(){
                if bicluster.contains(&(col+n)){
                    overlap.push(i);
                }
            }
            sum += overlap.len();
            if overlap.len() > 1 {
                println!("col {col}: {overlap:?}");
            }
        }

        println!("{sum} {:.3}", (sum as f64)/((n+m) as f64));
    }
    
}


