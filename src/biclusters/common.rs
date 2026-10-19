use ndarray::Array2;

use super::weighted_biadj::WeightedBiAdjacency;

pub struct AlgoStats {
    pub adjusted_error: f64,
    pub nb_operations: f64,
    pub nb_splits: f64,
    pub nb_deletions: f64,
    pub nb_additions: f64,
}



/// Return the transition matrix between rows
pub fn rows_transition_matrix(wadj: &WeightedBiAdjacency, n: usize, m: usize) -> Array2<f64> {
    let mut tm = Array2::zeros((n, n));

    let d = compute_degrees(wadj, n, m);

    // Compute tm[i][j] = probability i -> j
    // tm[i][j] = (sum wxi wxj) / (d[i] d[j])
    for i in 0..n {
        for j in 0..n {
            for (x, wxi) in wadj.row_adj[i].iter() {
                for (y, w) in wadj.row_adj[j].iter() {
                    if x == y {
                        tm[[i, j]] += (wxi * w) / (d[i + m] * d[*x]);
                        break;
                    }
                }
            }
        }
    }
    tm
}

/// Return the transition matrix between B vertices
///
pub fn transition_matrix_b(wadj: &WeightedBiAdjacency, n: usize, m: usize) -> Array2<f64> {
    let mut tm = Array2::zeros((m, m));

    let d = compute_degrees(wadj, n, m);

    // Compute tm[i][j] = probability i -> j
    // tm[i][j] = (sum wxi wxj) / (d[i] d[j])
    for i in 0..m {
        for j in 0..m {
            for (x, wxi) in wadj.col_adj[i].iter() {
                for (y, w) in wadj.col_adj[j].iter() {
                    if x == y {
                        tm[[i, j]] += (wxi * w) / (d[i] * d[x + m]);
                        break;
                    }
                }
            }
        }
    }
    tm
}

fn compute_degrees(wadj: &WeightedBiAdjacency, n: usize, m: usize) -> Vec<f64> {
    let mut degrees = vec![0.; n + m];
    for b in 0..m {
        degrees[b] = wadj.col_wdeg(b);
        // for (_, w) in wadj.iter(b) {
        //     degrees[b] += w;
        // }
    }
    for a in 0..n {
        degrees[a+m] = wadj.row_wdeg(a);
        // for (_, w) in wadj.iter(a + m) {
        //     degrees[a + m] += w;
        // }
    }
    degrees
}
