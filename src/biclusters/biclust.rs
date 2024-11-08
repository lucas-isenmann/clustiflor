use std::collections::HashSet;

use super::biclustering::Biclustering;

pub struct Biclust {
    n: usize,
    m: usize,
    rows_memberships: Vec<Vec<usize>>,
    cols_memberships: Vec<Vec<usize>>,
    biclusters: Vec<Vec<usize>>,
}

impl Biclust {
    pub fn new(n: usize, m: usize) -> Self {
        Biclust {
            n,
            m,
            rows_memberships: vec![vec![]; n],
            cols_memberships: vec![vec![]; m],
            biclusters: vec![],
        }
    }

    pub fn print(&self) {
        for bicluster in self.biclusters.iter() {
            println!("{:?}", bicluster);
        }
    }

    /// Matching score
    /// Defined in https://appliednetsci.springeropen.com/articles/10.1007/s41109-019-0180-x
    /// Preli et al 2006
    /// Eren et al. 2012
    pub fn matching_score(&self, other: &Biclust) -> f64{

        // For rows
        let mut sum = 0.;
        for bicluster in self.biclusters.iter(){
            let mut rows = vec![];
            for &x in bicluster {
                if x < self.n {
                    rows.push(x);
                }
            }
            let mut r = 0.;
            for bicluster2 in other.biclusters.iter() {
                let mut rows2 = vec![];
                for &x in bicluster2 {
                    if x < other.n {
                        rows2.push(x);
                    }
                }
                let y = jaccard_index(&rows, &rows2);
                if y > r {
                    r = y;
                }
            }
            sum += r;
        }
        let nb_biclusters = self.biclusters.len();
        let srows = sum/ (nb_biclusters as f64);

        // For columns
        let mut sum = 0.;
        for bicluster in self.biclusters.iter(){
            let mut cols = vec![];
            for &x in bicluster {
                if x >= self.n {
                    cols.push(x);
                }
            }
            let mut r = 0.;
            for bicluster2 in other.biclusters.iter() {
                let mut cols2 = vec![];
                for &x in bicluster2 {
                    if x >= other.n {
                        cols2.push(x);
                    }
                }
                let y = jaccard_index(&cols, &cols2);
                if y > r {
                    r = y;
                }
            }
            sum += r;
        }
        let nb_biclusters = self.biclusters.len();
        let scols = sum/ (nb_biclusters as f64);

        (scols*srows).sqrt()
    }
}

fn size_intersection(a: &Vec<usize>, b: &Vec<usize>) -> usize {
    let mut r = 0;
    for x in a.iter() {
        if b.contains(x) {
            r += 1;
        }
    }
    r
}

fn size_union(a: &Vec<usize>, b: &Vec<usize>) -> usize {
    let mut u = HashSet::new();
    for &x in a.iter() {
        u.insert(x);
    }
    for &x in b.iter() {
        u.insert(x);
    }
    u.len()
}

fn jaccard_index(a: &Vec<usize>, b: &Vec<usize>) -> f64 {
    (size_intersection(a, b) as f64) / (size_union(a, b) as f64)
}





impl Biclustering for Biclust {
    fn rows(&self) -> usize {
        self.n
    }

    fn cols(&self) -> usize {
        self.m
    }

    fn biclusters(&self) -> Vec<Vec<usize>> {
        self.biclusters.clone()
    }

    fn unclustered_rows(&self) -> Vec<usize> {
        let mut r = vec![];
        for (i,x) in self.rows_memberships.iter().enumerate() {
            if x.len() == 0 {
                r.push(i);
            }
        }
        r
    }

    fn unclustered_cols(&self) -> Vec<usize> {
        let mut r = vec![];
        for (i,x) in self.cols_memberships.iter().enumerate() {
            if x.len() == 0 {
                r.push(i);
            }
        }
        r
    }

    fn add_bicluster(&mut self, bicluster: Vec<usize>) {
        let bicluster_id = self.biclusters.len();
        for &x in bicluster.iter() {
            if x < self.n {
                self.rows_memberships[x].push(bicluster_id);
            }
            else {
                self.cols_memberships[x-self.n].push(bicluster_id);
            }
        }
        self.biclusters.push(bicluster.clone());
    }

}








#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matching_score_self_equals_1() {
        let mut biclust_a = Biclust::new(10, 10);
        biclust_a.add_bicluster(vec![0,10]);

        assert_eq!(biclust_a.matching_score(&biclust_a), 1.0);
    }
}
