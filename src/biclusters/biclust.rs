use std::{collections::HashSet, fs::File};
use std::io::Write;
use super::biclustering::Biclustering;

#[derive(Clone)]
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


    /// a_clusters values should be in [0,n-1]
    /// b_clusters values should be in [0,m-1]
    /// 
    pub fn from_separate_biclusters(n: usize, m: usize, a_clusters: &Vec<Vec<usize>>, b_clusters: &Vec<Vec<usize>>)  -> Self {
        let mut b = Biclust::new(n,m);

        for i in 0..a_clusters.len() {
            let mut c = a_clusters[i].clone();
            for x in b_clusters[i].iter() {
                c.push(n+*x);
            }
            b.add_bicluster(c);
        }
        b
    }

    pub fn from_biclusters(n: usize, m: usize, biclusters: &Vec<Vec<usize>>)  -> Self {
        let mut b = Biclust::new(n,m);
        for bicluster in biclusters.iter() {
            b.add_bicluster(bicluster.clone());
        }
        
        b
    }

    pub fn clear(&mut self) {
        self.biclusters.clear();
        for row in 0..self.n {
            self.rows_memberships[row].clear();
        }
        for col in 0..self.m {
            self.cols_memberships[col].clear();
        }
    }
    

    /// Remove biclusters consisting only of rows or of columns
    /// Add a 
    pub fn reduce_isolated(&mut self){
        let mut biclusters_non_trivial = vec![];
        let mut isolated_rows = vec![];
        let mut isolated_cols = vec![];

        for bicluster in self.biclusters.iter() {
            // Check if bicluster is isolated
            let mut has_row = false;
            let mut has_col = false;
            for &x in bicluster {
                if x < self.n {
                    has_row = true;
                } else {
                    has_col = true;
                }
            }
            if has_row && has_col {
                biclusters_non_trivial.push( bicluster.clone());
            } else {
                for &x in bicluster {
                    if x < self.n {
                        isolated_cols.push(x);
                    } else {
                        isolated_rows.push(x);
                    }
                }
            }
        }
        self.clear();
        for bicluster in biclusters_non_trivial.iter() {
            self.add_bicluster(bicluster.clone());
        }
        if isolated_cols.len() > 0 {
            self.add_bicluster(isolated_cols);
        } 
        if isolated_rows.len() > 0 {
            self.add_bicluster(isolated_rows);
        }

    }

    pub fn print(&self) {
        print!("[");
        for bicluster in self.biclusters.iter() {
            println!("{:?}", bicluster);
        }
        print!("]");
    }

    /// Matching score returning a float in [0,1]
    /// Geometric mean of the average of the maximum of the Jaccard index between a rows cluster of A and the rows cluster of B and the same for the cols
    /// 
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

    pub fn f_score(&self, other: &Biclust) -> f64 {
        let recall = proportion(&self.biclusters, &other.biclusters);
        let precision = proportion(&other.biclusters, &self.biclusters);
        2. * recall * precision / (recall + precision)
    }


    /// Return true iff there exist a bicluster containing i and j
    pub fn are_together(&self, i: usize, j: usize) -> bool{
        for bicluster in &self.biclusters {
            if bicluster.contains(&i) && bicluster.contains(&j){
                return true
            }
        }
        false
    }

    /// Return 1 if self and other agree on every (row,col) edges
    /// In general return the proportion of pairs in the rows x cols bipartition such that both biclusters agree on it.
    /// By "agree" we mean that either there exists in both a bicluster which contains row and col, either there exists no bicluster containing row and col in both
    pub fn accuracy(&self, other: &Biclust) -> f64 {
        let mut r = 0;
        for row in 0..self.n {
            for col in 0..self.m {
                if self.are_together(row, self.n + col) == other.are_together(row, self.n + col) {
                    r += 1
                }
            }
        }
        // for row1 in 0..self.n {
        //     for row2 in 0..row1 {
        //         if self.are_together(row1, row2) == other.are_together(row1, row2) {
        //             r += 1
        //         }
        //     }
        // }

        // r as f64 / (self.n * self.m + self.n*(self.n-1)/2 ) as f64
        r as f64 / (self.n * self.m) as f64
    }


    pub fn get_rows_overlapping(&self) -> f64 {
        let mut result = 0.;
        for cluster in self.biclusters.iter() {
            for &x in cluster {
                if x < self.n {
                    result += 1.;
                }
            }
        }
        result / (self.n as f64)
    }
    

    pub fn print_stats(&self,
        size_sensivity: f64, 
        split_threshold: f64, 
        markov_power: usize,
        labels_a: &Vec<String>, 
        labels_b: &Vec<String>, 
        file_path: Option<&str>) {
        
        

        let file_name = match file_path {
            Some(path) => path.to_string(),
            None => "biclusters.txt".to_string(),
        };

        let mut file = File::create(&file_name).expect("Failed to open file");
        
        writeln!(file, "# Hyperparameters").unwrap();
        writeln!(file, "- size sensivity: {size_sensivity}").unwrap();
        writeln!(file, "- split threshold: {split_threshold}").unwrap();
        writeln!(file, "- markov power: {markov_power}").unwrap();
        
        writeln!(file, "\n# Results").unwrap();

        
    // - Error: {eta:.6}
    // - Nb isolated A vertices: {}
    // - Nb isolated B vertices: {}
    // - Nb operations: {nb_operations}
    // - Nb splits: {nb_splits}
    // - Nb deletions: {nb_deletions}
    // - Nb additions: {nb_additions}

        writeln!(file, "- Number of biclusters: {}", self.biclusters.len()).unwrap();
        writeln!(file, "- A Overlapping: {:.3}", self.get_rows_overlapping()).unwrap();

        writeln!(file, "").unwrap();
        writeln!(file, "# Clusters\n").unwrap();

        for bicluster in self.biclusters.iter() {
            for &x in bicluster {
                if x < self.n {
                    write!(file, "{} ", labels_a[x]).unwrap();
                } else {
                    write!(file, "{} ", labels_b[x-self.n]).unwrap();
                }
            }
            writeln!(file, "").unwrap();
        }
        println!("Biclusters printed to {}", file_name);

    }

}

fn has_common(l1: &[usize], l2: &[usize]) -> bool {
    l1.iter().any(|&x| l2.contains(&x))
}

fn proportion(a1: &Vec<Vec<usize>>, a2: &Vec<Vec<usize>>) -> f64 {
    let mut nb_same_a1 = 0;
    let mut nb_same_a2 = 0;

    for i in 0..a1.len() {
        for j in 0..i {
            if has_common(&a1[i], &a1[j]) {
                nb_same_a1 += 1;
                if has_common(&a2[i], &a2[j]) {
                    nb_same_a2 += 1;
                }
            }
        }
    }

    if nb_same_a1 == 0 {
        0.0
    } else {
        nb_same_a2 as f64 / nb_same_a1 as f64
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
    let union_size = size_union(a, b);
    let inter_size = size_intersection(a,b);
    if union_size == 0 {
        1.
    } else {
        inter_size as f64 / union_size as f64
    }
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
    fn matching_score_tests() {
        // let mut biclust_a = Biclust::new(10, 10);
        // biclust_a.add_bicluster(vec![0,10]);

        // assert_eq!(biclust_a.matching_score(&biclust_a), 1.0);



        // // Same biclusters should give 1
        // let b1 = Biclust::from_biclusters(2, 2, 
        //     &vec![
        //         vec![0,2],
        //         vec![1,3]]);
        // let b2 = Biclust::from_biclusters(2, 2, 
        //     &vec![
        //         vec![0,2],
        //         vec![1,3]]);
        // assert_eq!(b1.matching_score(&b2), 1.0);  

        // // With one edge in more
        // let b1 = Biclust::from_biclusters(2, 2, 
        //     &vec![
        //         vec![0,1,2],
        //         vec![1,3]]);
        // let b2 = Biclust::from_biclusters(2, 2, 
        //     &vec![
        //         vec![0,2],
        //         vec![1,3]]);
        // // assert_eq!(b1.matching_score(&b2), 0.75);  

        // 
        let b1 = Biclust::from_biclusters(2, 2, 
            &vec![
                vec![0,2,3],
                vec![1,3]]);
        let b2 = Biclust::from_biclusters(2, 2, 
            &vec![
                vec![0,2],
                vec![0,1,3]]);
        assert_eq!(b1.matching_score(&b2), 0.75);  
        assert_eq!(b2.matching_score(&b1), 0.75);  
    }


    #[test]
    fn accuracy_test() {
        let b1 = Biclust::from_biclusters(2, 2, 
            &vec![
                vec![0,2],
                vec![1,3]]);
        let b2 = Biclust::from_biclusters(2, 2, 
            &vec![
                vec![0,2],
                vec![1,3]]);
        assert_eq!(b1.accuracy(&b2), 1.0);  

        let b1 = Biclust::from_biclusters(2, 2, 
            &vec![
                vec![0,1,2],
                vec![1,3]]);
        let b2 = Biclust::from_biclusters(2, 2, 
            &vec![
                vec![0,2],
                vec![1,3]]);
        assert_eq!(b1.accuracy(&b2), 0.75);  


        let b1 = Biclust::from_biclusters(2, 2, 
            &vec![
                vec![0,2,3],
                vec![1,3]]);
        let b2 = Biclust::from_biclusters(2, 2, 
            &vec![
                vec![0,2],
                vec![0,1,3]]);
        assert_eq!(b1.accuracy(&b2), 1.0);  

    }
}
