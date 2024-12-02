use std::{collections::HashMap, fs::File};
use std::io::Write;

use ndarray::Array2;
use rand::Rng;

use crate::common::print_matrix;

use super::biclust::Biclust;
use super::biclustering::Biclustering;


#[derive(Clone)]
pub struct WeightedBiAdjacency {
    wadj: Vec<HashMap<usize, f64>>, 
    n: usize,
    m: usize,
    labels_a: Vec<String>,
    labels_b: Vec<String>,
    nodes_a_map: HashMap<String, usize>,
    nodes_b_map: HashMap<String, usize>,
    ground_truth: Option<Biclust>
}


impl WeightedBiAdjacency {
    pub fn new(n: usize, m: usize) -> Self {
        let mut labels_a = vec![];
        let mut nodes_a_map = HashMap::new();
        for a in 0..n {
            labels_a.push(a.to_string());
            nodes_a_map.insert(a.to_string(), a);
        }
        let mut labels_b = vec![];
        let mut nodes_b_map = HashMap::new();
        for b in 0..m {
            labels_b.push(b.to_string());
            nodes_b_map.insert(b.to_string(), b);
        }

        WeightedBiAdjacency {
            wadj: vec![HashMap::new(); n+m],
            n,
            m,
            labels_a,
            labels_b,
            nodes_a_map,
            nodes_b_map,
            ground_truth: None
        }
    }


    pub fn write_to_file(&self, filename: &str)  {
        let mut file = File::create(filename).unwrap();

        writeln!(file, "# n={} m={}", self.n, self.m).unwrap();

        for a in 0..self.n {
            for (&b,_) in self.wadj[a+self.m].iter() {
                writeln!(file, "{a} {b}").unwrap();
            }
        }
       
    }

    pub fn get_n(&self) -> usize {
        self.n
    }

    pub fn get_m(&self) -> usize {
        self.m
    }

    pub fn get_labels(&self) -> (Vec<String>, Vec<String>, HashMap<String, usize>, HashMap<String, usize>) {
        (self.labels_a.clone(), self.labels_b.clone(), self.nodes_a_map.clone(), self.nodes_b_map.clone())
    }

    pub fn iter(&self, idx: usize) -> impl Iterator<Item = (&usize, &f64)> {
        self.wadj[idx].iter()
    }

    pub fn col_degree(&self, b: usize) -> usize {
        self.wadj[b].len()
    }

    pub fn has_edgee(&self, a: usize, b: usize) -> bool {
        self.wadj[b].contains_key(&a)
    }

    pub fn delete_edge(&mut self, a: usize, b: usize) {
        self.wadj[a+self.m].remove(&b);
        self.wadj[b].remove(&a);
    }
    
    pub fn add_edge(&mut self, a: usize, b: usize, weight: f64) {
        self.wadj[a+self.m].insert(b, weight);
        self.wadj[b].insert(a, weight);
    }

    fn toggle(&mut self, i: usize,j: usize) {
        if self.wadj[j].contains_key(&i) {
            self.delete_edge(i, j);
        } else {
            self.add_edge(i, j, 1.);
        }
    }

    pub fn print(&self) {
        for row in 0..self.n {
            println!("row {row} {:?}", self.wadj[self.m + row]);
        }
        for col in 0..self.m {
            println!("col {col} {:?}", self.wadj[col]);
        }
    }

    

    /// p is the noise
    pub fn rand(n: usize, m: usize, noise: f64, row_overlap: f64, row_separation: f64  ) -> Self{
        let mut wadj = WeightedBiAdjacency::new(n,m);
        let mut rng = rand::thread_rng();

        // For each A vertex, generate the number of biclusters which will contains it
        // At least 1. Using a exponential law
        let mut nb_biclusters = vec![0;n];
        let mut max_nb_bic = 1;
        let mut sum = 0;
        for a in 0..n {
            let c = generate_exponential_random(row_overlap);
            nb_biclusters[a] = c;
            sum += c;
            if c > max_nb_bic {
                max_nb_bic = c;
            }
        }
        let k = max_nb_bic+generate_binomial(row_separation, sum-max_nb_bic);
        println!("max_membership: {max_nb_bic}, nb_biclusters: {k}");

        let mut biclusters_indices = vec![0;k];
        for i in 0..k {
            biclusters_indices[i] = i;
        }
        let mut biclusters = vec![vec![]; k];
        // For each a vertex, choose nb_biclusters[a] among the k biclusters
        for a in 0..n {
            shuffle(&mut biclusters_indices);
            for j in 0..nb_biclusters[a] {
                biclusters[biclusters_indices[j]].push(a);
            }
        }

        // For each b vertex, choose a random bicluster among the k biclusters
        for b in 0..m {
            shuffle(&mut biclusters_indices);
            for &a in biclusters[biclusters_indices[0]].iter() {
                if a < n {
                    wadj.add_edge(a, b, 1.);
                }
            }
            biclusters[biclusters_indices[0]].push(b+n);
        }


        
        let mut gt = Biclust::from_biclusters(n, m, &biclusters);
        gt.print();
        gt.reduce_isolated();
        wadj.ground_truth = Some(gt);


        // Noise
        for i in 0..n  {
            for j in 0..m {
                if rng.gen_range(0.0..1.0) < noise {
                    wadj.toggle(i, j);
                }
            }
        }

        wadj
    }


    

    fn get_matrix(&self) -> Array2<f64> {
        let mut matrix = Array2::zeros((self.n,self.m));
        for i in 0..self.n {
            for j in 0..self.m {
                if let Some(w) = self.wadj[j].get(&i) {
                    matrix[[i,j]] = *w;
                }
            }
        }
        matrix
    }

    pub fn print_matrix(&self) {
        print_matrix(&self.get_matrix());
    }

    pub fn get_ground_truth(&self) -> Option<Biclust> {
        self.ground_truth.clone()
    }


    pub fn compute_ground_truth_noise(&self) -> f64 {
        if let Some(ground_truth) = &self.ground_truth {
            self.compute_noise(&ground_truth)
        } else {
            0.
        }
    }

    /// Density is a number in [0,1]
    /// Density of 1 if the graph is a biclique
    /// Density of 0 if the graph is empty
    pub fn density(&self) -> f64 {
        let mut s = 0;
        for col in 0..self.m {
            s += self.wadj[col].len()
        }
        s as f64 / (self.n * self.m) as f64
    }

    /// Noise is the ratio of false edges over the number of possible edges (n*m)
    pub fn compute_noise(&self, biclusters: &Biclust) -> f64 {
        let mut nb_erros = 0;

        for row in 0..self.n {
            for col in 0..self.m {
                let mut has_edge = false;
                for bicluster in biclusters.biclusters() {
                    if bicluster.contains(&row) && bicluster.contains(&(self.n + col)) {
                        has_edge = true;
                        break
                    }
                }
                if has_edge != self.has_edgee(row, col) {
                    nb_erros += 1;
                }
            }
        }

        nb_erros as f64 / (self.n * self.m) as f64
    }


    /// If weight is 0 or 1 then min error is 0
    /// If weight is 0.5 then min error is 0.5 for this edge. You need at least 0.5 operation to solve it
    pub fn compute_min_error(&self) -> f64 {
        let mut r = 0.;
        for row in 0..self.n {
            for col in 0..self.m {
                if let Some(&w) = self.wadj[col].get(&row) {
                    if w > 0.5 {
                        r += 1.- w;
                    } else {
                        r += w
                    }
                }
           } 
        }
        r / (self.n * self.m) as f64
    }

}




fn generate_exponential_random(p: f64) -> usize {
    let mut r = 1;

    let mut rng = rand::thread_rng();
    loop {
        let roll = rng.gen::<f64>();
        if p*roll < 1.  {
            return r;
        } else {
            r += 1;
        }
    }
}


fn generate_binomial(p: f64, n: usize) -> usize {
    let mut r = 0;
    let mut rng = rand::thread_rng();
    for _ in 0..n {
        let roll = rng.gen::<f64>();
        if roll < p  {
            r += 1;
        }
    }
    r
}


fn shuffle<T>(vec: &mut Vec<T>) {
    let mut rng = rand::thread_rng();
    for i in (1..vec.len()).rev() {
        let j = rng.gen_range(0..=i);
        vec.swap(i, j);
    }
}
