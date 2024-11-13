use std::collections::HashMap;

use ndarray::Array2;
use rand::Rng;

use crate::common::print_matrix;



pub struct WeightedBiAdjacency {
    wadj: Vec<HashMap<usize, f64>>, 
    n: usize,
    m: usize
}


impl WeightedBiAdjacency {
    pub fn new(n: usize, m: usize) -> Self {
        WeightedBiAdjacency {
            wadj: vec![HashMap::new(); n+m],
            n,
            m,
        }
    }


    fn delete_edge(&mut self, a: usize, b: usize) {
        self.wadj[a+self.m].remove(&b);
        self.wadj[b].remove(&a);
    }
    
    fn add_edge(&mut self, a: usize, b: usize, weight: f64) {
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

    pub fn rand(n: usize, m: usize, p: f64) -> Self{
        let mut wadj = WeightedBiAdjacency::new(n,m);
        let mut rng = rand::thread_rng();

        for i in 0..n  {
            for j in 0..m {
                if rng.gen_range(0.0..1.0) < p {
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

    pub fn print(&self) {
        print_matrix(&self.get_matrix());
    }

}