pub trait Biclustering {
    /// Get the number of rows in the biclustering
    fn rows(&self) -> usize;

    /// Get the number of columns in the biclustering
    fn cols(&self) -> usize;

    /// Get the biclusters represented as a Vec<Vec<usize>>
    fn biclusters(&self) -> Vec<Vec<usize>>;

    /// Get the number of unclustered rows
    fn unclustered_rows(&self) -> Vec<usize>;

    // Get the number of unclustered columns
    fn unclustered_cols(&self) -> Vec<usize>;

    fn add_bicluster(&mut self, bicluster: Vec<usize>);
}