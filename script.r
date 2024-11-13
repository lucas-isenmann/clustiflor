
# Create the biadjcency matrix from a file containing the list of the edges of a bipartite graph
create_graph_biadj_matrix <- function(filename) {
    lines <- readLines(filename)

    n = 0
    m = 0
    nodes_a = numeric()
    nodes_b = numeric()
    edges = numeric()
    row_names = c()
    col_names = c()

    for (line in lines) {
        if (grepl("^#", line)) {
            next
        }
        t = strsplit(line, "\\s+")[[1]]   
        a = t[1]
        b = t[2]

        if ( !(a %in% names(nodes_a))){
            n = n + 1
            nodes_a[a] = n
            row_names = c(row_names, a)
        }
        if ( !(b %in% names(nodes_b))){
            m = m + 1
            nodes_b[b] = m
            col_names = c(col_names, b)
        }
        edges = c( edges, nodes_a[a], nodes_b[b])
    }
        
    graph_matrix <- matrix(0, nrow=n, ncol=m, dimnames=list(row_names, col_names))
    for( i in 1:(length(edges)/2) ) {
        graph_matrix[edges[2*i-1], edges[2*i]] <- 1
    }
    print(graph_matrix)
  return(graph_matrix)
}

graph_matrix = create_graph_biadj_matrix("test.adj")

library(biclust)
res = biclust(x=graph_matrix, method=BCBimax())

writeBiclusterResults("results.txt", res,"Bimax biclusters", rownames(graph_matrix), colnames(graph_matrix))




