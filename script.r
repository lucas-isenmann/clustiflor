# library(biclust)

# edges <- cbind(c(1,2,2,3,3,4), c(2,3,3,4,4,4))
# graph_matrix <- matrix(0, nrow=4, ncol=4)
# for(i in 1:nrow(edges)) {
#   graph_matrix[edges[i,1], edges[i,2]] <- 1
# }
# graph_matrix

# res <- biclust(x=graph_matrix, method=BCBimax(), minr=2, minc=2, number=10)

# drawHeatmap(x=graph_matrix, bicResult=res, number=1)


create_graph_matrix <- function(filename) {
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

graph_matrix = create_graph_matrix("peter.adj")

library(biclust)
res = biclust(x=graph_matrix, method=BCBimax())

writeBiclusterResults("results.txt", res,"lol", rownames(graph_matrix), colnames(graph_matrix))




