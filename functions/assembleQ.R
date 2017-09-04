##' assembleQ
##'
##' Simple function that constructs a square matrix from a bipartite network
##'
##' @param graph bipartite network
##'
##' @return square adjacency matrix
##'
##' @details row species must come first
##' 
assembleQ <- function(graph)
    {
        n.sps <- dim(graph)
        rbind(cbind(diag(0, n.sps[1]), graph),
              cbind(t(graph), diag(0, n.sps[2])))
    }
