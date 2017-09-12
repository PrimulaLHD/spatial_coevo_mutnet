##' empiricalSpectralT
##'
##' Quantifies certain aspects of the spectral decomposition of T matrix at equilibrium,
##' assuming that empirical weighted bipartite graphs can be used as proxies for Q at
##' equilibrium.
##'
##' @param graph
##' @param g gene flow (symmetric between sites A and B) for all species
##' @param phi 'heritability' value
##' @param m.A vector of proportion of selection due to mutualism at site A
##' @param m.B vector of proportion of selection due to mutualism at site B
##'
##' @details
##'
##' @return
##'
empiricalSpectralT <- function(graph, g, phi, m.A, m.B)
    {
        f <- assembleQ(graph)

        ## Normalize rows
        Q <- f / rowSums(f)

        ## Build other matrices
        zeros <- diag(0, n.sp)
        I <- diag(1, n.sp * 2)
        
        ##
        M <- rbind(cbind(diag(rep(m.A, n.sp)), zeros),
                   cbind(zeros, diag(rep(m.B, n.sp))))
        
        G <- rbind(cbind(diag(1 - g, n.sp), diag(g, n.sp)),
                   cbind(diag(g, n.sp), diag(1 - g, n.sp)))

        Phi <- diag(rep(phi, 2))

        Ginv <- solve(G)

        ## assemble T
        T <- solve(Ginv - I + Phi %*% (I - M %*% Q)) %*% Phi %*% (I - M)

        T.eig <- eigen(T)
    }
    
