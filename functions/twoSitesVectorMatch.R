##' twoSitesVectorMatch
##'
##' Simulates the coevolutionary dynamics on a mutualistic network
##' of species interactions at two sites connected by gene flow.
##'
##' This function has some alternate output from a more standard simulation output;
##' see return statement for details.
##' 
##' @param network.data: list with two elements:
##'     graph: bipartite network
##'     n.sp: dims of graph (flower/pollinator species, for example)
##'
##' @param g gene flow (symmetric between sites A and B) for all species
##' @param phi 'heritability' value
##' @param alpha parameter alpha, sensitivity of selection to trait matching
##' @param theta.A vector of environmental optimum values for site A
##' @param theta.B vector of environmental optimum values for site B
##' @param m.A vector of proportion of selection due to mutualism at site A
##' @param m.B vector of proportion of selection due to mutualism at site B
##' @param init.A initial values for site A
##' @param init.B initial values for site B
##' @param epsilon value to determine when equilibrium is reached
##' @param t.max maximum number of timesteps allowed
##'     
##' @details 
##'     All vectors need to have first the row species attributes
##'     and then the column species attributes
##'     (e.g. c(row_sp[1], ..., row_sp[nrow], col_sp[1], ..., col_sp[ncol]))
##'
##'     This function assumes that all species have the same trait 'heritability'
##'     and that gene flow has the same rate for all species.
##' 
##' @return list with three values:
##'     Q.eval: matrix with eigenvalues of Q for all iterations
##'     Q.evec: array of eigenvectors of Q in equilibrium
##'     match.corr: vector correlation between traits and
##'     optimal trait matching (all traits in both sites are equal)
##'     final.traits: trait values at end of iterations
##' 
twoSitesVectorMatch <- function(graph, g, phi = 1, alpha, theta.A, theta.B,
                                m.A, m.B, init.A, init.B, epsilon = 1e-6, t.max = 100000)
    {

        Norm <- function(x) sqrt(sum (x * x))
        Normalize <- function(x) x / Norm (x)
        
        n.sp <- sum(dim(graph))
        
        ## square adj matrix
        f <- assembleQ(graph)
        
        ## same gene flow and 'heritability'
        g <- rep(g, times = n.sp)
        phi <- rep(phi, times = n.sp)
        
        Z <- array(0, c(t.max, n.sp, 2))
        dimnames(Z) [[3]] <- c('A', 'B')

        ## vector of ones (theoretical matching)
        theo.match <- rep(1, times = n.sp * 2)
        theo.match <- Normalize(theo.match)
        
        ## zero matrix
        zeros <- diag(0, n.sp)

        I <- diag(1, n.sp * 2)
        
        ## matrices
        M <- rbind(cbind(diag(rep(m.A, n.sp)), zeros),
                   cbind(zeros, diag(rep(m.B, n.sp))))
        
        G <- rbind(cbind(diag(1 - g, n.sp), diag(g, n.sp)),
                   cbind(diag(g, n.sp), diag(1 - g, n.sp)))

        Phi <- diag(rep(phi, 2))
        
        Ginv <- solve(G)
        
        ## initial trait values 
        Z [1, , 'A'] <- init.A
        Z [1, , 'B'] <- init.B
        
        ## simulation runs for a maximum of t.max timesteps
        for (t in 1:(t.max - 1))
        {
            ## current z values
            z.A <- Z[t, , 'A'] 
            z.B <- Z[t, , 'B']
            
            ## matrix with all trait differences
            z.dif.A <- t(f*z.A) - f*z.A 
            z.dif.B <- t(f*z.B) - f*z.B
            
            ## calculating matrix q
            q.A <- f*(exp(-alpha * (z.dif.A^2))) 
            q.B <- f*(exp(-alpha * (z.dif.B^2)))
            
            ## normalizing the matrix
            q.n.A <- q.A / apply(q.A, 1, sum) 
            q.n.B <- q.B / apply(q.B, 1, sum)
            
            ## assemble Q matrix
            Q <- cbind(rbind(q.n.A, zeros), rbind(zeros, q.n.B))

            ## assemble T
            T <- solve(Ginv - I + Phi %*% (I - M %*% Q))

            var.match <- T %*% theo.match
            
            ## multiplying each row i of matrix q by m[i]
            q.m.A <- q.n.A * m.A 
            q.m.B <- q.n.B * m.B
            
            ## calculating selection differentials
            sel.dif.A <- q.m.A * z.dif.A 
            sel.dif.B <- q.m.B * z.dif.B
            
            ## response to selection related to mutualism
            r.mut.A <- phi * apply(sel.dif.A, 1, sum) 
            r.mut.B <- phi * apply(sel.dif.B, 1, sum)
            
            ## response to selection related to the environment
            r.env.A <- phi * (1 - m.A) * (theta.A - z.A) 
            r.env.B <- phi * (1 - m.B) * (theta.B - z.B)

            z.next.A <-
                (1 - g) * (z.A + r.mut.A + r.env.A) +
                      g * (z.B + r.mut.B + r.env.B) 

            z.next.B <-
                (1 - g) * (z.B + r.mut.B + r.env.B) +
                      g * (z.A + r.mut.A + r.env.A) 
            
            ## updating z values
            Z[t+1, , 'A'] <- z.next.A
            Z[t+1, , 'B'] <- z.next.B
                    
            traits <- c(Z[t+1, , 'A'], Z[t+1, , 'B'])
            traits <- Normalize(traits)
            trait.match.corr [t] <- traits %*% theo.match
            
            ## computing the mean difference between old and new z values
            dif.A <- mean(abs(Z[t+1, , 'A'] - Z[t, , 'A'])) 
            dif.B <- mean(abs(Z[t+1, , 'B'] - Z[t, , 'B']))

                ## condition for stopping iterations
            if ((dif.A < epsilon) & (dif.B < epsilon))
                break
        }

        return(list('T.eq' = T,
                    'T.dist' = Tdist,
                    'match.corr' = trait.match.corr,
                    'initial.traits' = Z[1, , ], 
                    'final.traits' = Z[t+1, , ]))        
    }
