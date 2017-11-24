buildTmat <- function(graph, g, phi = 0.5, alpha = 0.2,
                      m.A, m.B, theta.A, theta.B, value.A, value.B)
    {
        Norm <- function(x) sqrt(sum (x * x))
        Normalize <- function(x) x / Norm (x)

        n.sps <- dim(graph)
        
        n.sp <- sum(n.sps)
        
        ## square adj matrix
        f <- assembleQ(graph)
        
        ## same gene flow and 'heritability'
        g <- rep(g, times = n.sp)
        phi <- rep(phi, times = n.sp)
        
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

        ## current z values
        z.A <- value.A
        z.B <- value.B
            
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

        ## laplacian
        Lap <- Ginv - I + Phi %*% (I - M %*% Q)
        
        ## assemble T
        T <- solve(Lap) %*% Phi %*% (I - M)

        Lap.eig <- eigen(Lap)
        
        T.eig <- eigen(T)
        
        ## metrics
        eval2 <- T.eig $ values [2]
        evar <- var(Re(T.eig $ values))
        perc <- 1 / sum(Re(T.eig $ values))

        Lap.evals <- Lap.eig $ values [1:2]
        
        complex <- as.numeric(any(abs(Im(T.eig $ values)) > .Machine$double.eps))
        
        eval.count <- table(cut(Re(T.eig $ values), seq(0, 1, 0.01)))
        eval.freq <- eval.count / sum(eval.count)

        eval.freq <- eval.freq[eval.freq != 0]
        evalS <- - sum(eval.freq * log(eval.freq, 2))

        tm.A <- mean(convMutNet(sum(n.sps), n.sps [1], n.sps [2],
                                z.A, 'exponential', alpha))

        tm.B <- mean(convMutNet(sum(n.sps), n.sps [1], n.sps [2],
                                z.B, 'exponential', alpha))

        sd.lines <- mean(apply(T, 2, sd))

        T.norm <- apply(T, 1, Normalize)

        T.corr <- T.norm %*% t(T.norm)

        T.mcorr <- mean(T.corr[lower.tri(T.corr)])
        
        c('traitmatch.A' = tm.A, 'traitmatch.B' = tm.B,
          'eval2' = eval2, 'evar' = evar, 'perc' = perc,
          'Lap.eval1' = Lap.evals [1], 'Lap.eval2' = Lap.evals [2],
          'sd.lines' = sd.lines, 'mcorr.lines' = T.mcorr,
          'evalS' = evalS, 'complex' = complex)
    }
