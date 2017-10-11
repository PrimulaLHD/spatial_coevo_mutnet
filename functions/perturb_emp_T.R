perturbEmpT <- function(graph, mA, mB, pert_g_seq,
                        theta_A_min = 0, theta_A_max = 10,
                        theta_B_min = 10, theta_B_max = 20,
                        n_theta = 30, n_rep_pert = 30)
    {
        mat <- graph
        n_row = nrow(mat)
        n_col = ncol(mat)
        n_sp = n_row + n_col

        ## Q and F matrices
        Q = rbind(cbind(matrix(0, n_row, n_row), mat),
                  cbind(t(mat), matrix(0, n_col, n_col)))
        
        F_matrix = Q

        F_matrix[F_matrix != 0] = 1

        Q_norm = Q/apply(Q, 1, sum)

        Q_2sites = rbind(cbind(Q_norm, matrix(0, n_sp, n_sp)),
                         cbind(matrix(0, n_sp, n_sp), Q_norm))

        ## identity matrix
        I = diag(rep(1, 2*n_sp))

        m_A = rep(mA, n_sp)
        m_B = rep(mB, n_sp)
        M = diag(c(m_A, m_B))

        out <-
            aaply(1:n_rep_pert, 1, function(i)
            {
                
                ## matrix G
                
                ## sampling species perturbation order
                sp_order = sample(1:n_sp, replace = FALSE)
                
                out.internal <- aaply(1:length(pert_g_seq), 1, function(j)
                {

                    ## perturbing matrix G
                    sp_index = round(pert_g_seq[j]*n_sp)
                    
                    g_mult = rep(g, n_sp)
                    
                    if (sp_index != 0) {
                        
                        pert_sp = sp_order[1:sp_index]
                        
                        g_mult[pert_sp] <- 0
                        
                    }
                    
                    G = rbind(cbind(diag(1 - g_mult), diag(g_mult)),
                              cbind(diag(g_mult), diag(1 - g_mult)))
                    
                    ## laplacian
                    Lap <- solve(G) -  M %*% Q_2sites
                    
                    eigen.Lap <- eigen(Lap) $ values
                    
                    eigen1.Lap <- eigen.Lap [1]
                    eigen2.Lap <- eigen.Lap [2]
                    
                    ## build T matrix
                    T_matrix = solve(Lap) %*% (I - M)
                    
                    ## T matrix eigenvalues
                    eigens = Re(eigen(T_matrix)$values)
                    
                    ## mean of eigenvalues
                    mean_eigens = mean(eigens)
                    
                    ## sd of eigenvalues
                    sd_eigens = sd(eigens)
                    
                    ## 1st and 2nd eigenvalues
                    eigen1 = eigens[1]
                    eigen2 = eigens[2]
                    
                    ## determinant
                    det_T = prod(eigens)

                    ## mean column sd
                    mean_col_sd = mean(apply(T_matrix, 2, sd))

                    prop_eigen_1 <- 1 / sum(eigens)
                    
                    ## theta vector to calculate trait matching
                    theta <-
                        rbind(matrix(runif(n_sp * n_theta,
                                           min = theta_A_min,
                                           max = theta_A_max), nrow = n_sp),
                              matrix(runif(n_sp * n_theta,
                                           min = theta_B_min,
                                           max = theta_B_max), nrow = n_sp))
                    
                    
                    ## equilibrium trait values
                    Zmat <- T_matrix %*% theta

                    match <- aaply(Zmat, 2, function(z)
                    {
                        match_A <-
                            MatchingMutNet(n_sp = n_sp,
                                           n_row = n_row, n_col = n_col,
                                           f = F_matrix, 
                                           z = z[1:n_sp],
                                           method = "exponential", alpha = 0.2) [[1]]

                        match_B <-
                            MatchingMutNet(n_sp = n_sp,
                                           n_row = n_row, n_col = n_col,
                                           f = F_matrix, 
                                           z = z[(n_sp+1):(2*n_sp)],
                                           method = "exponential", alpha = 0.2) [[1]]

                        c(match_A, match_B)
                    })

                    
                    return(c('eigen1.Lap' = eigen1.Lap, 'eigen2.Lap' = eigen2.Lap,
                             'mean_eigen' = mean_eigens, 'sd__eigen' = sd_eigens,
                             'eigen1' = eigen1, 'eigen2' = eigen2, 'det_T' = det_T,
                             'mean_col_sd' = mean_col_sd, 'prop_eigen' = prop_eigen_1, 
                             'matchA' = match [, 1], 'matchB' = match [, 2]))
                })
                
                return(out.internal)
            })
        
        return(out)

    }
