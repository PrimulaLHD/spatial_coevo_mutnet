AnalyticalJacobian <- function(z_A, z_B, alpha, m_A, m_B, h, f, g)
{
    z_dif_A = t(f * z_A) - f * z_A 
    z_dif_B = t(f * z_B) - f * z_B
    q_A = f * exp(-alpha * (z_dif_A^2))
    q_B = f * exp(-alpha * (z_dif_B^2))
    q_n_A = q_A / apply(q_A, 1, sum)
    q_n_B = q_B / apply(q_B, 1, sum)

    delQ.part.A <- 2 * alpha * z_A * t(z_dif_A)
    delQ.A <- (1 - delQ.part.A) * t(q_A) + delQ.part.A * t(q_n_A^2)

    delQ.part.B <- 2 * alpha * z_B * t(z_dif_B)
    delQ.B <- (1 - delQ.part.B) * t(q_B) + delQ.part.B * t(q_n_B^2)
    
    M <- P <- Q <- list()
    for(i in 1:length(z_A))
    {
        M [[i]] <- matrix(c(1 - g[i], g[i],
                            g[i], 1 - g[i]), byrow = T, ncol = 2)
        
        P [[i]] <- diag(2) - matrix(c(1 - m_A [i], 0,
                                      0, 1 - m_B [i]), byrow = T, ncol = 2)
        Q [[i]] <- list()
        
        for(j in 1:length(z_A))
        {
            if (i == j)
            {
                Q [[i]] [[j]] <- diag(2)
            }
            else
                {
                    Q [[i]] [[j]] <- matrix(c(delQ.A[i, j], 0,
                                              0, delQ.B[i, j]), byrow = T, ncol = 2)
                }
        }
    }
    J <- matrix(0, nrow = 2 * length(z_A), ncol = 2 * length(z_A))

    for (i in 1:length(z_A))
    {
        for (j in 1:length(z_A))
        {
### block diagonal elements
            if(i == j)
            {
                J[2*i + c(-1, 0), 2*i + c(-1, 0)] <- (1 - h[i]) * M [[i]]
            }
### off block diagonal 
            else
            {
                J[2*i + c(-1, 0), 2*j + c(-1, 0)] <- h[i] * M[[i]] %*% P [[i]] %*% Q [[i]] [[j]]
            }
        }
    }
    even <- seq(2, length(z_A)*2, 2)
    odd <- even - 1
    
    return(J[c(odd, even), c(odd, even)])
}
