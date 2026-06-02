#' BEMA for the standard spiked covariance model
#'
#' Apply BEMA algorithm for spiked covariance proposed in the paper "Estimation of the number of spiked eigenvalues in a covariance matrix by bulk eigenvalue matching analysis""
#'
#' @param eigenvalue a list of eigenvalues to choose from
#' @param p dimension of the features
#' @param n number of samples
#' @param alpha a tuning parameter in the analysis, a default value is set to 0.2
#' @param beta a tuning parameter on computing quantile in Tracy-Widom, a default value is set to be 0.1
#'
#' @importFrom RMTstat qmp qtw
#'
#' @return The total number of spikes extracted, K
#'
#' @keywords Spike numbers
#'
#' @examples
#' x = matrix(rnorm(1000, 100), nrow = 1000)
#' eigen_x = svd(x)
#' eigen_out = list(eigenvalue = eigen_x$d^2 / 100, p = 1000, n = 100)
#' BEMA(eigen_out$eigenvalue, p = 1000, n = 100)
#'
#' @export

BEMA <- function(eigenvalue, p, n, alpha = 0.2, beta = 0.1){
    # p = min(n, p)
    # n = max(n, p)

    p_tilde = min(n, p)
    gamma_n = p / n
    k = seq(ceiling(p_tilde * alpha), floor(p_tilde * (1 - alpha)), 1)
    # print(k / p_tilde)
    k_quantile = marchenko_pastur_quantile(k / p_tilde, n = n, p = p)

    sigma = sum(k_quantile * eigenvalue[k]) / sum(k_quantile * k_quantile)
    # print(sigma)
    K = sum(eigenvalue > sigma * ((1 + sqrt(gamma_n))^2 + qtw(1 - beta) * n^(-2/3) * gamma_n^(-1/6) * (1 + sqrt(gamma_n))^(4/3)))

    if(is.na(K)){
        return(0)
    }else{
        return(K)
    }
}
