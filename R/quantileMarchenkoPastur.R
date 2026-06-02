#' Quantile for Marchenko Pastur Distribution
#'
#' Compute the q quantile for Marchenko Pastur Distribution with n and p
#'
#' @param q quantile to be chosen
#' @param n number of samples to be considered
#' @param p dimension of the covariance matrix
#' @param step_size step length to do numerical integration
#'
#' @return Quantile location
#'
#' @keywords Marchenko-Pastur, quantile
#'
#' @examples
#' out = marchenko_pastur_quantile(0.2, 50, 500)
#'
#' @export

marchenko_pastur_quantile <- function(q, n, p, step_size = 10e-5){
    lambda = p / n
    lambda_minus = (1 - sqrt(lambda))^2
    lambda_plus = (1 + sqrt(lambda))^2
    x = seq(lambda_minus, lambda_plus, step_size)
    y = 1 / (2 * pi * min(lambda,1) * x) * sqrt((lambda_plus - x) * (x - lambda_minus))
    cumulative_y = cumsum(y * step_size)

    out = c()
    for(i in 1 : length(q)){
        out = c(out, x[which.min(abs(1 - cumulative_y - q[i]))])
    }
    return(out)
}

