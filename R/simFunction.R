#' Simulated Component Rotation
#'
#' Rotated the simulated component generated
#'
#' @param comp A matrix of components
#' @param angle Rotation angle of the component
#'
#' @return Matrix of rotated component
#'
#' @keywords rotation, component
#'
#' @examples
#' component = svd(matrix(rnorm(100 * 200), nrow = 200))$u[, 1 : 2]
#' rotate_component(component, pi / 6)
#'
#' @export

rotate_component <- function(comp, angle = 0){
    rot_mat = matrix(c(cos(angle), sin(angle), -1 * sin(angle), cos(angle)), nrow = 2, ncol = 2)
    return(comp %*% rot_mat)
}

#' Random Score Generation
#'
#' Gererate random scores for the simulation
#'
#' @param dim_score dimension of the scores
#' @param num_score number of scores
#' @param score_variance variance of the score
#'
#' @importFrom MASS mvrnorm
#'
#' @return a matrix of generated random scores, with dim_score * num_score
#'
#' @keywords random score
#'
#' @examples
#' score_generation(2, 10, c(1,2))
#'
#' @export

score_generation <- function(dim_score, num_score, score_variance){
    sim_cov = diag(score_variance)
    sim_mean = rep(0, dim_score)
    return(t(mvrnorm(num_score, sim_mean, sim_cov)))
}

#' Configuration for simulated data generation
#'
#' generate the configuration for the data
#'
#' @param featureNum number of features
#' @param DataNum number of data
#' @param commonlySharedNum number of common component
#' @param partiallySharedNum number of partial shared component
#' @param individualSharedNum number of individual component
#' @param noiseVariance variance of noise
#'
#' @keywords configure the setting for the generation of data
#'
#' @examples
#' configuration_setting_generation()
#'
#' @export

configuration_setting_generation <- function(featureNum = 50,
                                             DataNum = c(100, 100, 100, 100),
                                             commonlySharedNum = 2,
                                             partiallySharedNum = c(2,2,2,2),
                                             individualSharedNum = c(2,2,2,2),
                                             noiseVariance = c(1,1,1,1)){

    ## SVD of random matrix
    tempRandomMatrix = matrix(rnorm(featureNum * (sum(commonlySharedNum, partiallySharedNum, individualSharedNum))), nrow = featureNum)
    svdRandomMatrix = svd(tempRandomMatrix)

    ## Common Component
    commonComponent = svdRandomMatrix$u[, 1 : commonlySharedNum]

    ## Partially shared component
    partialComponent = list()
    partialComponent[[1]] = svdRandomMatrix$u[, (commonlySharedNum + 1) : (commonlySharedNum + partiallySharedNum[1])]
    partialComponent[[2]] = svdRandomMatrix$u[, (commonlySharedNum + partiallySharedNum[1] + 1) : (commonlySharedNum + partiallySharedNum[1] + partiallySharedNum[2])]
    partialComponent[[3]] = svdRandomMatrix$u[, (commonlySharedNum + partiallySharedNum[1] + partiallySharedNum[2] + 1) : (commonlySharedNum + partiallySharedNum[1] + partiallySharedNum[2] + partiallySharedNum[3])]
    partialComponent[[4]] = svdRandomMatrix$u[, (commonlySharedNum + partiallySharedNum[1] + partiallySharedNum[2] + partiallySharedNum[3] + 1) : (commonlySharedNum + partiallySharedNum[1] + partiallySharedNum[2] + partiallySharedNum[3] + partiallySharedNum[4])]

    ## Individual Component
    individualComponent = list()
    individualComponent[[1]] = svdRandomMatrix$u[, (sum(commonlySharedNum, partiallySharedNum) + 1) : (sum(commonlySharedNum, partiallySharedNum) + individualSharedNum[1])]
    individualComponent[[2]] = svdRandomMatrix$u[, (sum(commonlySharedNum, partiallySharedNum) + individualSharedNum[1] + 1) : (sum(commonlySharedNum, partiallySharedNum) + individualSharedNum[1] + individualSharedNum[2])]
    individualComponent[[3]] = svdRandomMatrix$u[, (sum(commonlySharedNum, partiallySharedNum) + individualSharedNum[1]+ individualSharedNum[2] + 1 ): (sum(commonlySharedNum, partiallySharedNum) + individualSharedNum[1] + individualSharedNum[2] + individualSharedNum[3])]
    individualComponent[[4]] = svdRandomMatrix$u[, (sum(commonlySharedNum, partiallySharedNum) + individualSharedNum[1]+ individualSharedNum[2] + individualSharedNum[3] + 1 ): (sum(commonlySharedNum, partiallySharedNum) + individualSharedNum[1] + individualSharedNum[2] + individualSharedNum[3] + individualSharedNum[4])]


    return(list(featureNum = featureNum,
                DataNum = DataNum,
                commonlySharedNum = commonlySharedNum,
                partiallySharedNum = partiallySharedNum,
                individualSharedNum = individualSharedNum,
                noiseVariance = noiseVariance,
                commonComponent = commonComponent,
                partialComponent = partialComponent,
                individualComponent = individualComponent))
}

#' Simulated Data Generation
#'
#' Generate simulation data
#'
#' @param configuration_setting setting for the configuration
#' @param amplitude The amplitude of the score variance
#' @param heterogeneousNoise Whether the noise for each dataset to be heteregeneous
#'
#' @keywords simulated data generation
#'
#' @examples
#' configuration_setting = configuration_setting_generation()
#' simulated_data_generation(configuration_setting)
#'
#' @export

simulated_data_generation <- function(configuration_setting,
                                      amplitude = 1,
                                      heterogeneousNoise = FALSE){
    ## Setting up configurations
    featureNum = configuration_setting$featureNum
    DataNum = configuration_setting$DataNum
    commonlySharedNum = configuration_setting$commonlySharedNum
    partiallySharedNum = configuration_setting$partiallySharedNum
    individualSharedNum = configuration_setting$individualSharedNum
    noiseVariance = configuration_setting$noiseVariance
    commonComponent = configuration_setting$commonComponent
    partialComponent = configuration_setting$partialComponent
    individualComponent = configuration_setting$individualComponent

    ## Rotate each component
    A11 = rotate_component(commonComponent, 0)
    A12 = rotate_component(commonComponent, pi / 6)
    A21 = rotate_component(commonComponent, pi / 3)
    A22 = rotate_component(commonComponent, pi / 2)

    B11 = rotate_component(partialComponent[[1]], 0)
    B12 = rotate_component(partialComponent[[1]], pi / 4)
    B21 = rotate_component(partialComponent[[2]], 0)
    B22 = rotate_component(partialComponent[[2]], pi / 4)

    C11 = rotate_component(partialComponent[[3]], 0)
    C12 = rotate_component(partialComponent[[4]], 0)
    C21 = rotate_component(partialComponent[[3]], pi / 4)
    C22 = rotate_component(partialComponent[[4]], pi / 4)

    D11 = rotate_component(individualComponent[[1]], 0)
    D12 = rotate_component(individualComponent[[2]], 0)
    D21 = rotate_component(individualComponent[[3]], 0)
    D22 = rotate_component(individualComponent[[4]], 0)

    F1 = score_generation(commonlySharedNum, DataNum[1], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[1])
    F2 = score_generation(commonlySharedNum, DataNum[2], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[2])
    F3 = score_generation(commonlySharedNum, DataNum[3], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[3])
    F4 = score_generation(commonlySharedNum, DataNum[4], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[4])

    G1 = score_generation(partiallySharedNum[1], DataNum[1], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[1])
    G2 = score_generation(partiallySharedNum[1], DataNum[2], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[2])
    G3 = score_generation(partiallySharedNum[2], DataNum[3], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[3])
    G4 = score_generation(partiallySharedNum[2], DataNum[4], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[4])

    H1 = score_generation(partiallySharedNum[3], DataNum[1], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[1])
    H2 = score_generation(partiallySharedNum[4], DataNum[2], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[2])
    H3 = score_generation(partiallySharedNum[3], DataNum[3], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[3])
    H4 = score_generation(partiallySharedNum[4], DataNum[4], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[4])

    K1 = score_generation(individualSharedNum[1], DataNum[1], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[1])
    K2 = score_generation(individualSharedNum[2], DataNum[2], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[2])
    K3 = score_generation(individualSharedNum[3], DataNum[3], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[3])
    K4 = score_generation(individualSharedNum[4], DataNum[4], c(runif(1, 1, 2), runif(1, 1, 2)) * amplitude * noiseVariance[4])

    E1 <- matrix(rnorm(featureNum * DataNum[1], 0, sqrt(noiseVariance[1])), nrow = featureNum)
    E2 <- matrix(rnorm(featureNum * DataNum[2], 0, sqrt(noiseVariance[2])), nrow = featureNum)
    E3 <- matrix(rnorm(featureNum * DataNum[3], 0, sqrt(noiseVariance[3])), nrow = featureNum)
    E4 <- matrix(rnorm(featureNum * DataNum[4], 0, sqrt(noiseVariance[4])), nrow = featureNum)

    if (heterogeneousNoise){
        E1 = diag(runif(featureNum, 0.5, 1)) %*% E1
        E2 = diag(runif(featureNum, 0.5, 1)) %*% E2
        E3 = diag(runif(featureNum, 0.5, 1)) %*% E3
        E4 = diag(runif(featureNum, 0.5, 1)) %*% E4
    }

    #
    Y1 = A11 %*% F1 + B11 %*% G1 + C11 %*% H1 + D11 %*% K1 + E1
    Y2 = A12 %*% F2 + B12 %*% G2 + C12 %*% H2 + D12 %*% K2 + E2
    Y3 = A21 %*% F3 + B21 %*% G3 + C21 %*% H3 + D21 %*% K3 + E3
    Y4 = A22 %*% F4 + B22 %*% G4 + C22 %*% H4 + D22 %*% K4 + E4

    return(list(Y1, Y2, Y3, Y4))
}
