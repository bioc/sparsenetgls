#'The path_result_for_roc() function
#'@description The path_result_for_roc function is designed to produce the 
#'Reciever Operative Characteristics (ROC) Curve for visualizing the prediction
#'accuracy of a Gaussian Graphical model (GGM) to the true graph structure.
#'The GGM must use a l-p norm regularizations (p=1,2) with the series of
#'solutions conditional on the regularization parameter.
#'
#'@param PREC_for_graph It is the known precision matrix which is used to
#'assess the estimated precision matrix from GGM.
#'
#'@param OMEGA_path It is an estimated precision matrix from a GGM model
#'using the penalized path with a range of values
#'(i.e. \eqn{\lambda,\in [0,1]}).
#'
#'@param pathnumber It represents the number of values (i.e. \eqn{\lambda})
#'used in the penalized GGM. This will result in the number of co-ordinates
#'used to form the ROC curve.
#'
#'@return Return the list of assessment results for a serie of precision
#'matrices. The results include sensitiviy/specificity/NPV/PPV
#'
#'
#'@examples
#'prec1 <- matrix(c(0,2,3,1,0,0.5,0,0,0.4),nrow=3,ncol=3)
#'Omega_est <- array(dim=c(3,3,3))
#'Omega_est[,,1] <- matrix(c(0,1,2,1,0.5,0.2,0,1,1),nrow=3,ncol=3)
#'Omega_est[,,2] <- matrix(c(0,1,0,1,0.5,0.2,0,1,1),nrow=3,ncol=3)
#'Omega_est[,,3] <- matrix(c(0,1,0,1,0,0.2,0,1,1),nrow=3,ncol=3)
#'rocpath <- path_result_for_roc(PREC_for_graph=prec1,OMEGA_path=Omega_est,
#'pathnumber=3)
#'
#'@export

path_result_for_roc <- function(PREC_for_graph, OMEGA_path, 
    pathnumber) {
    result_assessment_path <- array(dim = c(pathnumber, 
        4))
    pNO <- dim(OMEGA_path)[1]
    
    for (k in seq_len(pathnumber)) {
        
        print(assess_direct(PREC_for_graph, OMEGA_path[, 
            , k], p = pNO)[3:6])
        result <- assess_direct(PREC_for_graph, OMEGA_path[, 
            , k], p = pNO)
        
        result_assessment_path[k, 1] <- result$sensitivity
        result_assessment_path[k, 2] <- result$specificity
        result_assessment_path[k, 3] <- result$NPV
        result_assessment_path[k, 4] <- result$PPV
    }
    return(result_assessment_path)
}
