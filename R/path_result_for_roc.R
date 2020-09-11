#'The path_result_for_roc() function
#'@description The path_result_for_roc function is designed to evaluate the 
#'the prediction accuracy of a series Gaussian Graphical models (GGM) comparing 
#'to the true graph structure.
#'The GGM must use a l-p norm regularizations (p=1,2) with the series of
#'solutions conditional on the regularization parameter.
#'
#'@param PREC_for_graph It is the known precision matrix which is used to
#'assess the estimated precision matrix from GGM.
#'
#'@param OMEGA_path It is a matrix comprising of a series estimated precision 
#'matrices from a GGM model using a penalized path based on a range of structure
#'parameters (i.e. \eqn{\lambda,\in [0,1]}).
#'
#'@param pathnumber It represents the number of graph models 
#'(i.e. \eqn{\lambda}) for the evaluation.The value of pathnumber can be the 
#'same number used in a penalized path.  
#'
#'@return Return the list of assessment results for a series of precision
#'matrices. The results include sensitivity/specificity/NPV/PPV
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
