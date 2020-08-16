#'The assess_direct() function
#'@description  The assess_direct function is designed to 
#'evaluate the prediction accuracy of a Gaussian Graphical model(GGM) 
#'comparing with the true graph structure with a known precision matrix.
#'
#'@param PREC_for_graph It is the known precision matrix which is used to
#'assess the estimated precision matrix from GGM.
#'
#'@param OMEGA_for_graph It is the estimated precision matrix from a GGM.
#'
#'@param p It is an integer representing the number of dimension of both the
#'known and estimated precision matrix.
#'
#'@return  Return the list of assessment results including sensitivity,
#'specificity, NPV(test negative), PPV(test positive), true positive
#'and true negative. 
#'
#'
#'@examples
#'prec1 <- matrix(c(0,2,3,1,0,0.5,0,0,0.4),nrow=3,ncol=3)
#'prec0 <- matrix(c(0,1,2,1,0.5,0.2,0,1,1),nrow=3,ncol=3)
#'
#'assessresult <- assess_direct(prec1,prec0,p=3)
#'
#'@export


assess_direct <- function(PREC_for_graph, OMEGA_for_graph, p) {
    truepos = trueneg = 0
    testpos = testneg = 0
    corrpos = corrneg = 0
    
    PREC_vals <- PREC_for_graph[upper.tri(PREC_for_graph)]
    OMEGA_vals <- OMEGA_for_graph[upper.tri(OMEGA_for_graph)]
    
    ## add up TRUE's of various conditionals
    truepos <- sum(PREC_vals != 0)
    testpos <- sum(abs(OMEGA_vals) > 0)
    corrpos <- sum((PREC_vals != 0) & (abs(OMEGA_vals) != 
        0))
    trueneg <- sum((PREC_vals == 0))
    testneg <- sum((OMEGA_vals == 0))
    corrneg <- sum((PREC_vals == 0) & (OMEGA_vals == 
        0))
    
    sensitivity = corrpos/truepos
    specificity = corrneg/trueneg
    NPV = corrneg/testneg
    PPV = corrpos/testpos
    
    return(list(truepos = truepos * 2, trueneg = trueneg * 
        2, sensitivity = sensitivity, specificity = specificity, 
        NPV = NPV, PPV = PPV))
}
