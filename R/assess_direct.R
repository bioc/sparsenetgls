#'The assess_direct() function is designed to 
#'produce the prediction accuracy of a Gaussian Graphical model(GGM) 
#'to the true graph structure with a known precision matrix.
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


assess_direct <- function(PREC_for_graph,OMEGA_for_graph,p)

{
    truepos=trueneg=0
    testpos=testneg=0
    corrpos=corrneg=0
    
    for ( i in seq_len(p))
        for ( j in i:p)  {if ((i!=j) & (PREC_for_graph[i,j]!=0))
                            truepos=truepos+1}
    
    for ( i in seq_len(p))
        for ( j in i:p)  {if ((i!=j) & abs(OMEGA_for_graph[i,j])>0)
                            testpos=testpos+1}
    
    for ( i in seq_len(p))
        for ( j in i:p)  {if ((i!=j) & (PREC_for_graph[i,j]!=0) &
                                (OMEGA_for_graph[i,j]!=0)) corrpos=corrpos+1}
    
    #Counting the true negatives
    for ( i in seq_len(p))
        for ( j in i:p)  {if ((i!=j) & (PREC_for_graph[i,j]==0))
                            trueneg=trueneg+1}
    for ( i in seq_len(p))
        for ( j in i:p)  {if ((i!=j) & (OMEGA_for_graph[i,j]==0))
                            testneg=testneg+1}
    
    for ( i in seq_len(p))
        for ( j in i:p)  {if ((i!=j) & (PREC_for_graph[i,j]==0)
                            & (OMEGA_for_graph[i,j]==0))
                            corrneg=corrneg+1}
    
    sensitivity = corrpos/truepos
    specificity = corrneg/trueneg
    NPV = corrneg/testneg
    PPV = corrpos/testpos
    
    return(list(truepos=truepos*2,trueneg=trueneg*2,sensitivity=sensitivity,
                specificity=specificity,NPV=NPV,PPV=PPV))
}
