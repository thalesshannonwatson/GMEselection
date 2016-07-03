#' A variable selection procedure for regression models based on Generalized Maximum Entropy estimation (A. Golan, J. Judge and D. Miller, Maximum Entropy Econometrics, Wiley, 1996; Chapter 10). The procedure is useful for well- and ill-posed regression models, namely in models exhibiting small sample sizes, collinearity and non-normal errors.
#' 
#' @param y -- vector -- of size k.
#' @param x -- matrix of size (n,k) -- n samples of size k; The supports for all unknown parameters should be symmetric and') uniformly distributed around zero.
#' @param ssi (resp) -- TRUE/FALSE -- use the same support interval for all the unknown.
#' @param si (int1) -- c(float,float) -- support interval, in the form c(-|lower|,upper), common for all unknown parameters.
#' @param m -- integer -- the number of points in each parameter support. Usually the estimation is performed with five points in the parameter supports. Naturally, you can define a higher value.
#' @param intervals -- list of lists -- ex. list(c(-l1,s1), c(-l2,s2),...), list of pairs (.,.) constituting a different support for each parameter. Note that, in this model, you must specify the same number of pairs as size of vector y.
#' @param es (resp1) -- TRUE/FALSE -- error supports using an estimate of the error standard deviation from the OLS residuals.
#' @param error.csi (int2) -- c(-e,e) -- the support interval, [-e,e], for the error component.
#' @param j -- integer -- number of points in each error support. Usually the estimation is performed with three points in the error supports. Naturally, you can define a higher value.
#' @return 
#' \code{b} -- vector -- estimate of the unknown parameters;
#' \code{nepk} -- float -- normalized entropy for the intercept (if it exists) and for each variable;
#' \code{nep} -- vector -- normalized entropy for the signal.

# @examples
# GMEselection(1, 1)
# add(10, 1)
# Reference:
# Golan, A., Judge, G. and Miller, D. (1996). Maximum Entropy Econometrics: 
# Robust Estimation With Limited Data. Wiley, Chichester. (Chapter 10).


#' @export
GMEselection <- function(y,x,ssi=TRUE,si=c(-1,1),intervals=NULL,m=5,es=TRUE,error.csi=c(-1,1),j=3) {

    if (!requireNamespace("nloptr", quietly = TRUE)) {
        stop("nloptr package needed for this function to work. Please install it.",call. = FALSE)
    }

    # The supports for all unknown parameters should be symmetric and') uniformly distributed around zero.
    n <- nrow(x)
    k <- ncol(x)
    #cat("a matriz x tem dimensao:",n,k,"\n")

    #Verify:
    stopifnot( ! missing(x) )
    stopifnot( ! missing(y) )
    stopifnot( ! missing(si) )
    
    # m > 1
    # n > 
    # k >

    #disp('Do you want to specify the same support interval for all the unknown')
    #resp=input('parameters of the model (yes/no)? [yes] ','s');

    #TODO: what is Z ?
    Z <- matrix(0,k,k*m) # = zeros(k,k*m)

    if (ssi) {  
        #case: same support interval for all the unknown

        inc  <- (si[2]-si[1])/(m-1)
        s1   <- seq(from=si[1], to=si[2], by=inc)
        for (i in 1:k) {
            pos <- (i-1)*m+1
            ##cat("pos=",pos,"m=",m,"s1=",s1,"\n")
            Z[i,pos:(pos+m-1)] <- s1
        }
        
    } else {
        #case: NOT the same support interval for all the unknown
    
        #TODO: verify
        #size of intervals list is k
        #num_of_intervals <- len(intervals)
        #assert( num_of_intervals == k )

        s1 <- matrix(0, k, m)
        for (i in 1:k) {
            inc <- (intervals[i][2]-intervals[i][1])/(m-1)
            s1[i,1:m] <- seq(from=intervals[i][1],to=intervals[i][2],by=inc)
        }
        for (i in 1:k) {
            pos <- (i-1)*m+1
            Z[i,pos:(pos+m-1)] <- s1[i,1:m]
        }
    }

    sdy <- sd(y) #sdy is float, the standard deviation of vector y

    st1=floor(-3*sdy); if (st1==0) st1=-1
    st2=ceiling(3*sdy); if (st2==0) st2=1
    st3=floor(-4*sdy); if (st3==0) st3=-1
    st4=ceiling(4*sdy); if (st4==0) st4=1

    #QR decomposition here ???
    #assume x is NOT sparse (TODO: what could be this in R matrix ?)
    #write.table(x)
    res <- qr(x) #upper.tri(x, diag = TRUE) #TODO: 
    res$qr[ lower.tri( res$qr ) ] <- 0
    R <- res$qr[1:k,1:k]
    #cat("QR\n")
    #write.table(R)
    ##write.table(R)

    #if issparse(X)
    #    R = qr(X);
    #else
    #    R=triu(qr(X)); 
    #end
    #b_ols=R\(R'\(X'*Y));


    #disp('INFORMATION: the supports for the error component are usually defined by');
    #disp('the 3-sigma or 4-sigma rules, with sigma being the standard deviation of')
    #disp('the noisy observations, or an estimate of the error standard deviation');
    #disp('from the OLS estimation.');

    #ok
    #bb <- t(x) %*% y
    ##write.table(bb)


    #b_ols=R\(R'\(X'*Y)); %
    b_ols <- qr.solve(R, solve(R, t(x) %*% y))

    #msres=((Y-X*b_ols)'*(Y-X*b_ols))/(n-k);
    msres <- ( t(y-x %*% b_ols) %*% (y-x%*%b_ols)) / (n-k)

    st1a <- floor(-3*sqrt(msres)); if (st1a==0) st1a=-1
    st2a <- ceiling(3*sqrt(msres)); if (st2a==0) st2a=1
    st3a <- floor(-4*sqrt(msres)); if (st3a==0) st3a=-1
    st4a <- ceiling(4*sqrt(msres)); if (st4a==0) st4a=1



    #TODO: avisar o user da escolha?
    #TODO: como e' que resp1 afeta as contas seguintes?
#if isempty(resp1)
#    disp('The error supports can be defined by:')
#    disp(sprintf('[%d,%d] (using the 3-sigma rule);',st1a,st2a))
#    disp(sprintf('[%d,%d] (using the 4-sigma rule).',st3a,st4a))
#else
#    disp('Alternatively, using the standard deviation of the noisy observations,')
#    disp('the error supports can be defined by:')
#    disp(sprintf('[%d,%d] (using the 3-sigma rule);',st1,st2))
#    disp(sprintf('[%d,%d] (using the 4-sigma rule).',st3,st4))
#end



    inc2 <- (error.csi[2] - error.csi[1])/(j-1)
    s2   <- seq(from=error.csi[1],to=error.csi[2],by=inc2)
    V    <- matrix(0,n,n*j)
    for (i in 1:n) {
       pos <- (i-1)*j+1
       V[i,pos:(pos+j-1)] <- s2
    }


    p  <- (1/m)*matrix(1,k*m,1)
    w  <- (1/j)*matrix(1,n*j,1)
    pw <- cbind( t(p),t(w) )
    dp <- length(p)
    lb <- matrix(1e-10,nrow(pw),ncol(pw))
    ub <- matrix(1,nrow(pw),ncol(pw))

    #matrixadditivity1=kron(eye(k,k),ones(1,m));
    matrixadditivity1 <- kronecker(diag(k),matrix(1,1,m))

    #matrixadditivity2=kron(eye(n,n),ones(1,j));
    matrixadditivity2 <- kronecker(diag(n),matrix(1,1,j))

    #Aeq=[X*Z,V;matrixadditivity1,zeros(k,n*j);zeros(n,k*m),matrixadditivity2];
    # x dims are: n, k
    # Z dims are: k, k*m
    # x * Z  dims are: n, k*m
    #cat("x%*%Z=\n")
    #write.table(x%*%Z)
    # V dims are: n, n*j
    #cat("V\n")
    #write.table(V)
    # cbind(x%*%Z,V) are: n, k*m + n*j
    #cat("cbind(x%*%Z,V)\n")
    #write.table(cbind(x%*%Z,V))

    Aeq <- rbind( cbind(x%*%Z,V),  #dims = n,  k*m + n*j
                  cbind(matrixadditivity1,matrix(0,k,n*j)), # dims = k*1, k*m + n*j
                  cbind(matrix(0,n,k*m),matrixadditivity2)) # dims = n*1, n*j + k*m

    #Aeq dims are n+k+n , n*j + k*m

    #Beq=[Y;ones(n+k,1)];
    #cat("y\n")
    #write.table(y) #dims: n,1
    Beq=rbind( matrix(y,n,1), matrix(1,n+k,1) )
    #cat("Beq dims are ",nrow(Beq),ncol(Beq),"\n")
    #write.table(Beq)


    #TODO: 
    #%sempre a mesma FunGMES
    #o comando seguinte resolve:  Aeq x = Beq sendo lb <= x <= ub
    # nonlinear inequalities (dp is a vector)
    # x = fmincon(fun,      x0, A,  b,  Aeq, beq, lb, ub, nonlcon,options)
    # a = fmincon('FunGMES',pw, [], [], Aeq, Beq, lb, ub, [], [], dp);

    #Package nlopt
    #fmincon equivalente no R:
    #   http://www.ucl.ac.uk/~uctpjyy/nloptr.html
    # nloptr (09/11/2013 Version 0.9.4 updated NLopt to version 2.4.) does not work on R 3.3.1!
    #WORKS BUT has to be downloaded from britol-uk
    # O interface nlopt/R parou em 2014 mas o package nlopt ainda esta' em desenvolvimento.
    #
    # Steven G. Johnson, The NLopt nonlinear-optimization package, http://ab-initio.mit.edu/nlopt
    #

    #"optim" function in R
    # https://stat.ethz.ch/R-manual/R-patched/library/stats/html/optim.html

    #"optimx" (replacement of optim)
    # https://cran.r-project.org/web/packages/optimx/index.html

    #Optimization in R
    #https://cran.r-project.org/web/views/Optimization.html

    #Classifi#cation by Subject
    #https://cran.r-project.org/web/views/Optimization.html#Classifi#cationBySubject

    #cat("main dp=",dp,"\n")

    # Objective funtion
    eval_f <- function(pw) {  #this is FunGMES

        #TODO: 
        # dp is in the scope of caller function: must declare here?

        #cat("\npw=",pw,"\npw has length ",length(pw),"\n")
        #cat("local dp=",dp,"\n")

        #p=pw(1:dp)';
        p <- t(pw[1:dp])
        #cat("p=",p,"\n")

        #cat("t(p) %*% log(p)=",t(p) %*% log(p),"\n")

        #w=pw(dp+1:end)'
        w   <- t(matrix( pw[(dp+1):length(pw)] ))
        #cat("\nw=",w,"\n")

        #cat("dim p=",dim(p),"\ndim t(log(p))=",dim(t(log(p))),"\n")
        #cat("dim w=",dim(w),"\ndim t(log(w))=",dim(t(log(w))),"\n")

        f   <- -(-p %*% t(log(p)) - w %*% t(log(w)))

        df <- cbind( log(p), log(w)) + 1

        return( list( "objective" = f, "gradient" = df) )
    }


    # Equality constraints    
    eval_g_eq <- function(pw) {

        #cat("dim Aeq=",dim(Aeq),"\n")
        #cat("dim pw=",dim(pw),"\n")

        #cat("pw=",pw,"\n")

        constr <- Aeq %*% matrix(pw) - Beq 

        return( list( "constraints"=constr, "jacobian"=Aeq ) )
    }

 
    # Set optimization options.
    # Options for NLOPT_LD_AUGLAG
    local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                        "xtol_rel"  = 1.0e-7 )
    opts <- list( "algorithm"   = "NLOPT_LD_AUGLAG",
                  "xtol_rel"    = 1.0e-7,
                  "maxeval"     = 1000,
                  "local_opts"  = local_opts,
                  "print_level" = 0 )
    
    #opts = list( "algorithm" = "NLOPT_LN_COBYLA",
    #             "xtol_rel"  = 1.0e-8)


    # Do optimization.
    res <- nloptr( x0          = pw, 
                   eval_f      = eval_f, 
                   lb          = lb, 
                   ub          = ub, 
                   #eval_g_ineq = eval_g_ineq, 
                   eval_g_eq   = eval_g_eq, 
                   opts        = opts )


    a <- res$solution

    p <- a[1:dp]
    b <- Z %*% p

    nepk <- matrix(0,k,1)
    for (i in 1:k) {
        pos <- (i-1)*m+1
        nepk[i,1] <- -t(p[pos:pos+m-1]) %*% log(p[pos:pos+m-1])/log(m)
    }
    nepk <- t(nepk)
    nep  <- (-t(p)%*%log(p))/(k*log(m))

    return( list( "b" = b, "nepk" = nepk, "nep" = nep ) )

} #end of GMEselection function



