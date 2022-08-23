## magic
present<-function(y) y
value<-function(x) x

each<-function(x) x

## ordinary functions
wide<-function(mr_object) as.matrix(mr_object)

has<-function(mr_object, value) mr_object %has% value
hasonly<-function(mr_object, value) mr_object %hasonly% value
hasany<-function(mr_object, value) mr_object %hasany% value
hasall<-function(mr_object, value) mr_object %hasall% value

long_expand<-function(mr_object,name){
    y<-as.matrix(as.mr(mr_object))
    Y<-as.vector(y)
    id<-rep(1:nrow(y), ncol(y))
    level<-rep(levels(mr_object),each=nrow(y))
    rval<-data.frame(id=id, value=level, present=Y)
   
    names(rval)[1]<-deparse(bquote(id(.(as.name(name)))))
    names(rval)[2]<-deparse(bquote(value(.(as.name(name)))))
    names(rval)[3]<-deparse(bquote(present(.(as.name(name)))))
    rval
}


vcov.mrglm<-function(model,...) model$vcov

mrglm_estfun<-function(model) model.matrix(model)*model$weights*resid(model,"working")

mrglm<-function (formula, family = gaussian, data, weights, subset, 
    na.action, start = NULL, etastart, mustart, offset, control = list(...), 
    model = TRUE,  x = FALSE, y = TRUE, singular.ok = TRUE, 
    contrasts = NULL, ...) 
{
    method<-"glm.fit"
    cal <- match.call()
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
         
    if (missing(data)) 
        stop("'data' must be specified")
    special <- c("each","value","present")
    tform <-  terms(formula, special, data = data)
    
    if (!(all(all.vars(formula) %in% names(data))))
        stop("all variables must be in data= argument")
    
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mfterms<-attr(mf,"terms")

    haspresent<-!is.null(attr(tform,"specials")$present)
    hasvalue<-!is.null(attr(tform,"specials")$value)

 
    
    if (haspresent || hasvalue ){
        whichlong<-c(attr(tform,"specials")$present,attr(tform,"specials")$value)
        v<-unique(unlist(lapply(attr(tform,"variables"), all.vars)[-1][whichlong]))
        if (length(v)>1) stop("Only one variable may have present()/value() terms")
        d<-long_expand(mf[,whichlong[1]], v)
        longmf<-mf[d[[1]],,drop=FALSE]
        pos<-match(names(mf),names(d))
        longmf[,whichlong]<-d[,pos[!is.na(pos)],drop=FALSE]
        id1<-d[[1]]
        mf<-longmf
    } else {
        id1<-NULL
    }

    
    if (!is.null(attr(tform,"specials")$each)){
        ## need to expand the data
        whichlong<-attr(tform,"specials")$each
         nms<-sapply(names(mf)[whichlong],as.name)
        names(nms)<-NULL
        longmr<- eval(bquote(with(mf, mr_stack(..(nms))), splice=TRUE))
        longmf<-mf[longmr$id,,drop=FALSE]
        longmf[,whichlong]<-longmr[,!(names(longmr) %in% "id"),drop=FALSE]
        id2<-longmr$id
        if (!is.null(id1)){
            id<-id1[id2]
            mf<-longmf
        } else {
            id<-id2
            mf<-longmf
        }
    } else {
        id<-id1
    }
    attr(mf,"terms")<-mfterms
    
    control <- do.call("glm.control", control)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    fit <- eval(call(if (is.function(method)) "method" else method, 
        x = X, y = Y, weights = weights, start = start, etastart = etastart, 
        mustart = mustart, offset = offset, family = family, 
        control = control, intercept = attr(mt, "intercept") > 
            0L, singular.ok = singular.ok))
    if (length(offset) && attr(mt, "intercept") > 0L) {
        fit2 <- eval(call(if (is.function(method)) "method" else method, 
            x = X[, "(Intercept)", drop = FALSE], y = Y, weights = weights, 
            offset = offset, family = family, control = control, 
            intercept = TRUE))
        if (!fit2$converged) 
            warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
        fit$null.deviance <- fit2$deviance
    }
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL

    
    rval<-structure(c(fit, list(call = cal, formula = formula, terms = mt, 
        data = data, offset = offset, control = control, method = method, id=id,
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf))),
        class = c(fit$class, c("mrglm", "glm", "lm")))

  
     infl<-mrglm_estfun(rval)%*%summary(rval)$cov.unscaled
     
     if (!is.null(id))
         V<-crossprod(rowsum(infl, id))  ##FIXME: more than one id to merge
     else
         V<-crossprod(infl)
     
     rval$naive.cov<-summary(rval)$cov.scaled
     rval$vcov<-V
     rval
}



mrloglin<-function(formula, data,...){

    if (!all(all.vars(formula) %in% names(data)))
        stop("all variables must be in data= argument")

    mf<-model.frame(formula, data)
    nms<-sapply(names(mf),as.name)
    names(nms)<-NULL
    longmf<-eval(bquote(with(mf, mr_stack(..(nms))),splice=TRUE))

    longdes<-survey::svydesign(id=~id, data=longmf, prob=~1)
    model<-survey::svyloglin(formula, design=longdes)

    model$call<-sys.call()
    class(model)<-c("mrloglin",class(model))
    model
}

anova.mrloglin<-function(object, object1, ..., integrate = FALSE) {
    rval<-NextMethod()
    class(rval)<-c("anova.mrloglin", class(rval))
    rval
}



anova.mrglm<-function (object, object2 = NULL, test = c("F", "Chisq"), method = c("LRT", 
    "Wald"), tolerance = 1e-05, ...) 
{
    test <- match.arg(test)
    method <- match.arg(method)
    if (is.null(object2)) 
        return(oneanova.mrglm(object, test, method))
    t1 <- attr(terms(object), "term.labels")
    t2 <- attr(terms(object2), "term.labels")

    X <- model.matrix(object)
    Z <- model.matrix(object2)
    if (nrow(X) != nrow(Z)) 
        stop("models have different numbers of observations")
    if (ncol(X) > ncol(Z)) {
        tmp <- X
        X <- Z
        Z <- tmp
        bigger <- 1
    }
    else bigger <- 2
    resids<-suppressWarnings(summary(lm(X ~ Z)))
    if (NCOL(X)==1){
        if (suppressWarnings(summary(lm(X ~ Z)))$sigma/(tolerance+SD(X))>tolerance)
            stop("models not nested")
    } else if (any(sapply(suppressWarnings(summary(lm(X ~ Z))), "[[", 
                          "sigma")/(tolerance + SD(X)) > tolerance)) {
        stop("models not nested")
    }
    XX <- matrix(nrow = nrow(Z), ncol = ncol(Z))
    xform <- lm(Z[, 1] ~ X + 0)
    XX[, 1] <- resid(xform)
    for (i in 2:ncol(Z)) {
        XX[, i] <- resid(xform <- lm(Z[, i] ~ X + Z[, 1:(i - 
            1)] + 0))
    }
    colkeep <- colMeans(abs(XX))/(tolerance + colMeans(abs(Z))) > 
        tolerance
    XX <- XX[, colkeep, drop = FALSE]
    index <- ncol(X) + (1:ncol(XX))
    mu <- if (bigger == 1) 
        fitted(object)
    else fitted(object2)
    eta <- if (bigger == 1) 
        object$linear.predictors
    else object2$linear.predictors
    offset <- if (bigger == 1) 
        object$offset
    else object2$offset
    if (is.null(offset)) 
        offset <- 0
    y <- object$y
    pweights <- rep(1, length(y))

    ywork <- eta - offset + (y - mu)/object$family$mu.eta(eta)
    wwork <- ((pweights * object$family$mu.eta(eta)^2)/object$family$variance(mu))
    wlm <- lm.wfit(cbind(X, XX), ywork, wwork)
    p1 <- 1:wlm$rank
    Ainv <- chol2inv(wlm$qr$qr[p1, p1, drop = FALSE])
    estfun <- cbind(X, XX) * wwork * ((y - mu)/object$family$mu.eta(eta))
    design <- object$survey.design

    V<-crossprod(rowsum(estfun%*%Ainv, object$id))
    
    V <- V[index, index]
    df <- min(object$df.residual, object2$df.residual)
    if (method == "LRT") {
        V0 <- Ainv[index, index]
        chisq <- if (bigger == 1) 
            deviance(object2) - deviance(object)
        else deviance(object) - deviance(object2)
        misspec <- eigen(solve(V0) %*% V, only.values = TRUE)$values
        if (test == "Chisq") 
            p <- survey::pchisqsum(chisq, rep(1, length(misspec)), misspec, 
                method = "sad", lower.tail = FALSE)
        else p <- survey::pFsum(chisq, rep(1, length(misspec)), misspec, 
            ddf = df, method = "sad", lower.tail = FALSE)
        rval <- list(call = sys.call(), chisq = chisq, df = length(index), 
            p = p, lambda = misspec, ddf = df, mcall = if (bigger == 
                1) object$call else object2$call, test.terms = if (bigger == 
                1) c(setdiff(t1, t2), "-", setdiff(t2, t1)) else c(setdiff(t2, 
                t1), "-", setdiff(t1, t2)))
        class(rval) <- "regTermTestLRT"
    }
    else {
        beta <- wlm$coefficients[index]
        chisq <- crossprod(beta, solve(V, beta))
        if (test == "Chisq") {
            p <- pchisq(chisq, df = length(index), lower.tail = FALSE)
        }
        else {
            p <- pf(chisq/length(index), df1 = length(index), 
                df2 = df, lower.tail = FALSE)
        }
        rval <- list(call = sys.call(), Ftest = chisq/length(index), 
                     df = length(index),
                     p = p, ddf = df,
                     mcall = if (bigger == 1) object$call else object2$call,
                     test.terms = if (bigger == 1) c(setdiff(t1, t2), "-", setdiff(t2, t1)) else c(setdiff(t2, t1), "-", setdiff(t1, t2)))
        class(rval) <- "regTermTest"
    }
    rval
}


oneanova.mrglm<-function (object, test, method) 
{
    tt <- terms(object)
    tlbls <- attr(tt, "term.labels")
    nt <- length(tlbls)
    if (nt < 2) 
        return(NULL)
    seqtests <- vector("list", nt)
    if (test == "F") 
        ddf <- NULL
    else ddf <- Inf
    lastmodel<-thismodel <- object
    n<-length(object$y)
    if (!("formula") %in% names(thismodel$call)) 
        names(thismodel$call)[[2]] <- "formula"
    thisformula<-formula(thismodel)
    for (i in nt:1) {
        thisterm <- tlbls[i]
        dropformula <- survey::make.formula(thisterm)[[2]]
        thisformula<-eval(bquote(update(thisformula, .~.-(.(dropformula)))))
        thismodel<- eval(bquote(update(thismodel, .(thisformula))))
        seqtests[[i]] <-  anova(thismodel, lastmodel)
        lastmodel<-thismodel
        if (length(thismodel$y)!=n) stop("Data sets are not the same size: missing values?")
    }
    class(seqtests) <- "seqanova.mrglm"
    attr(seqtests, "method") <- method
    attr(seqtests, "test") <- test
    seqtests
}

print.seqanova.mrglm<-function (x, ...) 
{
    isWald <- attr(x, "method") == "Wald"
    isF <- attr(x, "test") == "F"
    cat("Anova table: ")
    if (isWald) 
        cat("(Wald tests)\n")
    else cat(" (Rao-Scott LRT)\n")
    print(x[[length(x)]]$mcall)
    terms <- sapply(x, "[[", "test.terms")
    stats <- if (isF && isWald) 
        sapply(x, "[[", "Ftest")
    else sapply(x, "[[", "chisq")
    if (!isWald) 
        stats <- cbind(stats, DEff = sapply(x, function(xi) mean(xi$lambda)))
    df <- sapply(x, "[[", "df")
    p <- sapply(x, "[[", "p")
    if (!isF) {
        rval <- cbind(stats, df, p)
    }
    else {
        ddf <- sapply(x, "[[", "ddf")
        rval <- cbind(stats, df, ddf, p)
    }
    rownames(rval) <- terms[1,]
    printCoefmat(rval, tst.ind = 1, zap.ind = 2:3, has.Pvalue = TRUE)
    invisible(x)
}


SD<-function (x)  if (NCOL(x) > 1) apply(x, 2, sd) else sd(x)






mrmultinom<-function(formula, data, family=multinomial(),...){

    if (!all(all.vars(formula) %in% names(data)))
        stop("all variables must be in data= argument")

    mf<-model.frame(formula, data)
    nms<-sapply(names(mf),as.name)
    names(nms)<-NULL
    longmf<-eval(bquote(with(mf, mr_stack(..(nms))),splice=TRUE))

    longdes<-survey::svydesign(id=~id, data=longmf, prob=~1)
    model<-svyVGAM::svy_vglm(formula, design=longdes,family=family)

    model$call<-sys.call()
    class(model)<-c("mrmultinom",class(model))
    model
}

