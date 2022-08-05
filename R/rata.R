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

    if ((!is.null(attr(tform,"specials")$present)) || (!is.null(attr(tform,"specials")$value))){
        whichlong<-c(attr(tform,"specials")$present,attr(tform,"specials")$value)
        v<-all.vars(attr(tform,"variables")[-1][whichlong])
        if (length(v)>1) stop("There can be only one")
        d<-long_expand(mf[,whichlong[1]], v)
        longmf<-mf[d[[1]],]
        pos<-match(names(mf),names(d))
        longmf[,whichlong]<-d[,pos[!is.na(pos)]]
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
        longmf<-mf[longmr$id,]
        longmf[,whichlong]<-longmr[,!(names(longmr) %in% "id")]
        id2<-longmr$id
        mf<-longmf
        if (!is.null(id1)){
          ## id merge thing somehow
            stop("write id merge thing here")
        } else {
            id<-id2
        }
    } else {
        id<-id1
    }
    
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
        data = data, offset = offset, control = control, method = method, 
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
