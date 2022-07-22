long<-function(x) x
wide<-function(mr_object) as.matrix(mr_object)
has<-function(mr_object, value) mr_object %has% value
hasonly<-function(mr_object, value) mr_object %hasonly% value
hasany<-function(mr_object, value) mr_object %hasany% value
hasall<-function(mr_object, value) mr_object %hasall% value

mrglm<-function(formula, data, family,...){
 	
    special <- c("long")
    tform <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
    
    if (!(all(all.vars(formula) %in% names(data))))
	    stop("all variables must be in data= argument")

	if (!is.null(attr(tform,"specials")$long)){
		## need to expand the data
		whichlong<-attr(tform,"specials")$long
		mf<-model.frame(formula, data)
		nms<-sapply(names(mf)[whichlong],as.name)
		longmr<- eval(bquote(with(mf, mr_stack(..(nms))), splice=TRUE))
		longmf<-mf[longmr$id,]
		longmf[,whichlong]<-longmr[,-1]
	}



}


mrglm_estfun<-function(model) model.matrix(model)*model$weights*resid(model,"working")
