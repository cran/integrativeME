# Copyright (C) 2009 
# Kim-Anh Lê Cao, ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
# and Queensland Facility for Advanced Bioinformatics, The University of Queensland, Australia
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.





# -- -----------------to count results ----------------------
counts = function(x, type, n){
	err = sum(x != type, na.rm = TRUE)/(n - sum(is.na(x)))
	return(err)
}	

# ---------------------integrativeME main function ----------

integrativeME = function(
	data.cat,
	data.cont,
	type,
	select = c('RF', 'student', 'sPLS'),
	method = c('logreg', 'indep', 'loc'),
	loc.ind =  NULL,
	keepX = 5,
	ng = 2,
	mode.sPLS = NULL,
	fold = 10,
	kcv =1                     # number of 10-fold cv    
){


# ------- put warnings --------------

if (length(dim(data.cont)) != 2)  stop("'data.cont' must be a numeric matrix.")
if (length(dim(data.cat)) != 2)  stop("'data.cat' must be a numeric matrix.")


if(length(type) != nrow(data.cont)) stop('vector type does not match the number of rows of data.cont')
if(length(type) != nrow(data.cat)) stop('vector type does not match the number of rows of data.cat')

if(any(data.cat == 0)) stop('0 values are not allowed in the categorical data')

if(any(is.na(data.cat))) stop('NA values in the categorical data')
if(any(is.na(data.cat))) stop('NA values in the continuous data')

if(select == 'sPLS' && is.null(mode.sPLS)) stop('mode.sPLS needed with sPLS selection')

if(method == 'loc' && is.null(loc.ind)) stop('location index missing')

# ---- unmap the data to apply sPLS and select variables------
if(select == 'sPLS' && mode.sPLS == 'canonical'){
nleveltotal=0
for(i in 1:ncol(data.cat)){
	nleveltotal = nleveltotal + nlevels(as.factor(data.cat[,i]))
}
nleveltotal  

# use library mclust to unmap the data
i=1
var.cat = data.cat[,i]
data.unmap = unmap(var.cat)

for(i in 2:ncol(data.cat)){
	var.cat = data.cat[,i]
	data.unmap = cbind(data.unmap,unmap(var.cat))
}
}  # end if spls


# ------------------- MAIN----------- -----------------

n = nrow(data.cont)             #number of samples (total)
nv = keepX                      #number of continuous var
nq = ncol(data.cat)             #number of discrete var


mat.predicted = matrix(data = NA, nrow = n, ncol = kcv)  

# -------------rerun 10 fold cross CV kcv TIMES ------------------------------
ktest=1
for(ktest in ktest:kcv){
cat('run number:', ktest, '\n')


# ---------------- sample the cross validation sets -----------------------
cv.again=TRUE                     # change here
while(cv.again==TRUE){
crossval  = sample(1:fold, n, replace=TRUE)

if(any(summary(as.factor(crossval)) <= 2)){cv.again=TRUE; cat('seems that variable fold is too large for the data set, we resample again', '\n'); break} else{cv.again=FALSE}
}


predicted = c(rep(NA,n))
mat.gum = matrix(nrow = n, ncol =ng)


# ======================= 10 CV begins =====================================

##continue = TRUE
for(jcross in 1:fold){  ;
cat('cross validation sample number:', jcross, '\n')
train = which(crossval != jcross)
test = which(crossval == jcross)

indclass = sum(type[train] == 0)      # number of indiv ranked in class 1 (first class)

X = data.cont[train,]

cat('selecting genes ...', '\n')


#-------------------------select with sPLS regression or canonical, only regression for cont! --------------
if(select == 'sPLS'){

if(mode.sPLS == 'regression') {spls.res = spls(X, type[train], ncomp=1, keepX = keepX, mode = mode.sPLS)} else {spls.res = spls(X, Y = data.unmap[train,], ncomp=1, keepX = keepX, mode = mode.sPLS, scaleY = FALSE)}

list.gene = names(which(spls.res$loadings$X[,1] != 0))

}


#-------------------------select with RF --------------
if(select == 'RF'){
library(randomForest)

rf = randomForest(X, as.factor(type[train]), ntree=1000, importance=TRUE)
list.gene = names(sort(rf$importance[,3], decreasing=TRUE)[1:keepX])   #accurancy imp
}

#-------------------------select with student --------------
if(select == 'student'){
pval=apply(t(X),1,function(x){oneway.test(x~as.factor(type[train]),var.equal=TRUE)$p.value})
names(pval) = colnames(X)

list.gene = names(sort(pval, decreasing=FALSE)[1:keepX])   #accurancy imp
}


# ------------- initialize with k-means ----------------
cat('kmeans \n')
res.kmeans = kmeans.init(data.cont[train, list.gene])
if(res.kmeans$continue ==TRUE){    # condition coming from K-means clusters


# ------------ choose the ME gating function ---------
cat('ngme full \n')
if(method == 'indep') res.ME = MEindep(
	jcross, 
	train, 
	test, 
	n, 
	nq, 
	nv, 
	ng, 
	indclass, 
	data.cat, 
	data.cont = data.cont[, list.gene], 
	prop.kmeans = res.kmeans$prop.kmeans, 
	means.kmeans = res.kmeans$means.kmeans, 
	var.kmeans = res.kmeans$var.kmeans)


if(method == 'logreg') res.ME = MElogreg(
	jcross, 
	train, 
	test, 
	n, 
	nq, 
	nv, 
	ng, 
	indclass, 
	data.cat, 
	data.cont = data.cont[, list.gene], 
	prop.kmeans = res.kmeans$prop.kmeans)

if(method == 'loc') res.ME = MEloc(
	jcross, 
	train, 
	test, 
	n, 
	nq, 
	nv, 
	ng, 
	indclass, 
	data.cat, 
	data.cont = data.cont[, list.gene], 
	prop.kmeans = res.kmeans$prop.kmeans, 
	means.kmeans = res.kmeans$means.kmeans, 
	var.kmeans = res.kmeans$var.kmeans,
	loc.ind)

#if(method == 'cont') source('NGMEcont.R')
} else {ktest=ktest-1; jcross=nlevels(crossval) +1}


# ---------- store the predictions -------------------------

	predicted[which(res.ME$mat.gum[,2]>=res.ME$mat.gum[,1])] = 1   # NORM W
	predicted[which(res.ME$mat.gum[,2]<res.ME$mat.gum[,1])] = 0    #NORM W

	if(res.kmeans$continue == TRUE){mat.predicted[test,ktest] = predicted[test]}else{ktest = ktest-1}

}   #end crossval


}  # end ktest



return(invisible(list(
mean.error = apply(mat.predicted, 2, counts, type = type, n=n),
mat.predicted = mat.predicted
)))


}



