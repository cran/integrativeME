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

#===================================================================================
MEindep = function(
	jcross, 
	train, 
	test, 
	n, 
	nq, 
	nv, 
	ng, 
	indclass, 
	data.cat, 
	data.cont, 
	prop.kmeans, 
	means.kmeans, 
	var.kmeans){

# number of categories
numcat = vector(length = nq)
for(i in 1:nq){
	numcat[i] = nlevels(as.factor(data.cat[,i]))
}
nc = max(numcat)   
  

# -- read proportions and means from kmeans output and variance
prop = prop.kmeans
mean.init = t(means.kmeans)
var.init = var.kmeans

# --- other initializations
lambda = array(0, dim = c(nq, max(numcat), ng))   # !!! I change the dim here
mean = array(0, dim = c(nv, ng))     # I change the dim here too
inv.mat = array(0, dim = c(nv, nv, ng)) 
det.val = vector(length = ng, mode = 'numeric')
mat.gum = matrix(nrow = n, ncol =ng)


# call main program
main.res = main.indep(n, nv, nq, ng,nc, indclass, prop, numcat, var.init, inv.mat, det.val, prop.kmeans, means.kmeans, data.cat, data.cont, train, test, lambda, mat.gum)


return(invisible(list(
prop = main.res$prop,
w = main.res$w,
loglik = main.res$vect.loglik,
mat.gum = main.res$mat.gum
)))

}



# ================= gating function and other functions for independence model =================

#------------------------norm function --------------------------------------------
norm.indep = function(x){
	 vect.norm = x/drop(sqrt(crossprod(x)))
	return(vect.norm)
}

#------------------------init function --------------------------------------------

initmain.indep = function(nv, nq, ng, numcat, var.init, inv.mat, det.val, means.kmeans, lambda){

#----fill lambda
for(i in 1:nq){
	for(k in 1:ng){
		lambda[i,1:numcat[i],k] = rep(1/numcat[i], numcat[i])
	}
}

# -- variance and means
var.mat = var.init
mean = means.kmeans


#--inverse matrix and regularize if needed
if (det(var.mat[,,1]) <10^(-10)) {var.mat[,,1] = var.mat[,,1] + diag(0.2, nv, nv);inv.mat[,,1] = solve(var.mat[,,1])} else {inv.mat[,,1] = solve(var.mat[,,1])}
if (det(var.mat[,,2]) <10^(-10)) {var.mat[,,2] = var.mat[,,2] + diag(0.2, nv, nv);inv.mat[,,2] = solve(var.mat[,,2])}  else {inv.mat[,,2] = solve(var.mat[,,2])}

det.val[1] = det(solve(inv.mat[,,1]))
det.val[2] = det(solve(inv.mat[,,2]))

#-- initialize weights
w = matrix(nrow = ng, ncol = nv+nq)
for(k in 1:ng){
	w[k,] = 2 * runif(nv+nq) -1
}
w = t(apply(w,1, norm.indep))   #   w is normed

return(invisible(list(
lambda = lambda,
var.mat = var.mat,
mean = mean,
inv.mat = inv.mat,
det.val = det.val,
w=w
)))
}


#============ gating function ==========================================================

main.indep = function(n, nv, nq, ng,nc, indclass, prop, numcat, var.init, inv.mat, det.val, prop.kmeans, means.kmeans, data.cat, data.cont, train, test, lambda, mat.gum){

##dyn.load('integrativeME.so')

vect.loglik = vector(length = 10, mode = 'numeric')  #  number of iount


cat('test:', test, '\n')

# --- training set
data.cat.learn = data.cat[-test,]       # remove test data
data.cont.learn = data.cont[-test,]


init = initmain.indep(nv, nq, ng, numcat, var.init, inv.mat, det.val, means.kmeans, lambda)
lambda = init$lambda
var.mat = init$var.mat
mean = init$mean
inv.mat = init$inv.mat
det.val = init$det.val
w=init$w



# -----EM ALGORITHM -----------------------------
iount = 1 
while(iount <=10){                            #lim iount

#--- E-step
sw = vector(mode = 'numeric', length = ng)                # is initialized in E-step
swc = array(0, dim = c(ng, nq, nc))  
swx = array(0, dim = c(ng, nv))  
tau = matrix(data = 0 , nrow = length(train), ncol = ng)  # is initialized in the E-step
loglik = 0

cat = matrix(nrow = length(train), ncol = ng, data =1)
for(jstep in 1:length(train)){
	for(kstep in 1:ng){
	for(nqstep in 1:nq){
		cat[jstep,kstep] = cat[jstep,kstep] * lambda[nqstep,(data.cat.learn[jstep,nqstep]), kstep]  ## +1?
} } }

res.estep = .Fortran('eindep',	
		n = as.integer(length(train)),   # careful here n = n-loo
		nv = as.integer(nv),
		ng = as.integer(ng),
		nc = as.integer(nc),
		nq = as.integer(nq),
		x = as.matrix(data.cont.learn),
		xc = as.matrix(data.cat.learn),
		tau = as.matrix(tau),
		w = as.matrix(w),
		det = as.double(det.val),
		inv = as.array(inv.mat),
		indclass = as.integer(indclass),
		mean = as.array(mean),
		lambda = as.array(lambda),
		t = as.double(prop),
		loglik = as.double(loglik),
		cat = as.array(cat)
	)

tau = res.estep$tau

for(m in 1:length(train)){
	if(length(which(is.na(tau[m,])))==ng){tau[m,] = c(0.5, 0.5)}else{
	if(is.na(tau[m,1])) tau[m,1] = 1- tau[m,2]
	if(is.na(tau[m,1])) tau[m,2] = 1- tau[m,1]
	}
	if(sum(tau[m,])==0){tau[m,] = c(0.5, 0.5)}
}

vect.loglik[iount] = res.estep$loglik

sw = apply(tau, 2, sum, na.rm = TRUE)
for(k in 1:ng){
	swx[k,] = apply(tau[,k] * data.cont.learn, 2,sum, na.rm=TRUE)

	for(j in 1:length(train)){
		for(i.nq in 1:nq){
			swc[k,i.nq, data.cat.learn[j, i.nq]] = swc[k,i.nq, data.cat.learn[j, i.nq]] + tau[j,k]
		}  #end i.nq
	}  #end j
}  #end k

if(sum(is.na(prop)) !=0){
prop = prop.kmeans        
#cat('na in prop, estep', '\n')
next
}

# --- M-step

res.mstep = .Fortran('mindep',
		n = as.integer(length(train)),   # careful here n = n-loo
		nv = as.integer(nv),
		ng = as.integer(ng),
		nc = as.integer(nc),
		nq = as.integer(nq),
		x = as.matrix(data.cont.learn),
		xc = as.matrix(data.cat.learn),
		sw = as.vector(sw),
		swx = as.array(swx),
		swc = as.array(swc),
		tau = as.matrix(tau),
		w = as.matrix(w),
		var = as.array(var.mat),
		mean = as.array(mean),
		lambda = as.array(lambda),
		numcat = as.vector(numcat),
		indclass = as.integer(indclass),
		t = as.double(prop)		
	)

lambda = res.mstep$lambda
var.mat = res.mstep$var
mean = res.mstep$mean
w =  res.mstep$w     
w = t(apply(w,1, norm.indep))   #   w is normed
             
prop = res.mstep$t

if(sum(is.na(var.mat)) !=0){     
init.mstep = initmain.indep(nv, nq, ng, numcat, var.init, inv.mat, det.val, means.kmeans, lambda)
lambda = init.mstep$lambda
var.mat = init.mstep$var.mat
mean = init.mstep$mean
inv.mat = init.mstep$inv.mat
det.val = init.mstep$det.val
prop = prop.kmeans         

w=init.mstep$w
#cat('na in var mat, mstep', '\n')
}

for(k in 1:ng){
	if(any(is.na(var.mat[,,k]))) var.mat[,,k] = var.init[,,k]
}

if(any(is.na(mean))) mean[which(is.na(mean), arr.ind = TRUE)] = init$mean[which(is.na(mean), arr.ind = TRUE)]


#--inverse matrix and regularize if needed
if (det(var.mat[,,1]) <10^(-10)) {var.mat[,,1] = var.mat[,,1] + diag(0.2, nv, nv);inv.mat[,,1] = solve(var.mat[,,1])} else {inv.mat[,,1] = solve(var.mat[,,1])}
if (det(var.mat[,,2]) <10^(-10)) {var.mat[,,2] = var.mat[,,2] + diag(0.2, nv, nv);inv.mat[,,2] = solve(var.mat[,,2])}  else {inv.mat[,,2] = solve(var.mat[,,2])}


det.val[1] = det(solve(inv.mat[,,1]))
det.val[2] = det(solve(inv.mat[,,2]))

if(any(is.na(det.val))) det.val = init$det.val

iount = iount+1
} # -- fin iount



#--subroutine check

#---- check for each sample in the test sample
indivtest=0
for(jtest in test){

indivtest = indivtest+1
gum = vector(mode = 'numeric',length=ng)


cat = matrix(nrow = n, ncol = ng, data =1)
for(jstep in 1:n){
	for(kstep in 1:ng){
	for(nqstep in 1:nq){
		cat[jstep,kstep] = cat[jstep,kstep] * lambda[nqstep,(data.cat[jstep,nqstep]), kstep]  ## +1 ??
} } }


res.check = .Fortran('checkindep',
		n = as.integer(n),   # careful here n = n total = num
		nv = as.integer(nv),
		ng = as.integer(ng),
		nc = as.integer(nc),
		nq = as.integer(nq),
		ind = as.integer(jtest),    #the loo obs
		xcj = as.matrix(data.cat),
		xj = as.matrix(data.cont),
		w = as.matrix(w),
		mean = as.array(mean),
		lambda = as.array(lambda),
		inv = as.array(inv.mat),
		det = as.double(det.val),
		t = as.double(prop),
		gum = as.double(gum),
		cat = as.array(cat)
	)
mat.gum[test[indivtest],] = res.check$gum
}  # end jtest


return(invisible(list(
prop = prop,
w = w,
##predicted = predicted,
vect.loglik = vect.loglik,
mat.gum = mat.gum
)))
}  #end function main



