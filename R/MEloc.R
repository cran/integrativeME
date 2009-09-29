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



MEloc = function(
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
	var.kmeans,
	loc.ind){

# number of categories
numcat = vector(length = nq)
for(i in 1:nq){
	numcat[i] = nlevels(as.factor(data.cat[,i]))
}
nc = nlevels(as.factor(data.cat[,loc.ind]))    #number of cat in the location 

# -- read proportions and means from kmeans output and variance
prop = prop.kmeans
means.init = means.kmeans

var.init = array(0, dim = c(nv, nv, ng)) 

var.init = var.kmeans

# --- other initializations
lambda = array(0, dim = c(nq, max(numcat), ng))   # !!! I change the dim here
var.mat = array(0, dim = c(nv, nv, ng)) 
inv.mat = array(0, dim = c(nv, nv, ng)) 
det.val = vector(length = ng, mode = 'numeric')
mat.gum = matrix(nrow = n, ncol =ng)


main.res = main.loc(n, nv, nq, ng,nc, indclass, prop, numcat, var.init, inv.mat, det.val, prop.kmeans, means.init, data.cat, data.cont, train, test, lambda, mat.gum, loc.ind)

return(invisible(list(
prop = main.res$prop,
w = main.res$w,
loglik = main.res$vect.loglik,
mat.gum = main.res$mat.gum
)))
}


# ================= gating function and other functions for independence model =================

#------------------------norm function --------------------------------------------
norm.loc = function(x){
	 vect.norm = x/drop(sqrt(crossprod(x)))
	return(vect.norm)
}

norm2.loc = function(x){
	vect.norm = x/sum(x)
	return(vect.norm)
}

#------------------------init function --------------------------------------------

initmain.loc = function(nv, nq, ng, nc, numcat, var.init, inv.mat, det.val, means.init, lambda){

#----fill lambda
for(i in 1:nq){
	for(k in 1:ng){
		lambda[i,1:numcat[i],k] = rep(1/numcat[i], numcat[i])
	}
}

# -- variance and means
var.mat = var.init

mean = array(0, dim = c(nc, nv, ng))     # I change the dim here too  mean = array(0, dim = c(max(numcat), nv, ng))  
for( i in 1:nc){
	for(nvj in 1:nv){
	mean[i,nvj,1] =  means.init[1,nvj]   #+ runif(1) - 0.5
	mean[i,nvj,2] =  means.init[2,nvj]   #+ runif(1) - 0.5
}
}

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
w = t(apply(w,1, norm.loc))   #   !!!! norm w


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

main.loc = function(n, nv, nq, ng,nc, indclass, prop, numcat, var.init, inv.mat, det.val, prop.kmeans, means.init, data.cat, data.cont, train, test, lambda, mat.gum, loc.ind){

vect.loglik = vector(length = 10, mode = 'numeric')  #  number of iount
mat.gum = matrix(data = 0, nrow = n , ncol = 2)  # 

cat('test:', test, '\n')

# --- training set
data.cat.learn = as.matrix(data.cat[-test,] )      # remove test data
data.cont.learn = as.matrix(data.cont[-test,])

idt = data.cat.learn[,loc.ind]  #valeur de la location sans loo


init = initmain.loc(nv, nq, ng, nc, numcat, var.init, inv.mat, det.val, means.init, lambda)
lambda = init$lambda
var.mat = init$var.mat
mean = init$mean
inv.mat = init$inv.mat
det.val = init$det.val
w=init$w

# ============== EM ALGORITHM ====================
iount = 1 
while(iount <=10){                            #lim iount

#--- E-step
sw = vector(mode = 'numeric', length = ng)        # is initialized in E-step
swc = array(0, dim = c(ng, nq, max(numcat)))  
swxc = array(0, dim = c(ng, nv, nc))  
tau = matrix(data = 0 , nrow = length(train), ncol = ng)  # is initialized in the E-step
loglik = 0

mat.cat = matrix(data = 1 , nrow = length(train), ncol = ng)

for(j in 1:length(train)){
for(k in 1:ng){
for(iq in 1:nq){
	    mat.cat[j,k]=mat.cat[j,k]*lambda[iq,data.cat.learn[j,iq], k]
}#end iq
}}


#cat('estep', iount)
res.estep = .Fortran('eloc',	
		n = as.integer(length(train)),   ### careful here n = n-loo
		nv = as.integer(nv),
		ng = as.integer(ng),
		nc = as.integer(nc),
		nq = as.integer(nq),
		maxcat = as.integer(max(numcat)),  #
		x = as.matrix(data.cont.learn),
		xc = as.matrix(data.cat.learn),
		tau = as.matrix(tau),
		w = as.matrix(w),
		det = as.double(det.val),
		inv = as.array(inv.mat),
		indclass = as.integer(indclass),
		mean = as.array(mean),
		lambda = as.array(lambda),
		idt = as.integer(idt),
		t = as.double(prop),
		loglik = as.double(loglik),
		cat = as.matrix(mat.cat)
	)

tau = res.estep$tau
for(m in 1:length(train)){
	if(length(which(is.na(tau[m,])))==ng){tau[m,] = c(0.5, 0.5)}else{
	if(is.na(tau[m,1])) tau[m,1] = 1- tau[m,2]
	if(is.na(tau[m,1])) tau[m,2] = 1- tau[m,1]
	}
	if(sum(tau[m,])==0){tau[m,] = c(0.5, 0.5)}
}


sw = apply(tau, 2, sum, na.rm = TRUE)
if(any(sw==0)){sw[which(sw==0)] = 0.01}

for(k in 1:ng){
for(j in 1: length(train)){
for(l in 1:nv){
		swxc[k,l,idt[j]] = sum(swxc[k,l,idt[j]], tau[j,k] * data.cont.learn[j,l], na.rm=TRUE)
} # end l
	for(iq in 1:nq){
		swc[k,iq, data.cat.learn[j,iq]] = sum(swc[k,iq, data.cat.learn[j,iq]], tau[j,k], na.m=TRUE)
	}

}  #end j
} #end k

vect.loglik[iount] = res.estep$loglik

if(sum(is.na(prop)) !=0){
##prop = prop.ngme
prop = prop.kmeans
#next
cat('na in prop, estep', '\n')
next
}

# --- M-step
#      subroutine mstep(n,nv,ng,nc, nq, x,sw,swc,swxc,
#     +  tau, w, var, mean, lambda, loc.ind)

#cat('mstep', iount)
res.mstep = .Fortran('mloc',
		n = as.integer(length(train)),   ### careful here n = n-loo
		nv = as.integer(nv),
		ng = as.integer(ng),
		nc = as.integer(nc),
		nq = as.integer(nq),
		maxcat = as.integer(max(numcat)),
		x = as.matrix(data.cont.learn),
		xc = as.matrix(data.cat.learn),
		sw = as.vector(sw),
		swc = as.array(swc),
		swxc = as.array(swxc),
		tau = as.matrix(tau),
		w = as.matrix(w),
		var = as.array(var.mat),
		mean = as.array(mean),
		lambda = as.array(lambda),
		loc.ind = as.integer(loc.ind),
		numcat = as.vector(numcat),
		indclass = as.integer(indclass),
		t = as.double(prop)		
	)

lambda = res.mstep$lambda
for(k in 1:ng){
if(nq > 1){lambda[,,k] = t(apply(lambda[,,k],1,norm2.loc))} else{lambda[,,k] = norm2.loc(lambda[,,k])} 
}                                            #norm lambda


var.mat = res.mstep$var
mean = res.mstep$mean
w =  res.mstep$w
w = t(apply(w,1, norm.loc))   #   !!!! norm w
                  
prop = res.mstep$t

mean[which(mean=='Inf', arr.ind = TRUE)] = 0
mean[which(mean=='-Inf', arr.ind = TRUE)] = 0


#--inverse matrix and regularize if needed
if (det(var.mat[,,1]) <10^(-10)) {var.mat[,,1] = var.mat[,,1] + diag(0.2, nv, nv);inv.mat[,,1] = solve(var.mat[,,1])} else {inv.mat[,,1] = solve(var.mat[,,1])}
if (det(var.mat[,,2]) <10^(-10)) {var.mat[,,2] = var.mat[,,2] + diag(0.2, nv, nv);inv.mat[,,2] = solve(var.mat[,,2])}  else {inv.mat[,,2] = solve(var.mat[,,2])}

det.val[1] = det(solve(inv.mat[,,1]))
det.val[2] = det(solve(inv.mat[,,2]))

iount = iount+1
} # -- fin iount


#--------------------------------------subroutine check
#      subroutine check(n,nv,ng,nc, nq, ind,xcj,xj,
#     +  w, mean,lambda,inv,det,loc.ind,t,indclass,total)


#---- check for each sample in the test sample
indivtest=0
for(jtest in test){
#cat('jtest, check', jtest)
indivtest = indivtest+1

gum = vector(mode = 'numeric',length=ng)
res.check = .Fortran('checkloc',
		n = as.integer(n),   ### careful here n = n total = num
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
		loc.ind = as.integer(loc.ind),
		t = as.double(prop),
		gum = as.double(gum),
		maxcat = as.integer(max(numcat))
	)

#predicted[jtest] = res.check$pred
#mat.gum[jtest,] = res.check$gum
mat.gum[test[indivtest],] = res.check$gum

}  # end jtest


return(invisible(list(
prop = prop,
w = w,
#predicted = predicted,
vect.loglik = vect.loglik,
mat.gum = mat.gum
)))
}  #end function main



