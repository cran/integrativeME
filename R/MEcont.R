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

MEcont = function(
	jcross, 
	train, 
	test, 
	n, 
	nv, 
	ng, 
	indclass, 
	data.cont, 
	prop.kmeans, 
	means.kmeans, 
	var.kmeans){

# -- read proportions and means from kmeans output and variance
prop = prop.kmeans
mean.init = t(means.kmeans)
var.init = var.kmeans

# --- other initializations
mean = array(0, dim = c(nv, ng))     # I change the dim here too
inv.mat = array(0, dim = c(nv, nv, ng)) 
det.val = vector(length = ng, mode = 'numeric')
mat.gum = matrix(nrow = n, ncol =ng)


# call main program
main.res = main.cont(n, nv, ng, indclass, prop, var.init, inv.mat, det.val, prop.kmeans, means.kmeans, data.cont, train, test, mat.gum)


return(invisible(list(
prop = main.res$prop,
w = main.res$w,
loglik = main.res$vect.loglik,
mat.gum = main.res$mat.gum
)))

}


# ================= gating function and other functions for cont model =================

#------------------------norm function --------------------------------------------
norm.cont = function(x){
	 vect.norm = x/drop(sqrt(crossprod(x)))
	return(vect.norm)
}

#------------------------init function --------------------------------------------

initmain.cont = function(nv, ng, var.init, inv.mat, det.val, means.kmeans){

# -- variance and means
var.mat = var.init
mean = means.kmeans


#--inverse matrix and regularize if needed
if (det(var.mat[,,1]) <10^(-10)) {var.mat[,,1] = var.mat[,,1] + diag(0.2, nv, nv);inv.mat[,,1] = solve(var.mat[,,1])} else {inv.mat[,,1] = solve(var.mat[,,1])}
if (det(var.mat[,,2]) <10^(-10)) {var.mat[,,2] = var.mat[,,2] + diag(0.2, nv, nv);inv.mat[,,2] = solve(var.mat[,,2])}  else {inv.mat[,,2] = solve(var.mat[,,2])}

det.val[1] = det(solve(inv.mat[,,1]))
det.val[2] = det(solve(inv.mat[,,2]))

#-- initialize weights
w = matrix(nrow = ng, ncol = nv)
for(k in 1:ng){
	w[k,] = 2 * runif(nv) -1
}
w = t(apply(w,1, norm.cont))   #   w is normed

return(invisible(list(
var.mat = var.mat,
mean = mean,
inv.mat = inv.mat,
det.val = det.val,
w=w
)))
}


#============ gating function ==========================================================

main.cont = function(n, nv, ng, indclass, prop, var.init, inv.mat, det.val, prop.kmeans, means.kmeans, data.cont, train, test, mat.gum){

##dyn.load('integrativeME.so')

vect.loglik = vector(length = 10, mode = 'numeric')  #  number of iount


cat('test:', test, '\n')

# --- training set
data.cont.learn = data.cont[-test,]


init = initmain.cont(nv, ng, var.init, inv.mat, det.val, means.kmeans)
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
swx = array(0, dim = c(ng, nv))  
tau = matrix(data = 0 , nrow = length(train), ncol = ng)  # is initialized in the E-step
loglik = 0

#econt(n,nv,ng, x,
#     +  tau, w, det, inv, indclass, mean, t, loglik)

res.estep = .Fortran('econt',	
		n = as.integer(length(train)),   # careful here n = n-loo
		nv = as.integer(nv),
		ng = as.integer(ng),
		x = as.matrix(data.cont.learn),
		tau = as.matrix(tau),
		w = as.matrix(w),
		det = as.double(det.val),
		inv = as.array(inv.mat),
		indclass = as.integer(indclass),
		mean = as.array(mean),
		t = as.double(prop),
		loglik = as.double(loglik)
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
}  #end k

if(sum(is.na(prop)) !=0){
prop = prop.kmeans        
#cat('na in prop, estep', '\n')
next
}

# --- M-step

#      subroutine mcont(n,nv,ng, x, sw, swx, 
#     +  tau, w, var, mean, indclass, t)

res.mstep = .Fortran('mcont',
		n = as.integer(length(train)),   # careful here n = n-loo
		nv = as.integer(nv),
		ng = as.integer(ng),
		x = as.matrix(data.cont.learn),
		sw = as.vector(sw),
		swx = as.array(swx),
		tau = as.matrix(tau),
		w = as.matrix(w),
		var = as.array(var.mat),
		mean = as.array(mean),
		indclass = as.integer(indclass),
		t = as.double(prop)		
	)

var.mat = res.mstep$var
mean = res.mstep$mean
w =  res.mstep$w     
w = t(apply(w,1, norm.cont))   #   w is normed
             
prop = res.mstep$t

if(sum(is.na(var.mat)) !=0){     
init.mstep = initmain.cont(nv, ng, var.init, inv.mat, det.val, means.kmeans)
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


#      subroutine checkcont(n,nv,ng, ind,xj,
#     +  w, mean,inv,det,t,
#     +  gum)


res.check = .Fortran('checkcont',
		n = as.integer(n),   # careful here n = n total = num
		nv = as.integer(nv),
		ng = as.integer(ng),
		ind = as.integer(jtest),    #the loo obs
		xj = as.matrix(data.cont),
		w = as.matrix(w),
		mean = as.array(mean),
		inv = as.array(inv.mat),
		det = as.double(det.val),
		t = as.double(prop),
		gum = as.double(gum)
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














