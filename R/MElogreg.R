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



MElogreg = function(
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
	prop.kmeans){


# -- read proportions from kmeans
prop = prop.kmeans


# ---------------- other initializations -------------------

mat.softmax.pred = matrix(nrow = length(test), ncol = ng)


main.res = main.logreg(n, nv, nq, ng,indclass, prop, prop.kmeans, data.cat, data.cont, train, test)
##?#mat.prop[jcross,] = main.res$prop
##?#mat.u[test,] = main.res$mat.u[test,] 

return(invisible(list(
prop = main.res$prop,
w = main.res$w,
loglik = main.res$vect.loglik,
mat.gum = main.res$mat.gum,
softmax = main.res$mat.softmax.pred
)))


}


# ================= gating function and other functions for independence model =================

#------------------------norm function --------------------------------------------
norm.logreg = function(x){
	 vect.norm = x/drop(sqrt(crossprod(x)))
	return(vect.norm)
}


#------------------------init function --------------------------------------------

initmain.logreg = function(nv, nq, ng){


#---- initialize weights
w = matrix(nrow = ng, ncol = nv+nq)
for(k in 1:ng){
	w[k,] = 2 * runif(nv+nq) -1
}

w = t(apply(w,1, norm.logreg))   #   w is normed
return(invisible(list(
w=w
)))

}


#============ gating function ==========================================================


main.logreg = function(n, nv, nq, ng,indclass, prop, prop.kmeans, data.cat, data.cont, train, test){

vect.loglik = vector(length = 20, mode = 'numeric')  #  number of iount
mat.gum = matrix(data = 0, nrow = n , ncol = 2)  # 
mat.u = matrix(data = 0, nrow = n , ncol = 2)  # 


# -- define training and test data  ---
                      
cat('test:', test, '\n')

# --- training set
data.cat.learn = data.cat[-test,]       # remove test data
data.cont.learn = data.cont[-test,]

#--- initialize the weights
w = initmain.logreg(nv, nq, ng)$w


#-----------logistic reg on training data-------------------------

mat.softmax = matrix(nrow = length(train), ncol = 2)

glm = glm(type[train]~as.matrix(cbind(data.cat.learn, data.cont.learn)), family=binomial)
mat.softmax = cbind((1-glm$fitted.values), glm$fitted.values)


# ============== EM ALGORITHM ====================
iount = 1 
while(iount <=20){                            #lim iount

#--- E-step

sw = vector(mode = 'numeric', length = ng)        # is initialized in E-step
swx = array(0, dim = c(ng, nv))  
tau = matrix(data = 0 , nrow = length(train), ncol = ng)  # is initialized in the E-step
loglik = 0

if(any(is.na(prop))) prop = prop.kmeans


res.estep = .Fortran('elogreg',	
		n = as.integer(length(train)),   # careful here n = n-loo
		nv = as.integer(nv),
		ng = as.integer(ng),
		nq = as.integer(nq),
		x = as.matrix(data.cont.learn),
		xc = as.matrix(data.cat.learn),
		sw = as.vector(sw),
		swx = as.array(swx),
		tau = as.matrix(tau),
		w = as.matrix(w),
		indclass = as.integer(indclass),
		t = as.double(prop),
		loglik = as.double(loglik),
		softmax = as.array(mat.softmax)
	)

tau = res.estep$tau

for(m in 1:length(train)){
	if(length(which(is.na(tau[m,])))==ng){tau[m,] = c(0.5, 0.5)}else{
	if(is.na(tau[m,1])) tau[m,1] = 1- tau[m,2]
	if(is.na(tau[m,1])) tau[m,2] = 1- tau[m,1]
	}
	if(sum(tau[m,])==0){tau[m,] = c(0.5, 0.5)}
}

#sw = res.estep$sw
#swx = res.estep$swx
vect.loglik[iount] = res.estep$loglik

sw = apply(tau, 2, sum, na.rm = TRUE)
for(k in 1:ng){
swx[k,] = apply(tau[,k] * data.cont.learn, 2,sum, na.rm=TRUE)
}  #end k


# --- M-step

res.mstep = .Fortran('mlogreg',
		n = as.integer(length(train)),   # careful here n = n-loo
		nv = as.integer(nv),
		ng = as.integer(ng),
		nq = as.integer(nq),
		x = as.matrix(data.cont.learn),
		xc = as.matrix(data.cat.learn),
		sw = as.vector(sw),
		swx = as.array(swx),
		tau = as.matrix(tau),
		w = as.matrix(w),
		indclass = as.integer(indclass),
		t = as.double(prop)		
	)
# extract weight and proportion
w =  res.mstep$w 
w = t(apply(w,1, norm.logreg))   #   w is normed
                 
prop = res.mstep$t

iount = iount+1
} # -- fin iount



#-----------subroutine check

# ----- computes softmax for test data---------------
data.test = cbind(data.cat[test,], data.cont[test,])

softmax.pred = 1/(1 + exp(-(as.matrix(data.test)%*%glm$coefficients[-1] + glm$coefficients[1])))

mat.softmax.pred = cbind((1-softmax.pred), softmax.pred)

#---- check for each sample in the test sample
indivtest=0
for(jtest in test){

indivtest = indivtest+1


soft.pred = mat.softmax.pred[indivtest,]
u = vector(length = ng, mode = 'numeric')
pred = 0
gum = vector(mode = 'numeric',length=ng)
res.check = .Fortran('checklogreg',
		n = as.integer(n),   # careful here n = n total = num
		nv = as.integer(nv),
		ng = as.integer(ng),
		nq = as.integer(nq),
		ind = as.integer(jtest),    #the loo obs
		xcj = as.matrix(data.cat),
		xj = as.matrix(data.cont),
		w = as.matrix(w),
		t = as.double(prop),
		u = as.double(u),
		gum = as.double(gum),
		softmax = as.double(soft.pred)
	)
mat.gum[test[indivtest],] = res.check$gum
mat.u[test[indivtest],] = res.check$u
}  # end jtest



return(invisible(list(
prop = prop,
w = w,
vect.loglik = vect.loglik,
mat.gum = mat.gum,
mat.softmax.pred = mat.softmax.pred,
mat.u = mat.u
)))
}  #end function main

