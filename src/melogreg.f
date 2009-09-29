C Copyright (C) 2009 
C Kim-Anh Lê Cao, ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
C and Queensland Facility for Advanced Bioinformatics, The University of Queensland, Australia
C
C This program is free software; you can redistribute it and/or
C modify it under the terms of the GNU General Public License
C as published by the Free Software Foundation; either version 2
C of the License, or (at your option) any later version.
C
C This program is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this program; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.




C fortran functions to be called by R
C independence model + binary classification

c ------ SUBROUTINE E STEP ---------------------------------------
cKA: for R code removed iount, xlogl, u
cka removed PARAMETERS and common blocks
cKA added as arguments: det inv indclass mean nq lambda idt t xc

cKA ngme2: removed idt, swxc
cKA changed: mean has 2D
cKA added: swx
ckA       subroutine estep(n,nv,ng,iount,x,sw,swx,swc,tau,u,w)

cka to implement log reg computed in R, 
cka removed lambda, cat, al, swc, inv, det, mean, nc
cka added softmax


cka      subroutine estep(n,nv,ng,nc, nq, x, xc,sw, swx,swc,
cka     +  tau, w, det, inv, indclass, mean, lambda, t,
cka     +  loglik)

      subroutine elogreg(n,nv,ng,nq, x, xc,sw, swx,
     +  tau, w, indclass, t, loglik, softmax)



      INTEGER n,nv,ng, indclass, nq,xc(n,nq)
      double precision sw(ng),gum, loglik,swx(ng,nv)
      double precision twopi
      double precision x(n,nv),tau(n,ng)
      double precision u(n,ng), w(ng,nv+nq),sumga(ng)
cka      double precision swc(ng,nv,nc)
cka      double precision cat(ng)
cka      double precision mean(nv,ng)
      double precision t(ng)
      double precision inv(nv,nv,ng)
cka      double precision lambda(nq,nc,ng)
      integer idt(n)
cka added softmax
      double precision softmax(n,ng)

      intrinsic log, sqrt, exp



c   obtain the sufficient stat
      twopi=2.*3.141592653
      loglik = 0

c ----------LOOP on the samples without loo, here n = num-1, here we
c work directly on the learning data sets

      do 800 j=1,n
	write(30,*) 'indiv', j

        gum=0.
	do 400 k=1,ng

           sumga(k)=0.
	   u(j,k) = 0.

cka	  cat(k)=1.
cka	  do 230 i=1,nq
cka	    cat(k)=cat(k)*lambda(i,xc(j,i), k)
cka230       continue
cc	print *, 'cat', cat(k)

cka          al(k)=0.
cka	  do 300 l=1,nv
cka	    do 250 m=1,nv
cka	      al(k)=al(k)+(x(j,l)-mean(l, k))*inv(l,m,k)*
cka     +	      (x(j,m)-mean(m,k))
cka250         continue
cka300       continue

c KA : I increased the limit from 175 to 250
cka	  if (al(k).gt.200.) then
c	    write(30,*) 'al(k) gt 175, k', al(k), k
c	    print *, 'al(k) gt 200, k', al(k), k
cka	    al(k)=0.
cka	    goto 400
cka	  endif

cka	  al(k)=exp(-.5*al(k))/(sqrt(det(k))*(twopi)**(nv*.5))
cka	  write(30,*) 'al(k)', al(k)

c KA I changed the limit from -30 to -40
cka350       if (t(k).lt.1.E-40.or.al(k).lt.1.E-40) then
350       if (t(k).lt.1.E-40) then
c		write(30,*) 't(k).lt.1.E-40.or.al(k).lt.1.E-40', t(k), al(k)
c		print *, 't(k).lt.1.E-40.or.al(k).lt.1.E-40', t(k), al(k)
		goto 400
	  endif

cc	  sumga(k)=0.  
	  do 355 l=1,nq
	    sumga(k)=sumga(k)+w(k,l)*xc(j,l)
355       continue
	  do 360 l=1,nv
	    sumga(k)=sumga(k)+w(k,nq+l)*x(j,l)
360       continue
	  u(j,k)=1./(1.+exp(sumga(k)))

c KA for those who develop metastases, 
	  if (j.gt.indclass) u(j,k)=u(j,k)*exp(sumga(k))
cka	  gum=gum+t(k)*cat(k)*al(k)*u(j,k)
	  gum=gum+t(k)*softmax(j,k)*u(j,k)
	  
c	  write(30,*) 'sumga(k)', sumga(k)
c	  write(30,*) 'u(j,k)', u(j,k)
400     continue
c	write(30,*) 'gum', gum

	if (gum.eq.0.) then
	  do 420 k=1,ng
	    tau(j,k)=0.
420       continue
c	  print *,'zero mixture density for entity in e-step',j
          do 200 k=1,ng
	  sw(k)=0.
cka	  do 80 i=1,nq
cka	  do 70 kk=1,nc
cka                swc(k,i,kk)=0.
cka70            continue
cka80        continue
200       continue
	  goto 800
        endif

	do 550 k=1,ng
cka	  tau(j,k)=t(k)*cat(k)*al(k)*u(j,k)/gum
	  tau(j,k)=t(k)*softmax(j,k)*u(j,k)/gum

	  sw(k)=sw(k)+tau(j,k)
c KA: compute loglikelihood
cka	  loglik = loglik + tau(j,k)* log(t(k)*al(k)*u(j,k))
	  loglik = loglik + tau(j,k)* log(softmax(j,k)*u(j,k))

cc	  print *, 'loglik', loglik
c KA: added more than 2 cat in swc
cka	  do 430 i=1,nq
cka		swc(k,i,xc(j,i))=swc(k,i,xc(j,i))+tau(j,k)
cka430       continue

	  do 500 l=1,nv
	    swx(k,l)=swx(k,l)+tau(j,k)*x(j,l)
500       continue


550     continue
800   continue
      return
      end

C---------------------- M STEP -----------------------------
cKA for R code 

cKA ngme2: removed lind, swxc
cKA changed: mean has 2D
cka added: swx

cka to implement log reg computed in R, 
cka removed lambda, numcat, ncat, swc, var, mean, nc

c sizehy and length of work have been set to N*(3*N +13)/2 + 10, N=nv+nq


cka      subroutine mstep(n,nv,ng,nc, nq, x, xc, sw, swx, swc,
cka     +  tau, w, var, mean, lambda, numcat, indclass, t)
      subroutine mlogreg(n,nv,ng,nq, x, xc, sw, swx,
     +  tau, w, indclass, t)


c KA added ncat
      INTEGER ng, indclass
      double precision sw(ng), sum, w(ng,nv+nq),swx(ng,nv)
c      double precision swc(ng,nv,nc)
      double precision x(n,nv)
cka      double precision mean(nv,ng)
      double precision t(ng)
cka      double precision var(nv,nv,ng)
      integer n,nv,gp
      double precision tau(n,ng)
cka      integer nq, numcat(nv)
cka      double precision lambda(nv,nc,ng)
c external function:
      external evalfw, HYBRD1
      integer sizehy,ifail, iflag,xc(n,nq)
c KA : be careful to allocate right dimension for work()
      double precision tem(nv+nq),xtol, taugp(n)
      double precision fvec(nv+nq)
      intrinsic exp
      double precision work((nv+nq)*(3*(nv+nq) +13)/2 + 10)



      ifail = 0

      sum=0.
      do 100 k=1,ng
	sum=sum+sw(k)
100   continue
c      print *, 'sum', sum

c  Update lambda  KA: added for more than 2 cat
cka      do 250 k=1,ng
cka	do 200 i=1,nq
cka	   ncat = numcat(i)
cka	   do 150 kk =1,ncat
cka		lambda(i,kk,k)=swc(k,i,kk)/sw(k)
cka150	   continue
cka200     continue
cka250   continue


c     Update t(k) and the mean and variance
      do 500 k=1,ng
	t(k)=sw(k)/sum
cka	do 400 l=1,nv
cka	  mean(l,k)=swx(k,l)/sw(k)

cka	  do 280 kk=1,nc
cka	    mean(kk,l,k)=swxc(k,l,kk)/swc(k,lind,kk)
cka280       continue
c	print *, 't(k)', t(k)
cka	  do 300 m=1,l
cka	    var(l,m,k)=0.
cka	  do 295 j=1,n 
cka	    do 290 kk=1,nc
cka	    var(l,m,k)=var(l,m,k)+(x(j,l)-mean(l,k))
cka     +	      *(x(j,m)-mean(m,k))*tau(j,k)
cka290         continue
cka295	continue
cka300       continue
cka400     continue
cka     	do 450 l=1,nv
cka	do 460 m=1,l
cka		var(l,m,k) = var(l,m,k)/sw(k)
cka		var(m,l,k)=var(l,m,k)
cla460     continue
cka450   	continue
500   continue
cc      call matinv(nv,ng,var,inv,det)   to inverse later !!
c  Update w

      do 600 gp=1,ng
	do 510 l=1,nv+nq
	  tem(l)=w(gp,l)
510     continue

	do 520 j=1,n
	  taugp(j)=tau(j,gp)
520     continue

c ---- call HYBRD1
	xtol=0.000000001
c KA put a warning here, sizehy > N*(3*N +13)/2, N=nv+nq
	sizehy=(nv+nq)*(3*(nv+nq) +13)/2 + 10
c	call HYBRD1(evalfw,nv+nq,tem,fvec,xtol,ifail,work,sizehy)

	call HYBRD1(nv+nq,tem,fvec,xtol,ifail,work,sizehy,
     +  n, x, xc, taugp, indclass, nv, nq)

c	if (ifail.ne.1) print *, 'ifail= ',ifail,gp
c	print *, (fvec(l),l=1,nv+nq)
	do 530 l=1,nv+nq
	  w(gp,l)=tem(l)
530     continue
600     continue
c       do 700 k=1,ng
c	write (30,*) 'swc: ', (swc(k,l,1),l=1,nq)
c	write (30,*) 'lambda: ', (lambda(k,l,1),l=1,nq)
c700   continue

      return 
      end


C ------------------------ EVALFW --------------------------------
C has already been defined somewhere else

c ---------------------- SUBROUTINE CHECK -----------------------     


cKA ngme2: removed lind, indclass
cKA changed: mean has 2D
cka added:

cka to implement log reg computed in R, 
cka removed lambda, numcat, ncat, swc, mean, inv, det, nc
cka added ntest

cka4: removed pred and added the output u
 

cc      subroutine check(n,nv,ng,nc,ind,w)
cka      subroutine check(n,nv,ng,nc, nq, ind,xcj,xj,
cka     +  w, mean,lambda,inv,det,t, pred,
cka     +  gum)

cka4      subroutine check(n,nv,ng, nq, ind,xcj,xj,
cka4     +  w, t, pred, gum, softmax)

      subroutine checklogreg(n,nv,ng, nq, ind,xcj,xj,
     +  w, t, u, gum, softmax)



cka4      INTEGER n,nv,ng,ind,out,error, pred
      INTEGER n,nv,ng,ind,out,error
      double precision gum(ng),twopi,al(ng)
cka4      double precision u(n,ng),w(ng,nv+nq),sumga(ng)
      double precision u(ng),w(ng,nv+nq),sumga(ng)

cka      double precision cat(ng)
cka      double precision mean(nv, ng),det(ng)
      double precision t(ng)
cka      double precision inv(nv,nv, ng)
      integer nq
cka      double precision lambda(nv,nc, ng)
      integer xcj(n,nv)
      double precision xj(n,nv)
cka added softmax
      double precision softmax(ng)

      intrinsic log,sqrt

      twopi=2.*3.141592653
      j=ind
c here deals only with 2 classes
      do 800 jj=1,2 
        out=jj-1
        gum(jj)=0.

	do 400 k=1,ng
cka          cat(k)=1.
cka	  do 200 i=1,nq
cka	     cat(k)=cat(k)*lambda(i,xcj(j,i),k)
cka200       continue
cc	print *, 'cat', cat(k)
cka          al(k)=0.
cka	  do 300 l=1,nv
c ka : should be xj instead of x ?!
c ka and I changed idt(j) to xcj(j,lind)
cka	    do 250 m=1,nv
cka	      al(k)=al(k)+(xj(j,l)-mean(l,k))*inv(l,m,k)*
cka     +	      (xj(j,m)-mean(m,k))
cka250         continue
cka300       continue

c KA : changed here from 175 to 200
cka	  if (al(k).gt.200.) then
c	    write(30,*) 'al(k) gt 175, k', al(k), k
c	    print *, 'al(k) gt 200, k in check', al(k), k
cka	    al(k)=0.
cka	    goto 400
cka	  endif

cka	  al(k)=exp(-.5*al(k))/(sqrt(det(k))*(twopi)**(nv*.5))
cka350       if (t(k).lt.1.E-30.or.al(k).lt.1.E-30) then
350       if (t(k).lt.1.E-30) then
c		write(30,*) 't(k).lt.1.E-40.or.al(k).lt.1.E-40', t(k), al(k)
c	        print *, 't(k)ltE-40.or.al(k)ltE-40 in check',t(k),al(k)		
		goto 400
	  endif
cc	print *, 'al(k)', al(k)
c ka: should be xcj and xj ?!
	  sumga(k)=0.  
	  do 355 l=1,nq
	    sumga(k)=sumga(k)+w(k,l)*xcj(j,l)
355       continue
	  do 360 l=1,nv
	    sumga(k)=sumga(k)+w(k,nq+l)*xj(j,l)
360       continue
c	print *, 'sumga', sumga(k)

cka4	  u(j,k)=1./(1.+exp(sumga(k)))
cka4	  if (out.eq.1) u(j,k)=u(j,k)*exp(sumga(k))
	  u(k)=1./(1.+exp(sumga(k)))
	  if (out.eq.1) u(k)=u(k)*exp(sumga(k))
cka	  gum(jj)=gum(jj)+t(k)*cat(k)*al(k)*u(j,k)
cka4	  gum(jj)=gum(jj)+t(k)*softmax(k)*u(j,k)
	  gum(jj)=gum(jj)+t(k)*softmax(k)*u(k)

400     continue
cc	print *, 'gum', gum(jj)

	if (gum(jj).eq.0.) then
c	  print *,'zero mixture density for entity ',j
	  goto 800
        endif
800   continue

cka4      pred=0
cka3      if (gum(2).gt. gum(1)) pred=1
cka4      if (gum(2).gt. 0.5) pred=1
cka4      error=0
c KA: indicate the class
cka      if (j .le. indclass .and. pred .eq. 1) error=1
cka      if (j .gt. indclass .and. pred .eq. 0) error=1
c      print *, j,gum(1),gum(2),pred,error
cc      print *, 'gum', gum(1),gum(2)
cc      total=total+error

cc      print *, 'total', total
      return 
      end

