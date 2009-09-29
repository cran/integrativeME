C Copyright (C) 2009 
C Shu-Kay Ng, School of Medecine, Logan Campus, Griffith University Meadowbrook, Australia
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

cka: removed the sw
c      subroutine estep(n,nv,ng,nc, nq, x, xc,sw, swx,swc,
c     +  tau, w, det, inv, indclass, mean, lambda, t,
c     +  loglik)

      subroutine eindep(n,nv,ng,nc, nq, x, xc,
     +  tau, w, det, inv, indclass, mean, lambda, t,
     +  loglik)


      INTEGER n,nv,ng,nc, indclass, nq
c      double precision sw(ng),gum, loglik,swx(ng,nv)
      double precision gum, loglik
      double precision twopi,al(ng)
      double precision x(n,nv),tau(n,ng)
      double precision u(n,ng), w(ng,nv+nq),sumga(ng)
c      double precision swc(ng,nq,nc)
      double precision cat(ng)
      double precision mean(nv,ng),t(ng),det(ng)
      double precision inv(nv,nv,ng)
      double precision lambda(nq,nc,ng)
      integer xc(n,nq)
      intrinsic log, sqrt, exp



c   obtain the sufficient stat
      twopi=2.*3.141592653
      loglik = 0

c ----------LOOP on the samples without loo, here n = num-1, here we
c work directly on the learning data sets

      do 800 j=1,n
        gum=0.
	do 400 k=1,ng

           sumga(k)=0.
	   u(j,k) = 0.

	  cat(k)=1.
	  do 230 i=1,nq
	    cat(k)=cat(k)*lambda(i,xc(j,i), k)
c	  print *,'lambda(i,xc(j,i), k)',lambda(i,xc(j,i), k)

230       continue

          al(k)=0.
	  do 300 l=1,nv
	    do 250 m=1,nv
	      al(k)=al(k)+(x(j,l)-mean(l, k))*inv(l,m,k)*
     +	      (x(j,m)-mean(m,k))
250         continue
300       continue

c KA : I increased the limit from 175 to 200
	  if (al(k).gt.200.) then
	    al(k)=0.
	    goto 400
	  endif

	  al(k)=exp(-.5*al(k))/(sqrt(det(k))*(twopi)**(nv*.5))

c KA I changed the limit from -30 to -40
350       if (t(k).lt.1.E-40.or.al(k).lt.1.E-40) then
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
	  gum=gum+t(k)*cat(k)*al(k)*u(j,k)
	  
400     continue

	if (gum.eq.0.) then
	  do 420 k=1,ng
	    tau(j,k)=0.
420       continue
c	  print *,'zero mixture density for entity in e-step',j
cc          do 200 k=1,ng
cc	  sw(k)=0.
cc	  do 80 i=1,nq
cc	  do 70 kk=1,nc
cc                swc(k,i,kk)=0.
cc70            continue
cc80        continue
cc200       continue
	  goto 800
        endif

	do 550 k=1,ng
	  tau(j,k)=t(k)*cat(k)*al(k)*u(j,k)/gum
cc	  sw(k)=sw(k)+tau(j,k)
c KA: compute loglikelihood
	  loglik = loglik + tau(j,k)* log(t(k)*al(k)*u(j,k))
cc	  print *, 'loglik', loglik
c KA: added more than 2 cat in swc
cc	  do 430 i=1,nq
cc		swc(k,i,xc(j,i))=swc(k,i,xc(j,i))+tau(j,k)
cc430       continue

cc	  do 500 l=1,nv
cc	    swx(k,l)=swx(k,l)+tau(j,k)*x(j,l)
cc500       continue


550     continue
800   continue
      return
      end

C---------------------- M STEP -----------------------------
cKA for R code 

cKA ngme2: removed lind, swxc
cKA changed: mean has 2D
cka added: swx
c sizehy and length of work have been set to N*(3*N +13)/2 + 10, N=nv+nq


      subroutine mindep(n,nv,ng,nc, nq, x, xc, sw, swx, swc,
     +  tau, w, var, mean, lambda, numcat, indclass, t)
c KA added ncat
      INTEGER ng,nc, ncat, indclass
      double precision sw(ng), sum, w(ng,nv+nq),swx(ng,nv)
      double precision swc(ng,nq,nc)
      double precision x(n,nv)
      double precision mean(nv,ng),t(ng)
      double precision var(nv,nv,ng)
      integer n,nv,gp
      double precision tau(n,ng)
c KA added numcat
      integer nq, numcat(nq)
      double precision lambda(nq,nc,ng)
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
      do 250 k=1,ng
	do 200 i=1,nq
	   ncat = numcat(i)
	   do 150 kk =1,ncat
		lambda(i,kk,k)=swc(k,i,kk)/sw(k)
150	   continue
200     continue
250   continue


c     Update t(k) and the mean and variance
      do 500 k=1,ng
	t(k)=sw(k)/sum
	do 400 l=1,nv
	  mean(l,k)=swx(k,l)/sw(k)

cka	  do 280 kk=1,nc
cka	    mean(kk,l,k)=swxc(k,l,kk)/swc(k,lind,kk)
cka280       continue
c	print *, 't(k)', t(k)
	  do 300 m=1,l
	    var(l,m,k)=0.
	  do 295 j=1,n 
	    do 290 kk=1,nc
	    var(l,m,k)=var(l,m,k)+(x(j,l)-mean(l,k))
     +	      *(x(j,m)-mean(m,k))*tau(j,k)
290         continue
295	continue
300       continue
400     continue
     	do 450 l=1,nv
	do 460 m=1,l
		var(l,m,k) = var(l,m,k)/sw(k)
		var(m,l,k)=var(l,m,k)
460     continue
450   	continue
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
c KA put a warning here, sizehy > N*(3*N +13)/2, N=nv+nq, and work(sizehy)
cka	sizehy=2000
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
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N)
C         ----------
C         CALCULATE THE FUNCTIONS AT X AND
C         RETURN THIS VECTOR IN FVEC.
C         ---------
C         RETURN
C         END

C I changed the names, ne to nv, fvecw to fvec
C removed the parameters
C added n, x, tau and indclass

c evalfw is changed here to deal with nv+nq
c added nvx, nqx and xc

      subroutine evalfw(nv, xw, fvec, iflag,
     + n, x, xc, taugp, indclass, nvx, nqx)
      INTEGER nv,iflag, indclass, n, nvx, xc(n,nqx)
c      double precision x(n,nv),taugp(n)
      double precision xw(nv),fvec(nv),tmp(nvx+nqx),sumxi,temp
      double precision x(n,nvx),taugp(n)

      intrinsic exp

      do 200 l=1,nvx+nqx
	tmp(l)=xw(l)
	fvec(l)=0.
200   continue
      do 500 j=1,n
        sumxi=0.
	do 300 l=1,nqx
          sumxi=sumxi+tmp(l)*xc(j,l)
300     continue
	do 400 l=1,nvx
          sumxi=sumxi+tmp(nqx+l)*x(j,l)
400     continue
	temp=exp(sumxi)
	sumxi=temp/(1.+temp)
	temp=0.
	if (j.gt.indclass) temp=1.
	do 420 l=1,nqx
cc          fvecw(l)=fvecw(l)+tau(j,gp)*(temp-sumxi)*xc(j,l)
          fvec(l)=fvec(l)+taugp(j)*(temp-sumxi)*xc(j,l)
420     continue
	do 450 l=1,nvx
cc          fvecw(nq+l)=fvecw(nq+l)+tau(j,gp)*(temp-sumxi)*x(j,l)
          fvec(nqx+l)=fvec(nqx+l)+taugp(j)*(temp-sumxi)*x(j,l)
450     continue
500   continue
      return
      end


c ---------------------- SUBROUTINE CHECK -----------------------     


cKA ngme2: removed lind, indclass
cKA changed: mean has 2D
cka removed: pred, error
 

cc      subroutine check(n,nv,ng,nc,ind,w)
      subroutine checkindep(n,nv,ng,nc, nq, ind,xcj,xj,
     +  w, mean,lambda,inv,det,t, gum)
      INTEGER n,nv,ng,ind,out, nc, nq
      double precision gum(2),twopi,al(ng)
      double precision u(n,ng),w(ng,nv+nq),sumga(ng)
      double precision cat(ng)
      double precision mean(nv, ng),t(ng),det(ng)
      double precision inv(nv,nv, ng)
      double precision lambda(nq,nc, ng)
      integer xcj(n,nq)
      double precision xj(n,nv)

      intrinsic log,sqrt

      twopi=2.*3.141592653
      j=ind
c here deals only with 2 classes
      do 800 jj=1,2 
        out=jj-1
        gum(jj)=0.

	do 400 k=1,ng
          cat(k)=1.
	  do 200 i=1,nq
	     cat(k)=cat(k)*lambda(i,xcj(j,i),k)
200       continue
cc	print *, 'cat', cat(k)
          al(k)=0.
	  do 300 l=1,nv
c ka : should be xj instead of x ?!
c ka and I changed idt(j) to xcj(j,lind)
	    do 250 m=1,nv
	      al(k)=al(k)+(xj(j,l)-mean(l,k))*inv(l,m,k)*
     +	      (xj(j,m)-mean(m,k))
250         continue
300       continue

c KA : changed here from 175 to 200
	  if (al(k).gt.200.) then
c	    write(30,*) 'al(k) gt 175, k', al(k), k
c	    print *, 'al(k) gt 200, k in check', al(k), k
	    al(k)=0.
	    goto 400
	  endif

	  al(k)=exp(-.5*al(k))/(sqrt(det(k))*(twopi)**(nv*.5))
350       if (t(k).lt.1.E-30.or.al(k).lt.1.E-30) then
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
cc	print *, 'sumga', sumga(k)

	  u(j,k)=1./(1.+exp(sumga(k)))
	  if (out.eq.1) u(j,k)=u(j,k)*exp(sumga(k))
	  gum(jj)=gum(jj)+t(k)*cat(k)*al(k)*u(j,k)
400     continue
cc	print *, 'gum', gum(jj)

	if (gum(jj).eq.0.) then
c	  print *,'zero mixture density for entity ',j
	  goto 800
        endif
800   continue

cka      pred=0
cka      if (gum(2).gt. gum(1)) pred=1
cka      error=0
c KA: indicate the class
cka      if (j .le. indclass .and. pred .eq. 1) error=1
cka      if (j .gt. indclass .and. pred .eq. 0) error=1
c      print *, j,gum(1),gum(2),pred,error
cc      print *, 'gum', gum(1),gum(2)
cc      total=total+error

cc      print *, 'total', total
      return 
      end

