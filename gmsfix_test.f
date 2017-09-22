
      program mjo
      
c     simple 1D mjo model

      implicit real(a-h,o-z)      
      parameter(ny=1000) !grid points
      parameter(pi=3.1415926)
      

      common /grid/yl,y(ny),yh(ny)
      common /greensfn/g(ny,ny)

c*****h is column water vapor, pt is "temperature"

      dimension h(ny,2),v(ny,2)
      dimension pt(ny,2),vpt(ny,2)
      dimension atendpt(ny),tendpt(ny)
      dimension evap(ny),xm(ny),prec(ny)
      dimension atendh(ny),difftendh(ny)
      dimension tendh(ny)
      dimension hold(ny),ptold(ny)
      dimension u(ny)
      dimension evap0(ny),evapfac(ny),vflux(ny)
      dimension q0(ny) !artificial forcing for kelvin wave

      real kh

      i_emean=0      !flag to decide whether mean evap constrained
      i_pertadv=1    !flag to decide if perturbations self-advect
      i_vflux=1      !flag to decide whether to add wind perturbations to flux
      i_kcpl=0         !flag to decide whether to couple Kelvin to WTG mode
      i_gmsvar=0     !is GMS variable or constant

c      set up grid - at moment distance measured in 
c     deformation radii
      yl=13.4
c      dy=yl/ny
      dy=2.*yl/float(ny)

      print*,'dy = ',dy
      do i=1,ny
         y(i)=-yl+(float(i)-0.5)*dy !this is the grid for h
                                    !u is staggered, but still
                                    !same dimension, since periodic
c        print*,y(i)
      end do

c      do i=1,ny
c         y(i)=(i-0.5)*dy
c      end do
      
      dt=0.001      !time step (days)
      nsteps=16000 !total time steps
      nrobert=10    !after how many time steps apply Robert filter
      afilt=0.2     !Robert filter parameter
      nout=80      !output after how many time steps
      kh=0.0001     !diffusivity
      r=0.1       !cloud-radiative feedback parameter
      r0=4.8         !clear-sky rad cooling, constant
      ec=4.8         !evap if constant
c      wmax=65.      !saturation water vapor path (mm)
      wmax=70.
      wsfc=68.      !sat wvp related to SST, for sfc fluxes
      xlsfc=1.      !length scale (still in l_d) for sfc fluxes
      adc=15.6      !constants from Bretherton et al. (2004)
      rdc=0.603     !relationship between precip and wvp
      xmc=0.1      !normalized gms if constant
      emean=5.2      !prescribed mean evaporation in mm/d!
      vmean=0.3     !prescribed background wind
c      vmean=0.
      ck=1.         !Free Kelvin wave phase speed
      epsk=0.1      !Free Kelvin - Moisture wave coupling parameter
      dampk=0.1     !Kelvin wave damping rate

      adv_fac=1   !multiply advective tendency by this (like V1b1)
      nshift=0      !number of grid points to shift u forward

      evap0b=87.5   !evaporation at zero wind speed in the warm pool
      evap0p=12.5
      evapfacb=6.25   !additional evaporation per m/s in the warm pool
      evapfacp=1.25

c*****GMS parameterization, triangle function
c     
      xmin=0.7 !value at which GMS is minimum
      xm0=0.4  !value of GMS at zero water vapor
      xm1=0.3 !value at saturation
      xm_min=-0.1 !minimum value of GMS
      s1=(xm0-xm_min)/xmin !slope to left of minimum
      s2=(xm1-xm_min)/(1.-xmin) !slope to right of minimum

c      for constant warm pool conditions, use evap0=100 
c     and evapfac=7.5

      do i=1,ny
         evap0(i)=100.
         evapfac(i)=7.5
c      evap0(i)=evap0b + evap0p*sin(2.*pi*y(i)/(2.*yl))
c      evapfac(i)=evapfacb + evapfacp*sin(2.*pi*y(i)/yl)
      end do
c      r=(xmc/(1-xmc))-0.02 !just past critical value
      print*,'r = ',r

c     construct projection function 
c     only need to do it once

      call greens(dy)

      jrec=1
      open(18,file='gfunc.dat',form='unformatted',
     &     access='direct',recl=ny*4)
        do i=1,ny
           write(18,rec=jrec)(g(i,j),
     &          j=1,ny)
           jrec=jrec+1
         end do
         close(18)
c     set v constant for practice problem
       do i=1,ny
          v(i,1)=vmean
          v(i,2)=vmean
          vflux(i)=vmean
          vpt(i,2)=ck
          vpt(i,1)=ck
       end do

c     set initial condition on h
       y1=0.4
       y2=0.5
       y3=0.6
       hmax=1.
       sigy=0.05
       do i=1,ny
c     triangle - very hard to do
c$$$          if(y(i).lt.y1.or.y(i).gt.y3)then
c$$$             h(i,1)=0
c$$$          else if(y(i).gt.y1.and.y(i).lt.y2)then
c$$$             h(i,1)=hmax*(y(i)-y1)/(y2-y1)
c$$$          else 
c$$$             h(i,1)=hmax-hmax*(y(i)-y2)/(y3-y2)
c$$$          end if

c     Gaussian
c          h(i,1)=hmax*exp(-(y(i)-y2)**2/(sigy**2))

c     tanh
c          h(i,1)=0.5*(tanh(10*(y(i)-y2))+1)

c     sine
c          h(i,1)=1.+0.01*sin(4*pi*y(i))

c          h(i,2)=0.

c     constant
          h(i,1)=50. + 4*sin(2*pi*y(i)/(2*yl));
c          H(i,1)=50. 
          h(i,2)=h(i,1)
          u(i)=0.
          pt(i,1)=0
          pt(i,2)=0
          q0(i)=1.*sin(pi*y(i)/yl);

c          print*,h(i,1)
       end do

       ny2=ny/2
c$$$       do i=ny2+1,ny
c$$$          h(i,1)=h(ny+1-i,1)
c$$$       end do

c     open output file
      
      open(12,file='h.dat',form='unformatted',
     &     access='direct',recl=ny*4)
      open(13,file='dt.dat',form='unformatted',
     &     access='direct',recl=ny*4)
      open(14,file='evap.dat',form='unformatted',
     &     access='direct',recl=ny*4)
      open(15,file='prec.dat',form='unformatted',
     &     access='direct',recl=ny*4)
      open(16,file='u.dat',form='unformatted',
     &     access='direct',recl=ny*4)
         open(17,file='pt.dat',form='unformatted',
     &     access='direct',recl=ny*4)
         open(18,file='atendh.dat',form='unformatted',
     &     access='direct',recl=ny*4)
           open(19,file='v.dat',form='unformatted',
     &     access='direct',recl=ny*4)

      irec=1
      write(12,rec=irec)(h(i,1),
     &     i=1,ny)
      write(14,rec=irec)(evap(i),
     &     i=1,ny)
      write(15,rec=irec)(prec(i),
     &     i=1,ny)
      write(16,rec=irec)(u(i),
     &     i=1,ny)
      write(17,rec=irec)(pt(i,1),
     &     i=1,ny)
      write(18,rec=irec)(atendh(i),
     &     i=1,ny)
         write(19,rec=irec)(v(i,1),
     &     i=1,ny)
c     do one Euler step to get going:  call advection

      call advectuh(v,h,atendh,dy,1)
      call advectuh(vpt,pt,atendpt,dy,1)

      do i=1,ny
         tendh(i)=-atendh(i)
      end do
      do i=1,ny
         h(i,2)=h(i,1)+dt*tendh(i)
         pt(i,2)=pt(i,1)-dt*atendpt(i)
      end do

      irec=irec+1 

      time=0.

c*****MASTER TIME LOOP**************************************************
c*****now we are doing leapfrog time steps
      countfilt=0               !counter for the robert filter      
      countout=0                !counter for output


      do it=1,nsteps

         do k=1,2
            if(k.eq.1)then
               m=2
            else
               m=1
            end if

            if(i_pertadv.eq.1)then
               do i=1,ny
                  v(i,k)=vmean+u(i)
               end do
            end if
             if(i_vflux.eq.1)then
               do i=1,ny
                  vflux(i)=vmean+u(i)
               end do
            end if
               
         call advectuh(v,h,atendh,dy,m)
         call advectuh(vpt,pt,atendpt,dy,m)
        
         call diffuh(kh,h,difftendh,dy,k)

         do i=1,ny
c            evap(i)=v(i,k)*(wsfc-h(i,k))/xlsfc  !mm/d?
c            evap(i)=(((u(i)*u(i))+umsq)**(0.5))
c     &           *(wsfc-h(i,k))/xlsfc !mm/d?
            uspd=abs(vflux(i)*17)

c     pure linear regression formula from Maloney simulation
            evap(i)=(evap0(i)+evapfac(i)*uspd)/29
            prec(i)=exp(adc* ( (h(i,k)/wmax)-rdc))
            xm(i)=xmc
         end do

c     "correct" evaporation to insure specific mean value
c     AHS 11/09
         if(i_emean.eq.1)then
            evapmean=0.
            do i=1,ny
               evapmean=evapmean+evap(i)
            end do
         
            evapmean=evapmean/float(ny)
            diff=emean-evapmean

            do i=1,ny
               evap(i)=evap(i)+diff
            end do
         end if

c     compute winds using projection - so far used only
c     for fluxes, not yet advection 

         call project(prec,u,dy)

c     arbitrarily shift winds forward some number of grid points
c     for current parameters (jan. 10) dx=0.02

         call shiftu(u,nshift)

         if(i_gmsvar.eq.1)then
c     compute gross moist stability
            do i=1,ny
               xx=h(i,k)/wmax
c     xm(i)=s*(xx-alpha)*(xx-alpha) + gamma*xx + beta
               if(xx.le.xmin)then
                  xm(i)=xm0-s1*xx
               else
                  xm(i)=xm_min+s2*(xx-xmin)
               end if
            end do
 3       end if
c     compute total tendencies of h and pt

         
         do i=1,ny     
            if(r*prec(i).gt.r0)then
               rm=r0
            else
               rm=r*prec(i)
            end if
            tendh(i)=-adv_fac*atendh(i)
c     &           -(xm(i)*(1+r)-r)*prec(i)
     &           -xm(i)*prec(i)-rm*(xm(i)-1.)
     &           +evap(i)
     &           -(1-xm(i))*r0
     &           +difftendh(i)
            tendpt(i)=-atendpt(i)+epsk*(prec(i)-r0)
     &           -dampk*pt(i,k)           
         end do         

         if(i_kcpl.eq.1)then
            do i=1,ny
               tendh(i)=tendh(i)+epsk*atendpt(i)
            end do
         end if
          
c        print*,tendpt(500)
c     save old values for the R-A filter

         do i=1,ny
            hold(i)=h(i,k)
            ptold(i)=pt(i,k)
         end do

         
         do i=1,ny
            h(i,k)=h(i,k)+2*dt*tendh(i)
c            pt(i,k)=pt(i,k)+2*dt*tendpt(i)
            pt(i,k)=pt(i,k)-2*dt*atendpt(i)
         end do         
         

c     filter

         do i=1,ny
            h(i,m)=h(i,m)+afilt*(hold(i)-2.*h(i,m)+h(i,k))
            pt(i,m)=pt(i,m)+afilt*(ptold(i)-2.*pt(i,m)+pt(i,k))
            
         end do

         countfilt=countfilt+1
         countout=countout+1

         time=time+dt

c     write output
      if(countout.eq.nout)then
         write(12,rec=irec)(h(i,2),
     &        i=1,ny)
         write(13,rec=irec)(difftendh(i),
     &        i=1,ny)
         write(14,rec=irec)(evap(i),
     &     i=1,ny)
         write(15,rec=irec)(prec(i),
     &     i=1,ny)
         write(16,rec=irec)(u(i),
     &     i=1,ny)
         write(17,rec=irec)(pt(i,2),
     &     i=1,ny)
         write(18,rec=irec)(atendh(i),
     &     i=1,ny)
          write(19,rec=irec)(v(i,2),
     &     i=1,ny)
         irec=irec+1
         countout=0
         
      end if

      end do
      
      end do

      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)

      end

  
      subroutine advectuh(v,x,atend,dy,k)
      parameter(ny=1000) !grid points

!     simple 1st order upwind differencing!

      common /grid/yl,y(ny),yh(ny)
      real v(ny,2),x(ny,2),atend(ny)
      real dy


      do i=2,ny-1
      
      if(v(i,k).gt.0)then
         atend(i)=v(i,k)*(x(i,k)-x(i-1,k))/dy
      else if(v(i+1,k).lt.0)then
         atend(i)=v(i+1,k)*(x(i+1,k)-x(i,k))/dy
      else 
         atend(i)=0.
      end if

      end do

c     boundary points - periodic!

      if(v(2,k).lt.0)then
         atend(1)=v(2,k)*(x(2,k)-x(1,k))/dy
      else if(v(1,k).gt.0)then
         atend(1)=v(1,k)*(x(1,k)-x(ny,k))/dy
      else
         atend(1)=0
      end if

      if(v(ny,k).gt.0)then
         atend(ny)=v(ny,k)*(x(ny,k)-x(ny-1,k))/dy
      else if(v(ny,k).lt.0)then
         atend(ny)=v(ny,k)*(x(1,k)-x(ny,k))/dy
      else
         atend(ny)=0
      end if

      end

      subroutine diffuh(kx,x,difftend,dy,k)
      parameter(ny=1000) !grid points

      common /grid/yl,y(ny),yh(ny)
      real x(ny,2),difftend(ny)
      real dy,kx
      integer k,i

      do i=2,ny-1      
         difftend(i)=kx*(x(i+1,k)-2*x(i,k)+x(i-1,k))/(dy*dy)
      end do

      difftend(1)=kx*(x(2,k)-2.*x(1,k)+x(ny,k))/(dy*dy)

      difftend(ny)=kx*(x(1,k)-2.*x(ny,k)+x(ny-1,k))/(dy*dy)
      
      end

      subroutine advectv(v,atend,dy,k)
      parameter(ny=1000) !grid points

!     simple 1st order upwind differencing!

      common /grid/yl,y(ny),yh(ny)
      real v(ny,2),atend(ny)
      real dy
      integer k

      do i=2,ny-1
      
      if(v(i,k)>0)then
         atend(i)=0.5*(v(i,k)+v(i-1,k))*(v(i,k)-v(i-1,k))/dy
      else 
         atend(i)=0.5*(v(i,k)+v(i+1,k))*(v(i+1,k)-v(i,k))/dy
      end if

      end do

c     set boundary points advection to zero, since v=0 there 

      atend(1)=0
      atend(ny)=0

      end

   
      subroutine greens(dy)

c     1D projection function subroutine
c     construct projection function g (green's fn-type object)
c     convolve it with forcing function q to produce response y
c     domain is periodic

      implicit real(a-h,o-z)      
      parameter(ny=1000) !grid points
      parameter(nx=1000,n=nx,n2=n/2) !grid points
      parameter(pi=3.1415926)

      common /grid/yl,y(ny),yh(ny)
      common /greensfn/g(n,n)
      real xl,x(n),xp(n)

c     set up grid
      
      xl=yl
      dx=dy

      do i=1,n
         x(i)=y(i)
         xp(i)=y(i)+0.5*dy
      end do


C     now we construct the green's function;  j is index for x'

      do j=1,n

c     separate cases of positive and negative x';  it is just
c     easier to handle the periodic BC that way
         if(xp(j).gt.0)then

            do i=j+1,n
               call gfunc(g1,g2,x(i),xp(j))
               g(j,i)=g2
            end do
         
            do i=j-n2+1,j
               call gfunc(g1,g2,x(i),xp(j))
               g(j,i)=g1
            end do

c     now we handle wrap-around points
c     xb is fake wrap-around value of x

            xb=-2.*xl+xp(j)
         
            do i=1,j-n2
               call gfunc(g1,g2,x(i),xb)
               g(j,i)=g2
            end do

         else

            do i=1,j
               call gfunc(g1,g2,x(i),xp(j))
               g(j,i)=g1
            end do

            do i=j+1,j+n2
               call gfunc(g1,g2,x(i),xp(j))
               g(j,i)=g2
            end do

c     now wrap-around points for this case

            xb=2.*xl+xp(j)

            do i=j+n2+1,n
               call gfunc(g1,g2,x(i),xb)
               g(j,i)=g1
            end do

         end if

      end do


      end

      subroutine project(q,u,dx)

c     project precip "forcing" to obtain wind
c     "response" using previously computed projection
c     function g

      implicit real(a-h,o-z)      
      parameter(ny=1000) !grid points
      parameter(nx=1000,n=nx,n2=n/2) !grid points
      parameter(pi=3.1415926)

      common /grid/yl,y(ny),yh(ny)
      common /greensfn/g(n,n)
      real q(n),u(n)


c     do the convolution integral

   
      do i=1,n
         u(i)=0.
         do j=1,n
            u(i)=u(i)+g(j,i)*q(j)*dx
         end do
      end do

      end

      subroutine gfunc(g1,g2,x,xp)

c     This subroutine contains the projection (greens)
c     function itself.  We simply compute it and fill the
c     arrays g1, g2.

c      parameter(gamp=0.6*86400/1.5e6) !equatorial forcing
c      parameter(gamp=1.4*86400/1.5e6) !equatorial forcing
c      parameter(gamp=1.*86400/1.5e6) !equatorial forcing, l=1
      parameter(gamp=0.8*86400/1.5e6) !equatorial forcing,l=2

      !convert m/s to L_R/d
      real g1,g2,x,xp
      real fac,xsc,faci     
 
c     fac=5
c     xsc=3.

      fac=3.  !for equatorial Gill forcing
      xsc=1.

      g1=gamp*fac*exp(fac*(x-xp)/xsc)
      g2=-gamp*exp(-(x-xp)/xsc)

      end

      subroutine shiftu(u,nf)

c     this subroutine shifts u forward set number of grid points

      implicit real(a-h,o-z)      
      parameter(ny=1000) !grid points
      parameter(nx=1000,n=nx,n2=n/2) !grid points
      parameter(pi=3.1415926)

      common /grid/yl,y(ny),yh(ny)
      real u(n),ushift(n)

      do i=1,n
         ushift(i)=u(i)
      end do

      do i=1,n
         if(i.le.nf)then
            u(i)=ushift(n+i-nf)
         else
            u(i)=ushift(i-nf)
         end if
      end do

      end
