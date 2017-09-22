      program greens

c     test program for 1D projection function subroutine
c     construct projection function g (green's fn-type object)
c     convolve it with forcing function q to produce response y
c     domain is periodic

      implicit real(a-h,o-z)      
      parameter(nx=1001,n=nx-1,n2=n/2) !grid points
      parameter(xsc=1.) !decay scale for green's fn
      parameter(fac=3)  !difference in magnitude of e'lies & w'lies
      parameter(pi=3.1415926)

      common /grid/xl,x(n),xp(n)

      dimension g(n,n),q(n),y(n)
  
c     set up grid
      
      xl=10.
      dx=2*xl/float(n)

      do i=1,n
         x(i)=-xl+(float(i)-0.5)*dx
         xp(i)=x(i)
      end do
c      print*,x(1),x(n)
c     construct a forcing function
      sig=0.2
      do i=1,n
         if(abs(x(i)).lt.1)then
            q(i)=sin(pi*x(i))
         else
            q(i)=0.
         end if

          q(i)=10
c         q(i)=exp(-(x(i)*x(i))/(sig*sig))

      end do

c     now we construct the green's function;  j is index for x'

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

c     now we have the projection function, let's do the 
c     convolution integral

   
      do i=1,n
         y(i)=0.
         do j=1,n
            y(i)=y(i)+g(j,i)*q(j)*dx
         end do
      end do
      
      open(1,file='gf.dat',form='formatted')
      do i=1,n
         write(1,*)x(i),q(i),y(i)
      end do

      open(2,file='g.dat',form='unformatted',
     &     access='direct',recl=n*4)
      irec=1
      do j=1,n
         write(2,rec=irec)(g(j,i),
     &        i=1,n)
         irec=irec+1
      end do
      close(2)

      end

      subroutine gfunc(g1,g2,x,xp)

      parameter(gamp=8.)
      real g1,g2,x,xp
      real fac,xsc,faci
      fac=10.
      xsc=1.
      faci=1./fac
      
      g1=gamp*exp(fac*(x-xp))
      g2=-gamp*faci*exp(-(x-xp))

      end
