module head
    implicit none

    !   6 3 5
    !    \|/
    !   2-0-1
    !    /|\
    !   7 4 8

	integer,parameter::ntmax= 400
	real,parameter:: q=1.5 !p*g^q-1  !parameters of non-Newtonian fluid
	real:: p
	real(8)::  x(-1:ntmax),y(-1:ntmax)
	real(8)::  cenx(-1:ntmax),ceny(-1:ntmax)
	real(8)::  vol(-1:ntmax,-1:ntmax)
	real(8)::  lenx(-1:ntmax),leny(-1:ntmax)
	real(8)::  lenx1(-1:ntmax),leny1(-1:ntmax)
	real(8)::  cx(0:ntmax),cx1(0:ntmax),cx2(0:ntmax)
	real(8)::  cy(0:ntmax),cy1(0:ntmax),cy2(0:ntmax)

    real(8)::  facex(0:ntmax,0:ntmax,3),taox(0:ntmax,0:ntmax)
    real(8)::  facey(0:ntmax,0:ntmax,3),taoy(0:ntmax,0:ntmax)
    integer::  imax,jmax
	real(8)::  re,uu0,niu,niuNX(-1:ntmax,-1:ntmax),niuNY(-1:ntmax,-1:ntmax),Cs,dt,tmax
	real(8)::  s0,s1,s2,s3,s4,s5,s6,s7,s8,s9

	real(8)::  ww0(0:ntmax,0:ntmax,3),ini(0:ntmax,0:ntmax,3)
	real(8)::  flux0(0:ntmax,0:ntmax,3)
	real(8)::  limx(0:ntmax,0:ntmax,3),limy(0:ntmax,0:ntmax,3)
	real(8)::  tminv(9,9),sm(9),tm(9,9),stmiv(9,9)
	real(8)::  fmomX(1:9,-1:ntmax,-1:ntmax),fmomY(1:9,-1:ntmax,-1:ntmax) ! f in moment space with MRT model
	real(8)::  Sx(-1:ntmax,-1:ntmax,3),Sy(-1:ntmax,-1:ntmax,3),gama_dotX(-1:ntmax,-1:ntmax),gama_dotY(-1:ntmax,-1:ntmax)  ! gama_dot xx , xy, yy
	real(8)::  ex(9),ey(9) ! lattice velocities

end module

program main
   USE head
   implicit none
   integer::i,j
   real(8)::a,b

   re  = 1000.d0
   uu0 = 0.1d0
   ! non-Newtonian fluid
   p = uu0**(2.0-q)/re
   ! Newtonian fluid
   niu = uu0/re
   Cs = 1.0/sqrt(3.0)
   dt = 0.005d0
   tmax = 100000

   call fluxMesh

    !initialization
    do j = 0,jmax
        do i = 0,imax
            ww0(i,j,1) = 1.d0
            ww0(i,j,2) = 0.d0
            ww0(i,j,3) = 0.d0

            ini(i,j,1) = 1.d0
            ini(i,j,2) = 0.d0
            ini(i,j,3) = 0.d0
        end do
    end do
    call LBFS

    call outPuts

end program main

    subroutine LBFS
        USE head
        implicit none
        integer:: i,j,iter
        real(8):: res0,res1,wwt(0:ntmax,0:ntmax,3),ww1(0:ntmax,0:ntmax,3),ww2(0:ntmax,0:ntmax,3)

        open(1,file='res.dat')
        do iter = 1,tmax

          do j =1,jmax-1
            do i =1,imax-1
                wwt(i,j,:) = ww0(i,j,:)
            end do
          end do

          ! rung-kutta code

          call BoundaryCondition(ww0)
          call fluxff

          do j=1,jmax-1
              do i=1,imax-1
                    ww1(i,j,:) = wwt(i,j,:)+(-flux0(i,j,:)*1.d0/vol(i,j))*dt
              end do
          end do

          call BoundaryCondition(ww1)
          call fluxff

          do j=1,jmax-1
              do i=1,imax-1
                    ww2(i,j,:) = 3.d0/4.d0*wwt(i,j,:)+1.d0/4.d0*ww1(i,j,:)+1.d0/4.d0*(-flux0(i,j,:)*1.d0/vol(i,j))*dt
              end do
          end do

          call BoundaryCondition(ww2)
          call fluxff

          do j=1,jmax-1
              do i=1,imax-1
                    ww0(i,j,:) = 1.d0/3.d0*wwt(i,j,:)+2.d0/3.d0*ww2(i,j,:)+2.d0/3.d0*(-flux0(i,j,:)*1.d0/vol(i,j))*dt
              end do
          end do

          if(mod(iter,1) .eq. 0) then
            res0 = 0.d0
            res1 = 0.d0
              do j = 1,jmax-1
                do i = 1,imax-1
                      res0 = res0+sqrt((wwt(i,j,2)-ww0(i,j,2))**2+(wwt(i,j,3)-ww0(i,j,3))**2)
                      res1 = res1+sqrt(ww0(i,j,2)**2+ww0(i,j,3)**2)
                end do
              end do
              res1 = res0/res1
              write(*,*) 'the iteration:: ', iter , 'the residual:: ',res1
              write(*,*) '==================================================================='
              write(1,*) iter,res1
          end if
        end do

        close(1)

    end subroutine

    subroutine outPuts
        USE head
        implicit none
        integer::i,j
        integer::mid

        mid = (imax-1)/2
        open(112,file='cavity.dat')
        write(112,*) 'variables= "x","y","rho","u","v" '
        write(112,*) 'zone i=',imax-1 ,'j=',jmax-1 ,'F=POINT'
        do j=1,jmax-1
            do i=1,imax-1
                write(112,*) cenx(i),ceny(j),ini(i,j,1),ini(i,j,2),ini(i,j,3)
            end do
        end do
        close(112)
        open(10,file='slice.dat')
        write(10,*) 'variables= "u","y"'
        write(10,*) 'zone t="LBFS"'
        do j=1,jmax-1
            write(10,*) ceny(j),ini(mid,j,2)/uu0
        end do
        close(10)

    end subroutine

    subroutine feq(dd,uu,vv,ff)

       implicit none

        real(8):: dd,uu,vv,ff(0:8)
        real(8):: a1,b1,c1,uu2,vv2,uv,u2v2

        a1 = 1.d0/9.d0
        b1 = 1.d0/36.d0
        c1 = a1*4.d0

        uu2 = uu*uu
        vv2 = vv*vv
        u2v2= (uu2+vv2)*3.d0
        uv  = uu*vv*9.d0

        ff(1) = a1*dd*(1.d0+3.d0*uu+3.d0*uu2-1.5d0*vv2)
        ff(2) = a1*dd*(1.d0-3.d0*uu+3.d0*uu2-1.5d0*vv2)
        ff(3) = a1*dd*(1.d0+3.d0*vv+3.d0*vv2-1.5d0*uu2)
        ff(4) = a1*dd*(1.d0-3.d0*vv+3.d0*vv2-1.5d0*uu2)
        ff(5) = b1*dd*(1.d0+3.d0*uu+3.d0*vv+u2v2+uv)
        ff(6) = b1*dd*(1.d0-3.d0*uu+3.d0*vv+u2v2-uv)
        ff(7) = b1*dd*(1.d0-3.d0*uu-3.d0*vv+u2v2+uv)
        ff(8) = b1*dd*(1.d0+3.d0*uu-3.d0*vv+u2v2-uv)
        ff(0) = c1*dd*(1.d0-0.5d0*u2v2)

          return
      end subroutine feq

    subroutine feqi(n,dd,uu,vv,ffi)

        implicit none

        integer:: n
        real(8):: dd,uu,vv,ffi
        real(8):: a1,b1,c1,uu2,vv2,uv,u2v2

        a1 = 1.d0/9.d0
        b1 = 1.d0/36.d0
        c1 = a1*4.d0

        uu2 = uu*uu
        vv2 = vv*vv
        u2v2= (uu2+vv2)*3.d0
        uv  = uu*vv*9.d0

        select case (n)
            case (1)
                ffi = a1*dd*(1.d0 + 3.d0*uu + 3.d0*uu2 - 1.5d0*vv2)
            case (2)
                ffi = a1*dd*(1.d0 - 3.d0*uu + 3.d0*uu2 - 1.5d0*vv2)
            case (3)
                ffi = a1*dd*(1.d0 + 3.d0*vv + 3.d0*vv2 - 1.5d0*uu2)
            case (4)
                ffi = a1*dd*(1.d0 - 3.d0*vv + 3.d0*vv2 - 1.5d0*uu2)
            case (5)
                ffi = b1*dd*(1.d0 + 3.d0*uu + 3.d0*vv + u2v2 + uv)
            case (6)
                ffi = b1*dd*(1.d0 - 3.d0*uu + 3.d0*vv + u2v2 - uv)
            case (7)
                ffi = b1*dd*(1.d0 - 3.d0*uu - 3.d0*vv + u2v2 + uv)
            case (8)
                ffi = b1*dd*(1.d0 + 3.d0*uu - 3.d0*vv + u2v2 - uv)
            case (0)
                ffi = c1*dd*(1.d0 - 0.5d0*u2v2)
            case default
                print *, "Error: Invalid value of n"
        end select

        return
    end subroutine feqi

    subroutine fluxff
      USE head
      implicit none
        integer i,j,k,nn,m
          real*8  dd,uu,vv,tao,ss0,f1,f2,f3
          real*8  ldx1,ldy1,lux1,luy1,lvx1,lvy1,lwx1,lwy1
          real*8  ldx2,ldy2,lux2,luy2,lvx2,lvy2,lwx2,lwy2
          real*8  tempd,tempu,tempv,temd,temu,temv
          real*8  dx,dx1,dx2,dy,dy1,dy2
          real*8  dd1,uu1,vv1,dd2,uu2,vv2
          real*8  dd3,uu3,vv3,dd4,uu4,vv4
	      real*8  feq0(0:8),ff(0:8),df,pp,c11,c22

	    do j = 0,jmax
              do i = 0,imax
                flux0(i,j,1) = 0.d0
                flux0(i,j,2) = 0.d0
                flux0(i,j,3) = 0.d0
              end do
	    end do

          !flux in the x direction
	    do j = 1,jmax-1
              ss0 = leny(j)
	        do i = 1,imax
                  dy = facey(i,j,1)
                  dx = facex(i,j,1) ! the slen = e*dt : the streaming distance
                  dx1= facex(i,j,2) ! before the face distance
                  dx2= facex(i,j,3) ! after the face distance

                  dd1 = ini(i-1,j,1)
                  uu1 = ini(i-1,j,2)
                  vv1 = ini(i-1,j,3)

                  ldx1 = limx(i-1,j,1)
                  ldy1 = limy(i-1,j,1)*dy
                  lux1 = limx(i-1,j,2)
                  luy1 = limy(i-1,j,2)*dy
                  lvx1 = limx(i-1,j,3)
                  lvy1 = limy(i-1,j,3)*dy

                  dd2 = ini(i,j,1)
                  uu2 = ini(i,j,2)
                  vv2 = ini(i,j,3)

                  ldx2 = limx(i,j,1)
                  ldy2 = limy(i,j,1)*dy
                  lux2 = limx(i,j,2)
                  luy2 = limy(i,j,2)*dy
                  lvx2 = limx(i,j,3)
                  lvy2 = limy(i,j,3)*dy

                  !left three points
                  tempd = ldx1*dx1
                  tempu = lux1*dx1
                  tempv = lvx1*dx1

                  nn = 1
                  dd3 = dd1 + tempd
                  uu3 = uu1 + tempu
                  vv3 = vv1 + tempv
                  call feqi(nn,dd3,uu3,vv3,ff(nn))

                  nn = 5
                  dd = dd3 - ldy1
                  uu = uu3 - luy1
                  vv = vv3 - lvy1
                  call feqi(nn,dd,uu,vv,ff(nn))

                  nn = 8
                  dd = dd3 + ldy1
                  uu = uu3 + luy1
                  vv = vv3 + lvy1
                  call feqi(nn,dd,uu,vv,ff(nn))

                  !right three points
                  tempd = ldx2*dx2
                  tempu = lux2*dx2
                  tempv = lvx2*dx2

                  nn = 2
                  dd3 = dd2 + tempd
                  uu3 = uu2 + tempu
                  vv3 = vv2 + tempv
                  call feqi(nn,dd3,uu3,vv3,ff(nn))

                  nn = 7
                  dd = dd3 + ldy2
                  uu = uu3 + luy2
                  vv = vv3 + lvy2
                  call feqi(nn,dd,uu,vv,ff(nn))

                  nn = 6
                  dd = dd3 - ldy2
                  uu = uu3 - luy2
                  vv = vv3 - lvy2
                  call feqi(nn,dd,uu,vv,ff(nn))

                  !middle three points
                  dx1 =  lenx1(i-1)
                  dx2 = -lenx1(i)

                  dd3 = dd1 + ldx1*dx1
                  uu3 = uu1 + lux1*dx1
                  vv3 = vv1 + lvx1*dx1
                  dd4 = dd2 + ldx2*dx2
                  uu4 = uu2 + lux2*dx2
                  vv4 = vv2 + lvx2*dx2

                  nn = 0
                  dd3 = (dd3 + dd4)*0.5d0
                  uu3 = (uu3 + uu4)*0.5d0
                  vv3 = (vv3 + vv4)*0.5d0
                  call feqi(nn,dd3,uu3,vv3,ff(nn))

                  ldy1 = 0.5d0*(ldy1 + ldy2)
                  luy1 = 0.5d0*(luy1 + luy2)
                  lvy1 = 0.5d0*(lvy1 + lvy2)

                  nn  = 3
	              dd  = dd3 - ldy1
                  uu  = uu3 - luy1
	              vv  = vv3 - lvy1
                  call feqi(nn,dd,uu,vv,ff(nn))

                  nn  = 4
	              dd  = dd3 + ldy1
                  uu  = uu3 + luy1
	              vv  = vv3 + lvy1
                  call feqi(nn,dd,uu,vv,ff(nn))

                  !end of points calculating the macros
                  dd = ff(1)+ff(2)+ff(3)+ff(4)+ff(5)+ff(6)+ff(7)+ff(8)+ff(0)
	              uu = 1.d0/dd*(ff(1)-ff(2)+ff(5)-ff(6)-ff(7)+ff(8))
                  vv = 1.d0/dd*(ff(3)-ff(4)+ff(5)+ff(6)-ff(7)-ff(8))
                  call feq(dd,uu,vv,feq0)

		          tao = taox(i,j)
			      do k = 0,8
                      ff(k) = feq0(k)-(feq0(k)-ff(k))*1.2
	              end do

                  f1 = (feq0(1)-feq0(2)+feq0(5)-feq0(6)-feq0(7)+feq0(8))*ss0
                  f2 = (ff(1)+ff(2)+ff(5)+ff(6)+ff(7)+ff(8))*ss0
                  f3 = (ff(5)-ff(6)+ff(7)-ff(8))*ss0

                  flux0(i-1,j,1) = flux0(i-1,j,1)+f1
                  flux0(i-1,j,2) = flux0(i-1,j,2)+f2
                  flux0(i-1,j,3) = flux0(i-1,j,3)+f3

                  flux0(i,j,1) = flux0(i,j,1)-f1
                  flux0(i,j,2) = flux0(i,j,2)-f2
                  flux0(i,j,3) = flux0(i,j,3)-f3
              end do
          end do

          !flux in y direction
          do i = 1,imax-1
	        ss0 = lenx(i)
	        do j = 1,jmax
                  dx = facex(i,j,1)
                  dy = facey(i,j,1)
                  dy1= facey(i,j,2)
	              dy2= facey(i,j,3)

                  dd1 = ini(i,j-1,1)
	              uu1 = ini(i,j-1,2)
                  vv1 = ini(i,j-1,3)

                  ldx1 = limx(i,j-1,1)*dx
                  ldy1 = limy(i,j-1,1)
                  lux1 = limx(i,j-1,2)*dx
                  luy1 = limy(i,j-1,2)
                  lvx1 = limx(i,j-1,3)*dx
                  lvy1 = limy(i,j-1,3)

                  dd2 = ini(i,j,1)
	              uu2 = ini(i,j,2)
                  vv2 = ini(i,j,3)

                  ldx2 = limx(i,j,1)*dx
                  ldy2 = limy(i,j,1)
                  lux2 = limx(i,j,2)*dx
                  luy2 = limy(i,j,2)
                  lvx2 = limx(i,j,3)*dx
                  lvy2 = limy(i,j,3)

                  !bottom three points
                  tempd = ldy1*dy1
                  tempu = luy1*dy1
                  tempv = lvy1*dy1

                  nn = 3
                  dd3 = dd1 + tempd
                  uu3 = uu1 + tempu
                  vv3 = vv1 + tempv
                  call feqi(nn,dd3,uu3,vv3,ff(nn))

                  nn = 5
                  dd = dd3 - ldx1
                  uu = uu3 - lux1
                  vv = vv3 - lvx1
                  call feqi(nn,dd,uu,vv,ff(nn))

                  nn = 6
                  dd = dd3 + ldx1
                  uu = uu3 + lux1
                  vv = vv3 + lvx1
                  call feqi(nn,dd,uu,vv,ff(nn))

                  !top three points
                  tempd = ldy2*dy2
                  tempu = luy2*dy2
                  tempv = lvy2*dy2

                  nn = 4
                  dd3 = dd2 + tempd
                  uu3 = uu2 + tempu
                  vv3 = vv2 + tempv
                  call feqi(nn,dd3,uu3,vv3,ff(nn))

                  nn = 7
                  dd = dd3 + ldx2
                  uu = uu3 + lux2
                  vv = vv3 + lvx2
                  call feqi(nn,dd,uu,vv,ff(nn))

                  nn = 8
                  dd = dd3 - ldx2
                  uu = uu3 - lux2
                  vv = vv3 - lvx2
                  call feqi(nn,dd,uu,vv,ff(nn))

                  !middle three points
                  dy1 =  leny1(j-1)
                  dy2 = -leny1(j)

                  dd3 = dd1 + ldy1*dy1
                  uu3 = uu1 + luy1*dy1
                  vv3 = vv1 + lvy1*dy1
                  dd4 = dd2 + ldy2*dy2
                  uu4 = uu2 + luy2*dy2
                  vv4 = vv2 + lvy2*dy2

                  nn = 0
                  dd3 = (dd3 + dd4)*0.5d0
                  uu3 = (uu3 + uu4)*0.5d0
                  vv3 = (vv3 + vv4)*0.5d0
                  call feqi(nn,dd3,uu3,vv3,ff(nn))

                  ldx1 = 0.5d0*(ldx1 + ldx2)
                  lux1 = 0.5d0*(lux1 + lux2)
                  lvx1 = 0.5d0*(lvx1 + lvx2)

                  nn  = 1
	              dd  = dd3 - ldx1
                  uu  = uu3 - lux1
	              vv  = vv3 - lvx1
                  call feqi(nn,dd,uu,vv,ff(nn))

                  nn  = 2
	              dd  = dd3 + ldx1
                  uu  = uu3 + lux1
	              vv  = vv3 + lvx1
                  call feqi(nn,dd,uu,vv,ff(nn))

                  !end of points
                  dd = ff(1)+ff(2)+ff(3)+ff(4)+ff(5)+ff(6)+ff(7)+ff(8)+ff(0)
	              uu = 1.d0/dd*(ff(1)-ff(2)+ff(5)-ff(6)-ff(7)+ff(8))
	              vv = 1.d0/dd*(ff(3)-ff(4)+ff(5)+ff(6)-ff(7)-ff(8))
                  call feq(dd,uu,vv,feq0)

                  tao = taoy(i,j)
                  do k = 0,8
                    ff(k) = feq0(k)-(feq0(k)-ff(k))*1.2
                  end do

                  f1 = (feq0(3)-feq0(4)+feq0(5)+feq0(6)-feq0(7)-feq0(8))*ss0
                  f2 = (ff(5)-ff(6)+ff(7)-ff(8))*ss0
                  f3 = (ff(3)+ff(4)+ff(5)+ff(6)+ff(7)+ff(8))*ss0

                  flux0(i,j-1,1) = flux0(i,j-1,1)+f1
                  flux0(i,j-1,2) = flux0(i,j-1,2)+f2
                  flux0(i,j-1,3) = flux0(i,j-1,3)+f3

                  flux0(i,j,1) = flux0(i,j,1)-f1
                  flux0(i,j,2) = flux0(i,j,2)-f2
                  flux0(i,j,3) = flux0(i,j,3)-f3
	        end do
          end do

          return
    end subroutine fluxff

    subroutine fluxMesh

      USE head
      implicit none

      integer::i,j
      real(8)::dx,dy,slen,dx1,dy1

      call gridgen

      do i = 0,imax
          lenx(i)=x(i+1) - x(i)
      end do

      do j = 0,jmax
          leny(j)=y(j+1) - y(j)
      end do

      do i = 0,imax
          lenx1(i) = 0.5d0*lenx(i)
      end do

      do j = 0,jmax
          leny1(j) = 0.5d0*leny(j)
      end do

      do j = 0,jmax
        dy = leny(j)
        do i = 0,imax
            dx = lenx(i)
            vol(i,j) = (dx*dy)
        end do
      end do

      do j = 1,jmax-1
          dy = leny(j)
          do i = 1, imax
              dx = lenx(i-1)
              dx1= lenx(i)
              slen= 0.25d0*min(dx,dx1,dy) ! dx=dt for lattices
              facex(i,j,1)= slen
              facex(i,j,2)= 0.5d0*dx - slen
              facex(i,j,3)= slen - 0.5d0*dx1
              taox(i,j) = uu0/(re*slen*Cs**2)
          end do
      end do

      do i = 1,imax-1
          dx= lenx(i)
          do j = 1,jmax
              dy   = leny(j-1)
              dy1  = leny(j)
              slen = 0.25d0*min(dx,dy,dy1)! dx=dt for lattices
              facey(i,j,1)= slen
              facey(i,j,2)= 0.5d0*dy - slen
              facey(i,j,3)= slen - 0.5d0*dy1
              taoy(i,j) = uu0/(re*slen*Cs**2)
          end do
       end do

      do i = 2,imax-1
          cx(i) = 0.5d0*(lenx(i-1)+lenx(i))
      end do
      cx(1)    = lenx(1)*0.5d0
      cx(imax) = lenx(imax-1)*0.5d0

      do i = 1,imax-1
          cx1(i) = cx(i)/(cx(i+1)*(cx(i)+cx(i+1)))
          cx2(i) = cx(i+1)/(cx(i)*(cx(i)+cx(i+1)))
      end do

      do j = 2,jmax-1
          cy(j) = 0.5d0*(leny(j-1)+leny(j))
      end do
      cy(1)    = leny(1)*0.5d0
      cy(jmax) = leny(jmax-1)*0.5d0

      do j = 1,jmax-1
          cy1(j) = cy(j)/(cy(j+1)*(cy(j)+cy(j+1)))
          cy2(j) = cy(j+1)/(cy(j)*(cy(j)+cy(j+1)))
      end do

      open(112,file='cavity_mesh.dat')
      write(112,*) 'variables= "x","y", '
      write(112,*) 'zone i=',imax-1 ,'j=',jmax-1 ,'F=POINT'
      do j=1,jmax-1
        do i=1,imax-1
            write(112,*) cenx(i),ceny(j)
        end do
      end do
      close(112)

      return
    end subroutine fluxMesh

    subroutine gridgen

    USE head
    implicit none

      integer::i,j
      real(8) ::rp,pi,eta,rlg_x,rlg_y,disr,tmp

        ! Specify dimensions
        imax = 101   ! Width (number of points in x direction)
        jmax = 101   ! Height (number of points in y direction)

        pi  = 4.d0 * atan(1.d0)
        eta = 50000.0d0

        ! Two separate coefficients for x and y directions
        rlg_x = 0.5d0   ! Coefficient for x-direction
        rlg_y = 0.5d0   ! Coefficient for y-direction

        ! Generate x-coordinates (width direction) with rlg_x
        do i = 1, ceiling(real(imax)/real(2))
            tmp  = atan((1.d0 - real(i-1)/((imax-1)/2)) * dtan(1.d0 / eta))
            disr = rlg_x * (1.d0 - eta * tmp)
            x(i) = disr
        end do

        do i = imax, ceiling(real(imax)/real(2)) + 1, -1
            x(i) = 2.d0 * x(ceiling(real(imax)/real(2))) - x(imax + 1 - i)
        end do

        ! Generate y-coordinates (height direction) with rlg_y
        do j = 1, ceiling(real(jmax)/real(2))
            tmp  = atan((1.d0 - real(j-1)/((jmax-1)/2)) * dtan(1.d0 / eta))
            disr = rlg_y * (1.d0 - eta * tmp)
            y(j) = disr
        end do

        do j = jmax, ceiling(real(jmax)/real(2)) + 1, -1
            y(j) = 2.d0 * y(ceiling(real(jmax)/real(2))) - y(jmax + 1 - j)
        end do

        ! Boundary for ghost nodes (adjust for both x and y directions)
        x(0)       = 2.d0 * x(1) - x(2)
        x(imax+1)  = 2.d0 * x(imax) - x(imax-1)

        y(0)       = 2.d0 * y(1) - y(2)
        y(jmax+1)  = 2.d0 * y(jmax) - y(jmax-1)


      do i = 0,imax
          cenx(i) = (x(i)+x(i+1))*0.5d0
      end do
      do j = 0,jmax
          ceny(j) = (y(j)+y(j+1))*0.5d0
      end do

      return
  end subroutine gridgen

    subroutine gradients
    use head
    implicit none
    integer:: i,j

     !derivatives in the x and y directions
      do j = 1,jmax-1
        do i = 1,imax-1
              limx(i,j,:) = cx1(i)*(ini(i+1,j,:)-ini(i,j,:))+cx2(i)*(ini(i,j,:)-ini(i-1,j,:))
              limy(i,j,:) = cy1(j)*(ini(i,j+1,:)-ini(i,j,:))+cy2(j)*(ini(i,j,:)-ini(i,j-1,:))
        end do
      end do

      !boundary conditions for gradients in cell borders
      i = 0
      do j = 1,jmax-1

          limx(i,j,1) =    (ini(i+1,j,1)-ini(i,j,1))/lenx(i)
          limx(i,j,2) =    (ini(i+1,j,2)-ini(i,j,2))/lenx(i)
          limx(i,j,3) =    (ini(i+1,j,3)-ini(i,j,3))/lenx(i)

      end do

      i = imax
      do j = 1,jmax-1

          limx(i,j,1) =    (ini(i,j,1)-ini(i-1,j,1))/lenx(i)
          limx(i,j,2) =    (ini(i,j,2)-ini(i-1,j,2))/lenx(i)
          limx(i,j,3) =    (ini(i,j,3)-ini(i-1,j,3))/lenx(i)

      end do

      j = 0
      do i = 1,imax-1

          limy(i,j,1) =    (ini(i,j+1,1)-ini(i,j,1))/leny(j)
          limy(i,j,2) =    (ini(i,j+1,2)-ini(i,j,2))/leny(j)
          limy(i,j,3) =    (ini(i,j+1,3)-ini(i,j,3))/leny(j)

      end do

      j = jmax
      do i = 1,imax-1
          !ini(i,j,2)  =  2.d0*uu0 - ini(i,j-1,2)
          !ini(i,j,3)  =  - ini(i,j-1,3)

          limy(i,j,1) =   (ini(i,j,1)-ini(i,j-1,1))/leny(j)
          limy(i,j,2) =   (ini(i,j,2)-ini(i,j-1,2))/leny(j)
          limy(i,j,3) =   (ini(i,j,3)-ini(i,j-1,3))/leny(j)

      end do

      return

    end subroutine gradients

    subroutine BoundaryCondition(ww)

    USE head
    implicit none

    integer:: i,j
    real(8) :: dd
    real(8), intent(IN) :: ww(0:ntmax,0:ntmax,3)

	    do j = 1,jmax-1
              do i = 1,imax-1
                  dd = ww(i,j,1)
                  ini(i,j,1) = dd
                  ini(i,j,2) = 1.d0/dd*ww(i,j,2)
                  ini(i,j,3) = 1.d0/dd*ww(i,j,3)
              end do
          end do

          ! node base calculation we dont store any flux variable and all ini variables are node based
          ! boundary conditions for flow field variables
          ! Neumann for density and no slip BC for velocities on walls
          ! note that the BC has been implemented on the ghost nodes to find the fluxes
          ! note that we do not store the flux values in any variable

          i = 0
          do j = 1,jmax-1
              ini(i,j,1) = ini(i+1,j,1)
              ini(i,j,2)  =  0.0d0
              ini(i,j,3)  =  - ini(i+1,j,3)
          end do

          i = imax
          do j = 1,jmax-1
              ini(i,j,1) = ini(i-1,j,1)
              ini(i,j,2)  =  0.0d0
              ini(i,j,3)  =  - ini(i-1,j,3)
          end do

          j = 0
          do i = 1,imax-1
              ini(i,j,1) = ini(i,j+1,1)
              ini(i,j,2)  =  - ini(i,j+1,2)
              ini(i,j,3)  =  0.0d0
          end do

          j = jmax
          do i = 1,imax-1
              ini(i,j,1) = ini(i,j-1,1)
              ini(i,j,2)  =  2.d0*uu0 - ini(i,j-1,2)
              ini(i,j,3)  =  0.0d0
          end do

      ! calculating the field variables after each time step updating
      call gradients

      return
  end subroutine BoundaryCondition




