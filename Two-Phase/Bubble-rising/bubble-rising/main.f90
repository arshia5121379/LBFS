module head
    implicit none

    !   6 3 5             8 4 7
    !    \|/               \|/
    !   2-9-1 streaming   1-9-2 collision
    !    /|\               /|\
    !   7 4 8             5 3 7

	integer,parameter::ntmax= 500
    integer::ii
	real,parameter:: q=1.5 !p*g^q-1  !parameters of non-Newtonian fluid
	real:: p
	real(8)::  x(-1:ntmax),y(-1:ntmax)
	real(8)::  cenx(-1:ntmax),ceny(-1:ntmax)
	real(8)::  vol(-1:ntmax,-1:ntmax)
	real(8)::  lenx(-1:ntmax),leny(-1:ntmax)
	real(8)::  lenx1(-1:ntmax),leny1(-1:ntmax)
	real(8)::  cx(0:ntmax),cx1(0:ntmax),cx2(0:ntmax)
	real(8)::  cy(0:ntmax),cy1(0:ntmax),cy2(0:ntmax)

    real(8)::  facex(0:ntmax,0:ntmax,4),taoxH(0:ntmax,0:ntmax),taoxL(0:ntmax,0:ntmax),taogx(0:ntmax,0:ntmax)
    real(8)::  facey(0:ntmax,0:ntmax,4),taoyH(0:ntmax,0:ntmax),taoyL(0:ntmax,0:ntmax),taogy(0:ntmax,0:ntmax)
    integer::  imax,jmax,file_unit,iter
	real(8)::  uu0,niu,niuNX(-1:ntmax,-1:ntmax),niuNY(-1:ntmax,-1:ntmax),Cs,dt,tmax
	real(8)::  s0,s1,s2,s3,s4,s5,s6,s7,s8,s9
    ! 1=>phi , 2=>p, 3=>u , 4=>v, 5=>density , 6=>miu
	real(8)::  ww0(0:ntmax,0:ntmax,4),ini(0:ntmax,0:ntmax,6) !number 5 for density and number 6 for chemical function & ww0 for fluxes variables
	real(8)::  flux0(0:ntmax,0:ntmax,4)
	real(8)::  limx(0:ntmax,0:ntmax,6),limy(0:ntmax,0:ntmax,6)
	real(8)::  tminv(9,9),sm(9),tm(9,9),stmiv(9,9)
	real(8)::  fmomX(1:9,-1:ntmax,-1:ntmax),fmomY(1:9,-1:ntmax,-1:ntmax) ! f in moment space with MRT model
	real(8)::  Sx(-1:ntmax,-1:ntmax,3),Sy(-1:ntmax,-1:ntmax,3),gama_dotX(-1:ntmax,-1:ntmax),gama_dotY(-1:ntmax,-1:ntmax)  ! gama_dot xx , xy, yy
	real(8)::  ex(9),ey(9),wa(9) ! lattice velocities
	real(8)::  rho_facex(-1:ntmax,-1:ntmax),rho_facey(-1:ntmax,-1:ntmax) ! this is for the density on the cell interfaces through solving order parameter for old time step
	real(8)::  rho_facexO(-1:ntmax,-1:ntmax),rho_faceyO(-1:ntmax,-1:ntmax)
	real(8)::  ddL,ddH,niuH,niuL,miuL,miuH
	real(8)::  Eo,Re,miu_ratio,dd_ratio
    real(8)::  rr,dia,pi,dtt,M,beta,epsilon,sigma,etha,kapa,gy ! the parameters of chan-hilliard Eqn.
    real(8)::  Fsx(-1:ntmax,-1:ntmax),Fsy(-1:ntmax,-1:ntmax)
    real(8)::  limx2(-1:ntmax,-1:ntmax),limy2(-1:ntmax,-1:ntmax) ! laplacian scheme
    real(8)::  MG0
    real(8)::  Us(ntmax),Ys(ntmax),Ts(ntmax)
end module

module rho_0
    implicit none
      integer:: i,j,k,nn,kk
      real(8)::  dd,miu,uu,vv,phi,taof,taog,ss0,f1,f2,f3,f4
      real(8)::  ldx1,ldy1,lux1,luy1,lvx1,lvy1,lwx1,lwy1
      real(8)::  ldx2,ldy2,lux2,luy2,lvx2,lvy2,lwx2,lwy2
      real(8)::  lpx1,lpy1,lmx1,lmy1,lpx2,lpy2,lmx2,lmy2
      real(8)::  lphix1,lphiy1,lphix2,lphiy2
      real(8)::  lmiux1,lmiuy1,lmiux2,lmiuy2
      real(8)::  tempd,tempu,tempv,temd,temu,temv,tempmiu,tempp,tempphi
      real(8)::  dx,dx1,dx2,dy,dy1,dy2
      real(8)::  pp1,uu1,vv1,pp2,uu2,vv2,miu1,miu2,phi1,phi2,dd1,dd2
      real(8)::  pp3,uu3,vv3,pp4,uu4,vv4,miu3,miu4,phi3,phi4,dd3,dd4
      real(8)::  ddf1,ddf2,ddf3,ddf4,ddf
      real(8)::  feq0(9),ffi(9),ff(9),ggi(9),fff(9),df,pp,c11,c22
      real(8)::  tmp1,tmp2,tmp3,tmp4
end module

    program main
       USE head
       implicit none
       integer::i,j

       dt = 0.1d0
       tmax = 380000

       call gridgen
       call fluxMesh
       call initialize
       call LBFS
       call outPuts

    end program main

    subroutine initialize
        USE head
        USE rho_0
        implicit none
        integer:: in,ip,jp,jn
        real(8):: factor,cs2,VL


        wa(9)=4.d0/9.d0
        do i=1,4
            wa(i)=1.d0/9.d0
        end do
        do i=5,8
            wa(i)=1.d0/36.d0
        end do

        ! for normalizing in the lattices coordination :
        ! U_pys = 1 m/s and U_LBM = 0.1 for more stability
        ! and the D_phy = 1mm and the L_LBM = 320 => L_factor = 1.25e-5
        ! dt_factor = U_LBM/U_phy * L_factor
        ! ===========================================================================================================
        ! **gy and sigma have to estimated from both Re and Eo number that ensure the consistency between those numbers**
        ! ===========================================================================================================
        Re = 35.0
        Eo = 10.0
        ddH  = 1.d0
        miuH = 0.015d0
        dtt = maxval(lenx)
        etha= 10.d0
        rr=40.0d0*dtt  ! radius of the drop
        dia=rr*2.0
        Cs = 1.d0/sqrt(3.d0)
        pi=4.0*atan(1.0)
        epsilon=4.d0*dtt  ! interface thickness
        gy = -(1.d0/dia**3)*((Re*miuH)/(ddH))**2
        sigma = ddH*abs(gy)*dia**2/Eo  ! surface tension
        kapa=1.5*sigma*epsilon
        beta=12.0*sigma/epsilon
        miu_ratio=10.0  ! the viscosity ratio
        dd_ratio=20.0  ! density ratio
        miuL = miuH/miu_ratio
        ddL= ddH/dd_ratio
        niuH = miuH/ddH
        niuL = miuL/ddL

        do i=0,imax
            do j=0,jmax
                ! initial order parameter for "phi"
                factor = sqrt((x(i)-80.0d0)**2+(y(j)-80.0d0)**2)
                ini(i,j,1) = 0.5d0*(1.d0-tanh((dia/2.d0-factor)/(0.5*epsilon)))
                ww0(i,j,1) = ini(i,j,1)

                ini(i,j,3) = 0.0 ! initial velocity for u
                ini(i,j,4) = 0.0 ! initial velocity for v
                ww0(i,j,3) = ini(i,j,3)
                ww0(i,j,4) = ini(i,j,4)

                ini(i,j,5) = ddH*ini(i,j,1)+(1.d0-ini(i,j,1))*ddL ! initial density
                ini(i,j,2) = (1.d0-ini(i,j,1))*4.d0*sigma/dia  ! initial pressure
                ww0(i,j,2) = ini(i,j,2)

            end do
        end do


        ! calculating the MG0 at as a initial gas mass :
        VL = 0.0d0
        do i=1,imax-1
            do j=1,jmax-1
                if (ini(i,j,1) >= 0 .AND. ini(i,j,1) <= 0.5 ) then
                    VL = VL + vol(i,j)
                end if
            end do
        end do
        MG0 = ddL*VL

      ! calculating the chemical potential function MiuC:
      ! miuc = 2A(PHI)*(PHI-1)*(2PHI-1)-kD2PHI
      ! call Laplacian

          do i=1,imax-1
            do j=1,jmax-1
                ip=i-1
                in=i+1
                jp=j-1
                jn=j+1
                ini(i,j,6) = 2.d0*beta*(ini(i,j,1))*(ini(i,j,1)-1.d0)*(2.d0*ini(i,j,1)-1.d0)
                ini(i,j,6) = ini(i,j,6)-&
                &kapa*(ini(in,jn,1)+ini(in,jp,1)+ini(ip,jn,1)+ini(ip,jp,1)+4.0*ini(ip,j,1)+&
                & 4.0*ini(in,j,1)+4.0*ini(i,jp,1)+4.0*ini(i,jn,1)-20.0*ini(i,j,1))/6.0
            end do
           end do

        ! for boundary conditions:
          i = 0
          do j = 1,jmax-1
              ini(i,j,6) = ini(i+1,j,6)
          end do

          i = imax
          do j = 1,jmax-1
              ini(i,j,6) = ini(i-1,j,6)
          end do

          j = 0
          do i = 1,imax-1
              ini(i,j,6) = ini(i,j+1,6)
          end do

          j = jmax
          do i = 1,imax-1
              ini(i,j,6) = ini(i,j-1,6)
          end do


        call gradients
        ! calculate the density and kinematic viscosity  amount on the interfaces for feqm subroutine

        ! in x direction :
        do j=1,jmax-1
            do i=1,imax
              phi1= ini(i-1,j,1)
              phi2= ini(i,j,1)
              lphix1= limx(i-1,j,1)
              lphix2= limx(i,j,1)
              dx1 =  lenx1(i-1)
              dx2 = -lenx1(i)
              phi3= phi1+ lphix1*dx1
              phi4= phi2+ lphix2*dx2
              phi3= (phi3+phi4)*0.5d0
              ddf3= ddH*phi3+(1.d0-phi3)*ddL ! rho (aplha = 9)
              rho_facexO(i,j) = ddf3
            end do
        end do

        ! in y direction
        do i=1,imax-1
            do j=1,jmax
              phi1= ini(i,j-1,1)
              phi2= ini(i,j,1)
              lphiy1= limx(i,j-1,1)
              lphiy2= limx(i,j,1)
              dy1 =  leny1(j-1)
              dy2 = -leny1(j)
              phi3= phi1+ lphiy1*dy1
              phi4= phi2+ lphiy2*dy2
              phi3= (phi3+phi4)*0.5d0
              ddf3= ddH*phi3+(1.d0-phi3)*ddL ! rho (aplha = 9)
              rho_faceyO(i,j) = ddf3
            end do
        end do

        call updateTao
        call tensionForces

        ! calculating the M factor:
        M = (maxval(taogx)-0.5d0)*etha*maxval(facex(:,:,1))

    end subroutine initialize
    ! main program loop
    subroutine LBFS
        USE head
        implicit none
        integer:: i,j,imid,jmid
        real(8):: res0,res1,wwt(0:ntmax,0:ntmax,4),ww1(0:ntmax,0:ntmax,4),ww2(0:ntmax,0:ntmax,4)
        real(8):: Qq(ntmax,ntmax)
        character(len=5)::key

        write(*,*) "select the time discretization : Rk for Runge-Kutta and Eu for Euler explicit"
        read(*,*) key

        if (key=='Rk') then
            ii = 1
            file_unit = 1

            do iter = 1,tmax
              ! storing old variables
              do j =1,jmax-1
                do i =1,imax-1
                    wwt(i,j,1) = ww0(i,j,1) ! phi
                    wwt(i,j,2) = ww0(i,j,2) ! pressure
                    wwt(i,j,3) = ww0(i,j,3) ! rho*Cs^2*u
                    wwt(i,j,4) = ww0(i,j,4) ! rho*Cs^2*v
                end do
              end do

            ! First RK Stage: Update ww1
             call BoundaryCondition(ww0)
             call tensionForces
            !call mass(MG0,Qq)

             rho_facexO = rho_facex
             rho_faceyO = rho_facey

             call fluxff ! Calculate flux at new step

             do j = 1, jmax-1
                do i = 1, imax-1
                    ! Component 1 Update
                    ww1(i,j,1) = wwt(i,j,1) - (flux0(i,j,1) * (1.0d0 / vol(i,j)) * dt)

                    ! Component 2 Update
                    ww1(i,j,2) = wwt(i,j,2) - (flux0(i,j,2) * (1.0d0 / vol(i,j)) * dt) &
                                 + ((ini(i,j,3) * limx(i,j,5) + ini(i,j,4) * limy(i,j,5)) * Cs**2 * dt)

                    ! Component 3 Update
                    ww1(i,j,3) = wwt(i,j,3) - (flux0(i,j,3) * (1.0d0 / vol(i,j)) * dt) &
                                 + (Cs**2 * Fsx(i,j) * dt)

                    ! Component 4 Update
                    ww1(i,j,4) = wwt(i,j,4) - (flux0(i,j,4) * (1.0d0 / vol(i,j)) * dt) &
                                 + ((Cs**2 * Fsy(i,j) + gy * (ini(i,j,5) - ddH)) * dt)
                end do
             end do

            ! Second RK Stage: Update ww2
             call BoundaryCondition(ww1)
             call tensionForces
            !call mass(MG0,Qq)

             rho_facexO = rho_facex
             rho_faceyO = rho_facey

             call fluxff ! Calculate flux at new step

             do j = 1, jmax-1
                do i = 1, imax-1
                    ! Component 1 Update
                    ww2(i,j,1) = 0.75d0 * wwt(i,j,1) + 0.25d0 * ww1(i,j,1) + &
                                 0.25d0 * (-flux0(i,j,1) * (1.0d0 / vol(i,j)) * dt)

                    ! Component 2 Update
                    ww2(i,j,2) = 0.75d0 * wwt(i,j,2) + 0.25d0 * ww1(i,j,2) + &
                                 0.25d0 * (-flux0(i,j,2) * (1.0d0 / vol(i,j)) + &
                                         (ini(i,j,3) * limx(i,j,5) + ini(i,j,4) * limy(i,j,5)) * Cs**2) * dt

                    ! Component 3 Update
                    ww2(i,j,3) = 0.75d0 * wwt(i,j,3) + 0.25d0 * ww1(i,j,3) + &
                                 0.25d0 * (-flux0(i,j,3) * (1.0d0 / vol(i,j)) + Cs**2 * Fsx(i,j)) * dt

                    ! Component 4 Update
                    ww2(i,j,4) = 0.75d0 * wwt(i,j,4) + 0.25d0 * ww1(i,j,4) + &
                                 0.25d0 * (-flux0(i,j,4) * (1.0d0 / vol(i,j)) + &
                                         Cs**2 * Fsy(i,j) + gy * (ini(i,j,5) - ddH)) * dt
                end do
             end do

            ! Third RK Stage: Update ww0
             call BoundaryCondition(ww2)
             call tensionForces
            !call mass(MG0,Qq)

             rho_facexO = rho_facex
             rho_faceyO = rho_facey

             call fluxff ! Calculate flux at new step

             do j = 1, jmax-1
                do i = 1, imax-1
                    ! Component 1 Update
                    ww0(i,j,1) = (1.d0 / 3.d0) * wwt(i,j,1) + &
                                 (2.d0 / 3.d0) * ww2(i,j,1) + &
                                 (2.d0 / 3.d0) * (-flux0(i,j,1) * (1.0d0 / vol(i,j)) * dt)

                    ! Component 2 Update
                    ww0(i,j,2) = (1.d0 / 3.d0) * wwt(i,j,2) + &
                                 (2.d0 / 3.d0) * ww2(i,j,2) + &
                                 (2.d0 / 3.d0) * (-flux0(i,j,2) * (1.0d0 / vol(i,j)) + &
                                                  (ini(i,j,3) * limx(i,j,5) + ini(i,j,4) * limy(i,j,5)) * Cs**2) * dt

                    ! Component 3 Update
                    ww0(i,j,3) = (1.d0 / 3.d0) * wwt(i,j,3) + &
                                 (2.d0 / 3.d0) * ww2(i,j,3) + &
                                 (2.d0 / 3.d0) * (-flux0(i,j,3) * (1.0d0 / vol(i,j)) + Cs**2 * Fsx(i,j)) * dt

                    ! Component 4 Update
                    ww0(i,j,4) = (1.d0 / 3.d0) * wwt(i,j,4) + &
                                 (2.d0 / 3.d0) * ww2(i,j,4) + &
                                 (2.d0 / 3.d0) * (-flux0(i,j,4) * (1.0d0 / vol(i,j)) + &
                                                  Cs**2 * Fsy(i,j) + gy * (ini(i,j,5) - ddH)) * dt
                end do
              end do


              if(mod(iter,500) .eq. 0) then
                res0 = 0.d0
                res1 = 0.d0
                do j = 1,jmax-1
                    do i = 1,imax-1
                      res0 = res0+sqrt((wwt(i,j,3)-ww0(i,j,3))**2+(wwt(i,j,4)-ww0(i,j,4))**2)
                      res1 = res1+sqrt(ww0(i,j,3)**2+ww0(i,j,4)**2)
                    end do
                end do
                res1 = res0/res1
                write(*,*) 'the iteration:: ', iter , 'the residual:: ',res1
                write(*,*) '==================================================================='
                call output_time
              end if
            end do

            call outPuts

        elseif(key=='Eu') then
            ii = 1
            file_unit = 1

            do iter = 1,tmax
              ! storing old variables
              do j =1,jmax-1
                do i =1,imax-1
                    wwt(i,j,1) = ww0(i,j,1) ! phi
                    wwt(i,j,2) = ww0(i,j,2) ! pressure
                    wwt(i,j,3) = ww0(i,j,3) ! rho*Cs^2*u
                    wwt(i,j,4) = ww0(i,j,4) ! rho*Cs^2*v
                end do
              end do

              call BoundaryCondition(ww0)
              call tensionForces
              !call mass(MG0,Qq)

              rho_facexO = rho_facex
              rho_faceyO = rho_facey

              call fluxff ! calculate flux at new step

              do j=1,jmax-1
                  do i=1,imax-1

                    ww0(i,j,1) = wwt(i,j,1)-flux0(i,j,1)*(1.0d0/vol(i,j))*dt!+Qq(i,j)*dt
                    ww0(i,j,2) = wwt(i,j,2)-flux0(i,j,2)*(1.0d0/vol(i,j))*dt&
                    & +(ini(i,j,3)*limx(i,j,5)+ini(i,j,4)*limy(i,j,5))*dt*Cs**2
                    ww0(i,j,3) = wwt(i,j,3)-flux0(i,j,3)*(1.0d0/vol(i,j))*dt+Cs**2*Fsx(i,j)*dt
                    ww0(i,j,4) = wwt(i,j,4)-flux0(i,j,4)*(1.0d0/vol(i,j))*dt+(Cs**2*Fsy(i,j)+gy*(ini(i,j,5)-ddH))*dt

                  end do
              end do

              if(mod(iter,1000) .eq. 0) then
                res0 = 0.d0
                res1 = 0.d0
                do j = 1,jmax-1
                    do i = 1,imax-1
                      res0 = res0+sqrt((wwt(i,j,3)-ww0(i,j,3))**2+(wwt(i,j,4)-ww0(i,j,4))**2)
                      res1 = res1+sqrt(ww0(i,j,3)**2+ww0(i,j,4)**2)
                    end do
                end do
                res1 = res0/res1
                write(*,*) 'the iteration:: ', iter , 'the residual:: ',res1
                write(*,*) '==================================================================='
                call output_time
              end if
            end do
            call outPuts

        else
            write(*,*) "you've entered the wrong key word"
        end if

    end subroutine
    ! the final outputs
    subroutine outPuts
        USE head
        implicit none
        integer::i,j
        integer::mid

        mid = (imax-1)/2
        open(2,file='two_phase.dat')
        write(2,*) 'variables= "x","y","phi","P","U","V","rho","miu"'
        write(2,*) 'zone i=',imax-1 ,'j=',jmax-1 ,'F=POINT'
        do j=1,jmax-1
            do i=1,imax-1
                write(2,*) cenx(i),ceny(j),ini(i,j,1),ini(i,j,2),ini(i,j,3),ini(i,j,4),ini(i,j,5),ini(i,j,6)
            end do
        end do
        close(2)

        open(10,file='Us.dat')
        write(10,*) 'variables= "Ts" , "Us"'
        write(10,*) 'zone t="Us"'
        do i=1,ii-1
            write(10,*) Ts(i),Us(i)
        end do
        close(10)

        open(20,file='Ys.dat')
        write(20,*) 'variables= "Ts" , "Ys"'
        write(20,*) 'zone t="Ys"'
        do i=1,ii-1
            write(20,*) Ts(i),Ys(i)
        end do
        close(20)

    end subroutine outPuts

    subroutine output_time
        USE head
        USE rho_0
        implicit none
        character(len=100)::filename,directory
        directory = 'time_steps/'

        write(filename, '(A, "TimeStep_", I8.8, ".dat")') trim(directory), iter
        open(unit=file_unit, file=trim(filename))
        write(file_unit,*) 'TITLE = "Time-dependent bubble rising"'
        write(file_unit,*) 'variables= "x","y","phi"'
        write(file_unit,*) 'zone t=" time = ' , iter , ' " i=',imax-1 ,'j=',jmax-1 ,'F=POINT'
        do j=1,jmax-1
          do i=1,imax-1
            write(file_unit,*) cenx(i),ceny(j),ini(i,j,1)
          end do
        end do
        close(file_unit)


        Ts(ii) = iter*dt*sqrt(abs(gy)*dia)/dia

        tmp1=0.d0
        tmp2=0.d0

        ! to find the mass center of the bubble
        do i=1,imax-1
            do j=1,jmax-1
                if (ini(i,j,1) >= 0 .AND. ini(i,j,1) <= 0.5 )then
                    tmp1 = tmp1 + ceny(j)*vol(i,j)
                    tmp2 = tmp2 + vol(i,j)
                    if (tmp2<1e-6) then
                        tmp2 = 1e-6
                    end if
                    Ys(ii) = 1.d0/dia*(tmp1/tmp2)
                end if
            end do
        end do

        tmp1=0.d0
        tmp2=0.d0

        ! to find the velocity of mass center of the bubble
        do i=1,imax-1
            do j=1,jmax-1
                if (ini(i,j,1) >= 0 .AND. ini(i,j,1) <= 0.5 )then
                    tmp1 = tmp1 + ini(i,j,4)*vol(i,j)
                    tmp2 = tmp2 + vol(i,j)
                    if (tmp2<1e-6) then
                        tmp2 = 1e-6
                    end if
                    Us(ii) = 1.d0/(sqrt(abs(gy)*dia))*(tmp1/tmp2)
                end if
            end do
        end do

        ii = ii + 1

    end subroutine output_time

    ! the mass conversation term
    subroutine mass(M0,qs)
        USE head
        implicit none

        real(8),intent(in) ::M0
        real(8),intent(out)::qs(ntmax,ntmax)
        real(8):: mL,VL,VL2  ! higher for the liquid and lower for the gas
        real(8):: temp
        integer::i,j,in,ip,jn,jp

        VL = 0.0d0
        VL2 = 0.0d0
        temp = 0.0d0

        ! calculation of the volume of bubble (for 2-D case the area of bubble)
        do i=1,imax-1
            do j=1,jmax-1
                ! just for gas
                if (ini(i,j,1) >= 0 .AND. ini(i,j,1) <= 0.5 )then
                    VL = VL + vol(i,j)
                end if
            end do
        end do

        do i=1,imax-1
            do j=1,jmax-1
                ! bubble volume and its interface
                if (ini(i,j,1) >= 0 .AND. ini(i,j,1) <= 0.9 )then
                    VL2 = VL2 + vol(i,j)
                end if
            end do
        end do

        ! calculation of the mass of bubble gas in each iteration
        mL = ddL*VL

        ! find M*del2_PHI
        do i=1,imax-1
            do j=1,jmax-1
                ip=i-1
                in=i+1
                jp=j-1
                jn=j+1
                if (ini(i,j,1) >= 0 .AND. ini(i,j,1) <= 0.5 ) then
                    temp = temp + vol(i,j)*M*(ini(in,jn,6)+ini(in,jp,1)+ini(ip,jn,6)+ini(ip,jp,6)+4.0*ini(ip,j,6)+&
                         & 4.0*ini(in,j,6)+4.0*ini(i,jp,6)+4.0*ini(i,jn,6)-20.0*ini(i,j,6))/6.0
                end if
            end do
        end do

        do i=1,imax-1
            do j=1,jmax-1
                if(ini(i,j,1) >= 0.1 .AND. ini(i,j,1) <= 0.9) then
                    qs(i,j) = 1.d0/(VL2-VL)*((mL-M0)/((ddH-ddL)*dt) - temp)
                else
                    qs(i,j) = 0.d0
                end if
            end do
        end do


    end subroutine mass
    ! calculate the multiple relaxation time
    subroutine MRT
        USE head
        implicit none
        real(4)::a1
        integer::i,j,k

        ex(:) = (/1.0,-1.0,0.0,0.0,1.0,-1.0,-1.0,-1.0,0.0/)
        ey(:) = (/0.0,0.0,1.0,-1.0,1.0,1.0,-1.0,-1.0,0.0/)

        ! matrix M
        tm(1,:) = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
        tm(2,:) = (/-4.0,-1.0,-1.0,-1.0,-1.0,2.0,2.0,2.0,2.0/)
        tm(3,:) = (/4.0,-2.0,-2.0,-2.0,-2.0,1.0,1.0,1.0,1.0/)
        tm(4,:) = (/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
        tm(5,:) = (/0.0,-2.0,0.0,2.0,0.0,1.0,-1.0,-1.0,1.0/)
        tm(6,:) = (/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
        tm(7,:) = (/0.0,0.0,-2.0,0.0,2.0,1.0,1.0,-1.0,-1.0/)
        tm(8,:) = (/0.0,1.0,-1.0,1.0,-1.0,0.0,0.0,0.0,0.0/)
        tm(9,:) = (/0.0,0.0,0.0,0.0,0.0,1.0,-1.0,1.0,-1.0/)

        a1 = 1.0/36.0

        ! inversion of M matrix
        tminv(1,:) = (/4.0*a1,-4.0*a1,4.0*a1,0.0,0.0,0.0,0.0,0.0,0.0/)
        tminv(2,:) = (/4.0*a1,-a1,-2.0*a1,6.0*a1,-6.0*a1,0.0,0.0,9.0*a1,0.0/)
        tminv(3,:) = (/4.0*a1,-a1,-2.0*a1,0.0,0.0,6.0*a1,-6.0*a1,-9.0*a1,0.0/)
        tminv(4,:) = (/4.0*a1,-a1,-2.0*a1,-6.0*a1,6.0*a1,0.0,0.0,9.0*a1,0.0/)
        tminv(5,:) = (/4.0*a1,-a1,-2.0*a1,0.0,0.0,-6.0*a1,6.0*a1,-9.0*a1,0.0/)
        tminv(6,:) = (/4.0*a1,2.0*a1,a1,6.0*a1,3.0*a1,6.0*a1,3.0*a1,0.0,9.0*a1/)
        tminv(7,:) = (/4.0*a1,2.0*a1,a1,-6.0*a1,-3.0*a1,6.0*a1,3.0*a1,0.0,-9.0*a1/)
        tminv(8,:) = (/4.0*a1,2.0*a1,a1,-6.0*a1,-3.0*a1,-6.0*a1,-3.0*a1,0.0,9.0*a1/)
        tminv(9,:) = (/4.0*a1,2.0*a1,a1,6.0*a1,3.0*a1,-6.0*a1,-3.0*a1,0.0,-9.0*a1/)

        s1 = 1.0
        s2 = 1.4
        s3 = 1.4
        s4 = 1.0
        s5 = 1.2
        s6 = 1.0
        s7 = 1.2
        s8 = 1.0/(niu/(1.0*Cs**2)+0.5)
        s9 = s8

        sm(1) = s1
        sm(2) = s2
        sm(3) = s3
        sm(4) = s4
        sm(5) = s5
        sm(6) = s6
        sm(7) = s7
        sm(8) = s8
        sm(9) = s9

        ! calculate the updated M-1*S
        do i=1,9
            do j=1,9
                stmiv(i,j) = tminv(i,j)*sm(j)
            end do
        end do

    end subroutine
    ! calculate the equilibrium function using the macro on the interfaces
    subroutine feq(dd,pp,uu,vv,ff)
        USE head
        implicit none

        integer:: n
        real(8):: dd,uu,vv,pp,ff(9)
        real(8):: a1,b1,c1,uu2,vv2,uv,u2v2

        cs = 1.0d0/sqrt(3.0d0)
        a1 = 1.d0/9.d0
        b1 = 1.d0/36.d0
        c1 = a1*4.d0

        ff(1) = a1*(pp + dd*(uu+1.5d0*uu*uu-0.5d0*(uu*uu+vv*vv)))
        ff(2) = a1*(pp + dd*(-uu +1.5d0*uu*uu-0.5d0*(uu*uu+vv*vv)))
        ff(3) = a1*(pp + dd*(vv +1.5d0*vv*vv-0.5d0*(uu*uu+vv*vv)))
        ff(4) = a1*(pp + dd*(-vv +1.5d0*vv*vv-0.5d0*(uu*uu+vv*vv)))
        ff(5) = b1*(pp + dd*((uu+vv) +1.5d0*(uu+vv)**2-0.5d0*(uu*uu+vv*vv)))
        ff(6) = b1*(pp + dd*((-uu+vv) +1.5d0*(-uu+vv)**2-0.5d0*(uu*uu+vv*vv)))
        ff(7) = b1*(pp + dd*((-uu-vv) +1.5d0*(-uu-vv)**2-0.5d0*(uu*uu+vv*vv)))
        ff(8) = b1*(pp + dd*((uu-vv) +1.5d0*(uu-vv)**2-0.5d0*(uu*uu+vv*vv)))
        ff(9) = c1*(pp - dd*0.5d0*(uu*uu + vv*vv))

        return
    end subroutine
    ! calculate the equilibrium g function using the macro on the interfaces
    subroutine geq(phi,miu,uu,vv,gg)

        implicit none

        integer::n
        real(8)::phi,miu,uu,vv,gg(9)
        real(8)::a1,b1,c1,u2,uu2,vv2,etha

        etha = 10.0d0
        a1 = 1.d0/9.d0
        b1 = 1.d0/36.d0
        c1 = a1*4.d0

        gg(9) = phi + 3.0d0*miu*etha*(1.d0 - c1)
        gg(1) = 3.d0*a1*(miu*etha + phi*uu)
        gg(2) = 3.d0*a1*(miu*etha + phi*(-uu))
        gg(3) = 3.d0*a1*(miu*etha + phi*vv)
        gg(4) = 3.d0*a1*(miu*etha + phi*(-vv))
        gg(5) = 3.d0*b1*(miu*etha + phi*(uu + vv))
        gg(6) = 3.d0*b1*(miu*etha + phi*(-uu + vv))
        gg(7) = 3.d0*b1*(miu*etha + phi*(-uu - vv))
        gg(8) = 3.d0*b1*(miu*etha + phi*(uu - vv))

    end subroutine
    ! calculating the reformulated equilibrium of f on the faces use rho(alpha = 0)
    subroutine feqmi(n,dd,pp,uu,vv,ffi)

        implicit none

        integer:: n
        real(8):: dd,uu,vv,pp,ffi
        real(8):: a1,b1,c1,uu2,vv2,uv,u2v2
        real(8):: Cs

        cs = 1.0d0/sqrt(3.0d0)
        a1 = 1.d0/9.d0
        b1 = 1.d0/36.d0
        c1 = a1*4.d0

        select case (n)
            case (1)
                ffi = a1*(pp + dd*(uu+1.5d0*uu*uu-0.5d0*(uu*uu+vv*vv)))
            case (2)
                ffi = a1*(pp + dd*(-uu +1.5d0*uu*uu-0.5d0*(uu*uu+vv*vv)))
            case (3)
                ffi = a1*(pp + dd*(vv +1.5d0*vv*vv-0.5d0*(uu*uu+vv*vv)))
            case (4)
                ffi = a1*(pp + dd*(-vv +1.5d0*vv*vv-0.5d0*(uu*uu+vv*vv)))
            case (5)
                ffi = b1*(pp + dd*((uu+vv) +1.5d0*(uu+vv)**2-0.5d0*(uu*uu+vv*vv)))
            case (6)
                ffi = b1*(pp + dd*((-uu+vv) +1.5d0*(-uu+vv)**2-0.5d0*(uu*uu+vv*vv)))
            case (7)
                ffi = b1*(pp + dd*((-uu-vv) +1.5d0*(-uu-vv)**2-0.5d0*(uu*uu+vv*vv)))
            case (8)
                ffi = b1*(pp + dd*((uu-vv) +1.5d0*(uu-vv)**2-0.5d0*(uu*uu+vv*vv)))
            case (9)
                ffi = c1*(pp - dd*0.5d0*(uu*uu + vv*vv))
        end select

        return
    end subroutine feqmi
    ! calculating the equilibrium of g on the faces
    subroutine geqi(n,phi,miu,uu,vv,ggi)

        implicit none

        integer::n
        real(8)::phi,miu,uu,vv,ggi
        real(8)::a1,b1,c1,u2,uu2,vv2,etha

        etha = 10.0d0
        a1 = 1.d0/9.d0
        b1 = 1.d0/36.d0
        c1 = a1*4.d0

        u2  = uu*uu+vv*vv
        uu2 = uu*uu2
        vv2 = vv*vv

        select case (n)
            case (9)
                ggi = phi + 3.0d0*miu*etha*(c1 - 1.d0)
            case (1)
                ggi = 3.d0*a1*(miu*etha + phi*uu)
            case (2)
                ggi = 3.d0*a1*(miu*etha + phi*(-uu))
            case (3)
                ggi = 3.d0*a1*(miu*etha + phi*vv)
            case (4)
                ggi = 3.d0*a1*(miu*etha + phi*(-vv))
            case (5)
                ggi = 3.d0*b1*(miu*etha + phi*(uu + vv))
            case (6)
                ggi = 3.d0*b1*(miu*etha + phi*(-uu + vv))
            case (7)
                ggi = 3.d0*b1*(miu*etha + phi*(-uu - vv))
            case (8)
                ggi = 3.d0*b1*(miu*etha + phi*(uu - vv))
        end select

    end subroutine geqi
    ! update relaxation matrix
    subroutine relaxation(alpha,zeta)
        USE head
        implicit none
        ! alpha => niu , zeta => faceX => deltaT
        real(8),INTENT(IN)::alpha,zeta
        integer::i,j

        sm(8) = 1.0/(alpha/(beta*Cs**2))
        sm(9) = sm(8)

        ! calculate the updated M-1*S
        do i=1,9
            do j=1,9
                stmiv(i,j) = tminv(i,j)*sm(j)
            end do
        end do

    end subroutine
    ! calculating the fluxes variables using interpolation
    subroutine fluxff
!=====================================================================
      USE head
      implicit none

      integer:: i,j,k,nn,kk
      real(8)::  dd,pp,miu,uu,vv,phi,taof,taog,ss0,f1,f2,f3,f4
      real(8)::  ldx1,ldy1,lux1,luy1,lvx1,lvy1,lwx1,lwy1
      real(8)::  ldx2,ldy2,lux2,luy2,lvx2,lvy2,lwx2,lwy2
      real(8)::  lpx1,lpy1,lmx1,lmy1,lpx2,lpy2,lmx2,lmy2
      real(8)::  lphix1,lphiy1,lphix2,lphiy2
      real(8)::  lmiux1,lmiuy1,lmiux2,lmiuy2
      real(8)::  tempd,tempu,tempv,temd,temu,temv,tempmiu,tempp,tempphi
      real(8)::  dx,dx1,dx2,dy,dy1,dy2
      real(8)::  pp1,uu1,vv1,pp2,uu2,vv2,miu1,miu2,phi1,phi2,dd1,dd2
      real(8)::  pp3,uu3,vv3,pp4,uu4,vv4,miu3,miu4,phi3,phi4,dd3,dd4
      real(8)::  ddf1,ddf2,ddf3,ddf4,ddf
      real(8)::  feq0(9),geq0(9),ffi(9),ggi(9),fs(9),gs(9)
      real(8)::  tmp1,tmp2,tmp3,tmp4m,miufx(9),miufy(9)

      do j = 0,jmax
          do i = 0,imax
            flux0(i,j,1) = 0.d0
            flux0(i,j,2) = 0.d0
            flux0(i,j,3) = 0.d0
            flux0(i,j,4) = 0.d0
          end do
      end do

      !flux in the x direction
      do j = 1,jmax-1
          ss0 = leny(j)
        do i = 1,imax
              dy = facey(i,j,1)
              dx = facex(i,j,1)
              dx1= facex(i,j,2)
              dx2= facex(i,j,3)

              phi1= ini(i-1,j,1)
              miu1= ini(i-1,j,6)
              pp1 = ini(i-1,j,2)
              uu1 = ini(i-1,j,3)
              vv1 = ini(i-1,j,4)

              lphix1= limx(i-1,j,1)
              lphiy1= limy(i-1,j,1)*dy
              lmiux1= limx(i-1,j,6)
              lmiuy1= limy(i-1,j,6)*dy
              lpx1  = limx(i-1,j,2)
              lpy1  = limy(i-1,j,2)*dy
              lux1  = limx(i-1,j,3)
              luy1  = limy(i-1,j,3)*dy
              lvx1  = limx(i-1,j,4)
              lvy1  = limy(i-1,j,4)*dy

              phi2= ini(i,j,1)
              miu2= ini(i,j,6)
              pp2 = ini(i,j,2)
              uu2 = ini(i,j,3)
              vv2 = ini(i,j,4)

              lphix2= limx(i,j,1)
              lphiy2= limy(i,j,1)*dy
              lmiux2= limx(i,j,6)
              lmiuy2= limy(i,j,6)*dy
              lpx2  = limx(i,j,2)
              lpy2  = limy(i,j,2)*dy
              lux2  = limx(i,j,3)
              luy2  = limy(i,j,3)*dy
              lvx2  = limx(i,j,4)
              lvy2  = limy(i,j,4)*dy

              !left three points
              tempp  = lpx1*dx1
              tempphi= lphix1*dx1
              tempmiu= lmiux1*dx1
              tempu  = lux1*dx1
              tempv  = lvx1*dx1

              nn = 1
              pp3  = pp1 + tempp
              phi3 = phi1+ tempphi
              if(phi3>1.d0) phi3 = 1.d0
              if(phi3<0.d0) phi3 = 0.d0
              miu3 = miu1+ tempmiu
              uu3 = uu1 + tempu
              vv3 = vv1 + tempv
              dd3 = rho_facexO(i,j)
              ddf3= ddH*phi3+(1.d0-phi3)*ddL
              call feqmi(nn,dd3,pp3,uu3,vv3,ffi(nn))
              call geqi(nn,phi3,miu3,uu3,vv3,ggi(nn))
              miufx(nn) = miu3

              nn = 5
              pp = pp3  - lpy1
              miu= miu3 - lmiuy1
              phi= phi3 - lphiy1
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              uu = uu3  - luy1
              vv = vv3  - lvy1
              dd = rho_facexO(i,j)
              ddf= ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufx(nn) = miu

              nn = 8
              pp = pp3  + lpy1
              miu= miu3 + lmiuy1
              phi= phi3 + lphiy1
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              uu = uu3  + luy1
              vv = vv3  + lvy1
              dd = rho_facexO(i,j)
              ddf= ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufx(nn) = miu

              !right three points
              tempp  = lpx2*dx2
              tempphi= lphix2*dx2
              tempmiu= lmiux2*dx2
              tempu  = lux2*dx2
              tempv  = lvx2*dx2

              nn = 2
              pp3  = pp2 + tempp
              phi3 = phi2+ tempphi
              if(phi3>1.d0) phi3 = 1.d0
              if(phi3<0.d0) phi3 = 0.d0
              miu3 = miu2+ tempmiu
              uu3 = uu2 + tempu
              vv3 = vv2 + tempv
              dd3 = rho_facexO(i,j)
              ddf3= ddH*phi3+(1.d0-phi3)*ddL
              call feqmi(nn,dd3,pp3,uu3,vv3,ffi(nn))
              call geqi(nn,phi3,miu3,uu3,vv3,ggi(nn))
              miufx(nn) = miu3

              nn = 7
              pp = pp3  + lpy2
              miu= miu3 + lmiuy2
              phi= phi3 + lphiy2
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              uu = uu3  + luy2
              vv = vv3  + lvy2
              dd = rho_facexO(i,j)
              ddf= ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufx(nn) = miu

              nn = 6
              pp = pp3  - lpy2
              miu= miu3 - lmiuy2
              phi= phi3 - lphiy2
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              uu = uu3  - luy2
              vv = vv3  - lvy2
              dd = rho_facexO(i,j)
              ddf= ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufx(nn) = miu

              !middle three points
              dx1 =  lenx1(i-1)
              dx2 = -lenx1(i)

              uu3 = uu1 + lux1*dx1
              vv3 = vv1 + lvx1*dx1
              pp3 = pp1 + lpx1*dx1
              phi3= phi1+ lphix1*dx1
              if(phi3>1.d0) phi3 = 1.d0
              if(phi3<0.d0) phi3 = 0.d0
              miu3= miu1+ lmiux1*dx1

              uu4 = uu2 + lux2*dx2
              vv4 = vv2 + lvx2*dx2
              pp4 = pp2 + lpx2*dx2
              phi4= phi2+ lphix2*dx2
              if(phi4>1.d0) phi4 = 1.d0
              if(phi4<0.d0) phi4 = 0.d0
              miu4= miu2+ lmiux2*dx2

              nn = 9
              pp3 = (pp3 + pp4)*0.5d0
              miu3= (miu3+miu4)*0.5d0
              phi3= (phi3+phi4)*0.5d0
              if(phi3>1.d0) phi3 = 1.d0
              if(phi3<0.d0) phi3 = 0.d0
              uu3 = (uu3 + uu4)*0.5d0
              vv3 = (vv3 + vv4)*0.5d0
              dd3 = ddH*phi3+(1.d0-phi3)*ddL ! rho (aplha = 9)
              rho_facex(i,j) = dd3
              ddf3= ddH*phi3+(1.d0-phi3)*ddL ! rho (aplha = 9)
              call feqmi(nn,rho_facexO(i,j),pp3,uu3,vv3,ffi(nn))
              call geqi(nn,phi3,miu3,uu3,vv3,ggi(nn))
              miufx(nn) = miu3

              lpy1 = 0.5d0*(lpy1 + lpy2)
              lphiy1= 0.5d0*(lphiy1+lphiy2)
              lmiuy1= 0.5d0*(lmiuy1+lmiuy2)
              luy1 = 0.5d0*(luy1 + luy2)
              lvy1 = 0.5d0*(lvy1 + lvy2)

              nn  = 3
              pp  = pp3 - lpy1
              phi = phi3 - lphiy1
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              miu = miu3 - lmiuy1
              uu  = uu3 - luy1
              vv  = vv3 - lvy1
              dd  = rho_facexO(i,j)
              ddf = ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufx(nn) = miu

              nn  = 4
              pp  = pp3 + lpy1
              phi = phi3 + lphiy1
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              miu = miu3 + lmiuy1
              uu  = uu3 + luy1
              vv  = vv3 + lvy1
              dd  = rho_facexO(i,j)
              ddf = ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufx(nn) = miu

            ! calculate the macros on the faces :
            != = = = = = = = = = = = = = = = = = =
              pp = ffi(1)+ffi(2)+ffi(3)+ffi(4)+ffi(5)+ffi(6)+ffi(7)+ffi(8)+ffi(9)
              phi= ggi(1)+ggi(2)+ggi(3)+ggi(4)+ggi(5)+ggi(6)+ggi(7)+ggi(8)+ggi(9)
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              dd = ddL*(1.d0-phi)+phi*ddH
              uu = 1.d0/(dd*Cs**2)*(ffi(1)-ffi(2)+ffi(5)-ffi(6)-ffi(7)+ffi(8))
              vv = 1.d0/(dd*Cs**2)*(ffi(3)-ffi(4)+ffi(5)+ffi(6)-ffi(7)-ffi(8))
              miu = wa(1)*miufx(1)+wa(2)*miufx(2)+wa(3)*miufx(3)+wa(4)*miufx(4)+wa(5)*miufx(5)&
              &    +wa(6)*miufx(6)+wa(7)*miufx(7)+wa(8)*miufx(8)+wa(9)*miufx(9)
              if (abs(miu) < 1e-10 ) miu = 0.d0

              ! calculate the equilibrium functions on the center of interfaces:
              call feq(dd,pp,uu,vv,feq0)
              call geq(phi,miu,uu,vv,geq0)
            != = = = = = = = = = = = = = = = = = =
            ! for the Newtonian
            !==================
              taof = (taoxH(i,j)*taoxL(i,j))/(phi*taoxL(i,j)+(1.d0-phi)*taoxH(i,j))
              taog = taogx(i,j)

              do k=1,9
                fs(k) = feq0(k)-(taof)*(feq0(k)-ffi(k))
                gs(k) = geq0(k)-(taog)*(geq0(k)-ggi(k))
              end do

              f1 = (gs(1)-gs(2)+gs(5)-gs(6)-gs(7)+gs(8))*ss0
              f2 = (feq0(1)-feq0(2)+feq0(5)-feq0(6)-feq0(7)+feq0(8))*ss0
              f3 = (fs(1)+fs(2)+fs(5)+fs(6)+fs(7)+fs(8))*ss0
              f4 = (fs(5)-fs(6)+fs(7)-fs(8))*ss0
            !==================

              flux0(i-1,j,1) = flux0(i-1,j,1)+f1
              flux0(i-1,j,2) = flux0(i-1,j,2)+f2
              flux0(i-1,j,3) = flux0(i-1,j,3)+f3
              flux0(i-1,j,4) = flux0(i-1,j,4)+f4

              flux0(i,j,1) = flux0(i,j,1)-f1
              flux0(i,j,2) = flux0(i,j,2)-f2
              flux0(i,j,3) = flux0(i,j,3)-f3
              flux0(i,j,4) = flux0(i,j,4)-f4
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

              phi1= ini(i,j-1,1)
              miu1= ini(i,j-1,6)
              pp1 = ini(i,j-1,2)
              uu1 = ini(i,j-1,3)
              vv1 = ini(i,j-1,4)

              lphix1= limx(i,j-1,1)*dx
              lphiy1= limy(i,j-1,1)
              lmiux1= limx(i,j-1,6)*dx
              lmiuy1= limy(i,j-1,6)
              lpx1  = limx(i,j-1,2)*dx
              lpy1  = limy(i,j-1,2)
              lux1  = limx(i,j-1,3)*dx
              luy1  = limy(i,j-1,3)
              lvx1  = limx(i,j-1,4)*dx
              lvy1  = limy(i,j-1,4)

              phi2= ini(i,j,1)
              miu2= ini(i,j,6)
              pp2 = ini(i,j,2)
              uu2 = ini(i,j,3)
              vv2 = ini(i,j,4)

              lphix2= limx(i,j,1)*dx
              lphiy2= limy(i,j,1)
              lmiux2= limx(i,j,6)*dx
              lmiuy2= limy(i,j,6)
              lpx2  = limx(i,j,2)*dx
              lpy2  = limy(i,j,2)
              lux2  = limx(i,j,3)*dx
              luy2  = limy(i,j,3)
              lvx2  = limx(i,j,4)*dx
              lvy2  = limy(i,j,4)

              !bottom three points
              tempp  = lpy1*dy1
              tempphi= lphiy1*dy1
              tempmiu= lmiuy1*dy1
              tempu  = luy1*dy1
              tempv  = lvy1*dy1

              nn = 3
              pp3  = pp1 + tempp
              phi3 = phi1+ tempphi
              if(phi3>1.d0) phi3 = 1.d0
              if(phi3<0.d0) phi3 = 0.d0
              miu3 = miu1+ tempmiu
              uu3 = uu1 + tempu
              vv3 = vv1 + tempv
              dd3 = rho_faceyO(i,j)
              ddf3= ddH*phi3+(1.d0-phi3)*ddL
              call feqmi(nn,dd3,pp3,uu3,vv3,ffi(nn))
              call geqi(nn,phi3,miu3,uu3,vv3,ggi(nn))
              miufy(nn) = miu3

              nn = 5
              pp = pp3  - lpx1
              miu= miu3 - lmiux1
              phi= phi3 - lphix1
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              uu = uu3  - lux1
              vv = vv3  - lvx1
              dd = rho_faceyO(i,j)
              ddf = ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufy(nn) = miu

              nn = 6
              pp = pp3  + lpx1
              miu= miu3 + lmiux1
              phi= phi3 + lphix1
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              uu = uu3  + lux1
              vv = vv3  + lvx1
              dd = rho_faceyO(i,j)
              ddf = ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufy(nn) = miu

              !top three points
              tempp  = lpy2*dy2
              tempphi= lphiy2*dy2
              tempmiu= lmiuy2*dy2
              tempu  = luy2*dy2
              tempv  = lvy2*dy2

              nn = 4
              pp3  = pp2 + tempp
              phi3 = phi2+ tempphi
              if(phi3>1.d0) phi3 = 1.d0
              if(phi3<0.d0) phi3 = 0.d0
              miu3 = miu2+ tempmiu
              uu3 = uu2 + tempu
              vv3 = vv2 + tempv
              dd3 = rho_faceyO(i,j)
              ddf3 = ddH*phi3+(1.d0-phi3)*ddL
              call feqmi(nn,dd3,pp3,uu3,vv3,ffi(nn))
              call geqi(nn,phi3,miu3,uu3,vv3,ggi(nn))
              miufy(nn) = miu3

              nn = 7
              pp = pp3  + lpx2
              miu= miu3 + lmiux2
              phi= phi3 + lphix2
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              uu = uu3  + lux2
              vv = vv3  + lvx2
              dd = rho_faceyO(i,j)
              ddf = ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufy(nn) = miu

              nn = 8
              pp = pp3  - lpx2
              miu= miu3 - lmiux2
              phi= phi3 - lphix2
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              uu = uu3  - lux2
              vv = vv3  - lvx2
              dd = rho_faceyO(i,j)
              ddf = ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufy(nn) = miu

              !middle three points
              dy1 =  leny1(j-1)
              dy2 = -leny1(j)

              uu3 = uu1 + luy1*dy1
              vv3 = vv1 + lvy1*dy1
              pp3 = pp1 + lpy1*dy1
              phi3= phi1+ lphiy1*dy1
              if(phi3>1.d0) phi3 = 1.d0
              if(phi3<0.d0) phi3 = 0.d0
              miu3= miu1+ lmiuy1*dy1

              uu4 = uu2 + luy2*dy2
              vv4 = vv2 + lvy2*dy2
              pp4 = pp2 + lpy2*dy2
              phi4= phi2+ lphiy2*dy2
              if(phi4>1.d0) phi4 = 1.d0
              if(phi4<0.d0) phi4 = 0.d0
              miu4= miu2+ lmiuy2*dy2

              nn = 9
              pp3 = (pp3 + pp4)*0.5d0
              miu3= (miu3+miu4)*0.5d0
              phi3= (phi3+phi4)*0.5d0
              if(phi3>1.d0) phi3 = 1.d0
              if(phi3<0.d0) phi3 = 0.d0
              uu3 = (uu3 + uu4)*0.5d0
              vv3 = (vv3 + vv4)*0.5d0
              dd3 = ddH*phi3+(1.d0-phi3)*ddL
              rho_facey(i,j) = dd3
              ddf3= ddH*phi3+(1.d0-phi3)*ddL
              call feqmi(nn,rho_faceyO(i,j),pp3,uu3,vv3,ffi(nn))
              call geqi(nn,phi3,miu3,uu3,vv3,ggi(nn))
              miufy(nn) = miu3

              lpx1 = 0.5d0*(lpx1 + lpx2)
              lphix1= 0.5d0*(lphix1+lphix2)
              lmiux1= 0.5d0*(lmiux1+lmiux2)
              lux1 = 0.5d0*(lux1 + lux2)
              lvx1 = 0.5d0*(lvx1 + lvx2)

              nn  = 1
              pp  = pp3 - lpx1
              phi = phi3 - lphix1
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              miu = miu3 - lmiux1
              uu  = uu3 - lux1
              vv  = vv3 - lvx1
              dd  = rho_faceyO(i,j)
              ddf = ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufy(nn) = miu

              nn  = 2
              pp  = pp3 + lpx1
              phi = phi3 + lphix1
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              miu = miu3 + lmiux1
              uu  = uu3 + lux1
              vv  = vv3 + lvx1
              dd  = rho_faceyO(i,j)
              ddf = ddH*phi+(1.d0-phi)*ddL
              call feqmi(nn,dd,pp,uu,vv,ffi(nn))
              call geqi(nn,phi,miu,uu,vv,ggi(nn))
              miufy(nn) = miu

            ! calculate the macros on the faces :
            != = = = = = = = = = = = = = = = = = =
              pp = ffi(1)+ffi(2)+ffi(3)+ffi(4)+ffi(5)+ffi(6)+ffi(7)+ffi(8)+ffi(9)
              phi= ggi(1)+ggi(2)+ggi(3)+ggi(4)+ggi(5)+ggi(6)+ggi(7)+ggi(8)+ggi(9)
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              dd = ddL*(1.d0-phi)+phi*ddH
              uu = 1.d0/(dd*Cs**2)*(ffi(1)-ffi(2)+ffi(5)-ffi(6)-ffi(7)+ffi(8))
              vv = 1.d0/(dd*Cs**2)*(ffi(3)-ffi(4)+ffi(5)+ffi(6)-ffi(7)-ffi(8))
              miu = wa(1)*miufy(1)+wa(2)*miufy(2)+wa(3)*miufy(3)+wa(4)*miufy(4)+wa(5)*miufy(5)&
              &    +wa(6)*miufy(6)+wa(7)*miufy(7)+wa(8)*miufy(8)+wa(9)*miufy(9)
              if (abs(miu) < 1e-10 ) miu = 0.d0

              call feq(dd,pp,uu,vv,feq0)
              call geq(phi,miu,uu,vv,geq0)
            != = = = = = = = = = = = = = = = = = =
            ! for the Newtonian
            !==================
              taof = (taoyH(i,j)*taoyL(i,j))/(phi*taoyL(i,j)+(1.d0-phi)*taoyH(i,j))
              taog = taogy(i,j)

              do k=1,9
                fs(k) = feq0(k)-(taof)*(feq0(k)-ffi(k))
                gs(k) = geq0(k)-(taog)*(geq0(k)-ggi(k))
              end do

              f1 = (gs(3)-gs(4)+gs(5)+gs(6)-gs(7)-gs(8))*ss0
              f2 = (feq0(3)-feq0(4)+feq0(5)+feq0(6)-feq0(7)-feq0(8))*ss0
              f3 = (fs(5)-fs(6)+fs(7)-fs(8))*ss0
              f4 = (fs(3)+fs(4)+fs(5)+fs(6)+fs(7)+fs(8))*ss0
            !==================

              flux0(i,j-1,1) = flux0(i,j-1,1)+f1
              flux0(i,j-1,2) = flux0(i,j-1,2)+f2
              flux0(i,j-1,3) = flux0(i,j-1,3)+f3
              flux0(i,j-1,4) = flux0(i,j-1,4)+f4

              flux0(i,j,1) = flux0(i,j,1)-f1
              flux0(i,j,2) = flux0(i,j,2)-f2
              flux0(i,j,3) = flux0(i,j,3)-f3
              flux0(i,j,4) = flux0(i,j,4)-f4

        end do
      end do

      return
    end subroutine fluxff

    subroutine tensionForces
    ! calculating the gradients 2nd order
      USE head
      implicit none
      integer::i,j,ip,in,jp,jn

      ! F = -phi * grad miu
      do i=1,imax-1
        do j=1,jmax-1
            ip=i-1
            in=i+1
            jp=j-1
            jn=j+1
            !Fsx(i,j) = -ini(i,j,1)*((ini(in,j,6)-ini(ip,j,6))/3.0 &
            !&+(ini(in,jp,6)-ini(ip,jp,6))/12.0 +(ini(in,jn,6)-ini(ip,jn,6))/12.0)
            !Fsy(i,j) = -ini(i,j,1)*((ini(i,jn,6)-ini(i,jp,6))/3.0 &
            !&+(ini(ip,jn,6)-ini(ip,jp,6))/12.0 +(ini(in,jn,6)-ini(in,jp,6))/12.0)
            Fsx(i,j) = -ini(i,j,1)*limx(i,j,6)
            Fsy(i,j) = -ini(i,j,1)*limy(i,j,6)
        end do
      end do
    ! handling the boundary nodes :
        i = 0
        do j=1,jmax-1
            Fsx(i,j) = Fsx(i+1,j)
        end do

        i = imax
        do j=1,jmax-1
            Fsx(i,j) = Fsx(i-1,j)
        end do

        j=0
        do i=1,imax-1
            Fsy(i,j) = Fsy(i,j+1)
        end do

        j=jmax
        do i=1,imax-1
            Fsy(i,j) = Fsy(i,j-1)
        end do

  end subroutine tensionForces

    subroutine gradients
    USE head
    implicit none
    integer:: i,j,k

      !derivatives in the x and y directions 2nd order
    do j = 1,jmax-1
        do i = 1,imax-1
              limx(i,j,:) = cx1(i)*(ini(i+1,j,:)-ini(i,j,:))+cx2(i)*(ini(i,j,:)-ini(i-1,j,:))
              limy(i,j,:) = cy1(j)*(ini(i,j+1,:)-ini(i,j,:))+cy2(j)*(ini(i,j,:)-ini(i,j-1,:))
              do k=1,6
                  if (abs(limx(i,j,k)) < 1e-10) limx(i,j,k) = 0.d0
                  if (abs(limy(i,j,k)) < 1e-10) limy(i,j,k) = 0.d0
              end do
        end do
    end do

    !flow properties in ghost cells
    !boundary conditions for gradients in ghost cells
      i = 0
      do j = 1,jmax-1

          limx(i,j,1) =    (ini(i+1,j,1)-ini(i,j,1))/lenx(i)
          limx(i,j,2) =    (ini(i+1,j,2)-ini(i,j,2))/lenx(i)
          limx(i,j,3) =    (ini(i+1,j,3)-ini(i,j,3))/lenx(i)
          limx(i,j,4) =    (ini(i+1,j,4)-ini(i,j,4))/lenx(i)
          limx(i,j,5) =    (ini(i+1,j,5)-ini(i,j,5))/lenx(i)
          limx(i,j,6) =    (ini(i+1,j,6)-ini(i,j,6))/lenx(i)

      end do

      i = imax
      do j = 1,jmax-1

          limx(i,j,1) =    (ini(i,j,1)-ini(i-1,j,1))/lenx(i)
          limx(i,j,2) =    (ini(i,j,2)-ini(i-1,j,2))/lenx(i)
          limx(i,j,3) =    (ini(i,j,3)-ini(i-1,j,3))/lenx(i)
          limx(i,j,4) =    (ini(i,j,4)-ini(i-1,j,4))/lenx(i)
          limx(i,j,5) =    (ini(i,j,5)-ini(i-1,j,5))/lenx(i)
          limx(i,j,6) =    (ini(i,j,6)-ini(i-1,j,6))/lenx(i)

      end do

      j = 0
      do i = 1,imax-1

          limy(i,j,1) =    (ini(i,j+1,1)-ini(i,j,1))/leny(j)
          limy(i,j,2) =    (ini(i,j+1,2)-ini(i,j,2))/leny(j)
          limy(i,j,3) =    (ini(i,j+1,3)-ini(i,j,3))/leny(j)
          limy(i,j,4) =    (ini(i,j+1,4)-ini(i,j,4))/leny(j)
          limy(i,j,5) =    (ini(i,j+1,5)-ini(i,j,5))/leny(j)
          limy(i,j,6) =    (ini(i,j+1,6)-ini(i,j,6))/leny(j)

      end do

      j = jmax
      do i = 1,imax-1

          limy(i,j,1) =   (ini(i,j,1)-ini(i,j-1,1))/leny(j)
          limy(i,j,2) =   (ini(i,j,2)-ini(i,j-1,2))/leny(j)
          limy(i,j,3) =   (ini(i,j,3)-ini(i,j-1,3))/leny(j)
          limy(i,j,4) =   (ini(i,j,3)-ini(i,j-1,3))/leny(j)
          limy(i,j,5) =   (ini(i,j,3)-ini(i,j-1,3))/leny(j)
          limy(i,j,6) =   (ini(i,j,6)-ini(i,j-1,6))/leny(j)

      end do

      return
    end subroutine gradients
    ! calculating Laplacian for the miu potential chemical function
    subroutine Laplacian
        ! updating the relaxation time
        USE head
        implicit none
        integer :: i,j,ip,in,jn,jp
        real(8) :: h(ntmax),k(ntmax)

        do i=0,imax
            h(i) = x(i+1)-x(i)
        end do

        do j=0,jmax
            k(j) = y(j+1)-y(j)
        end do

        do i=1,imax-1
            ip = i+1
            in = i-1
            do j=1,jmax-1
                jp = j+1
                jn = j-1
                limx2(i,j) = (2.d0/(h(i)+h(i-1)))*((ini(ip,j,6)-ini(i,j,6))/h(i) - (ini(i,j,6)-ini(in,j,6))/h(i-1))
                limy2(i,j) = (2.d0/(k(j)+k(j-1)))*((ini(i,jp,6)-ini(i,j,6))/k(j)-(ini(i,j,6)-ini(i,jn,6))/k(i-1))
            end do
        end do

        ! handling Boundary Conditions for laplacian operator :
        i=0
        do j=1,jmax-1
            limx2(i,j) = (ini(i+2,j,6)-2.d0*ini(i+1,j,6)+ini(i,j,6))/lenx(i)
        end do

        i=imax
        do j=1,jmax-1
            limx2(i,j) = (ini(i,j,6)-2.d0*ini(i-1,j,6)+ini(i-2,j,6))/lenx(i)
        end do

        j=0
        do i=i,imax-1
            limy2(i,j) = (ini(i,j+2,6)-2.d0*ini(i,j+1,6)+ini(i,j,6))/leny(j)
        end do

        j=jmax
        do i=1,imax-1
            limy2(i,j) = (ini(i,j,6)-2.d0*ini(i,j-1,6)+ini(i,j-2,6))/leny(j)
        end do


        return
    end subroutine Laplacian

    subroutine updateTao
      USE head
      implicit none
      integer::i,j
      real(8)::dx,dy,slen,dx1,dy1

      do j = 1,jmax-1
          do i = 1, imax
              taoxH(i,j) = niuH/(Cs**2*facex(i,j,1))+0.5d0
              taoxL(i,j) = niuL/(Cs**2*facex(i,j,1))+0.5d0
              taogx(i,j) = 1.2d0!M/(etha*slen)
          end do
      end do

      do i = 1,imax-1
          do j = 1,jmax
              taoyH(i,j) = niuH/(Cs**2*facey(i,j,1))+0.5d0
              taoyL(i,j) = niuL/(Cs**2*facey(i,j,1))+0.5d0
              taogy(i,j) = 1.2d0!M/(etha*slen)
          end do
       end do

    end subroutine updateTao

    subroutine fluxMesh
      USE head
      implicit none

      integer::i,j
      real(8)::dx,dy,slen,dx1,dy1

      do i = 0,imax
          lenx(i)=x(i+1) - x(i)
      end do

      do j = 0,jmax
          leny(j)=y(j+1) - y(j)
      end do

      do i = 0,imax
          lenx1(i) = 0.5d0*lenx(i)
      end do

      do i = 0,jmax
          leny1(i) = 0.5d0*leny(i)
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
          end do
      end do

      do i = 1,imax-1
          dx= lenx(i)
          do j = 1,jmax
              dy   = leny(j-1)
              dy1  = leny(j)
              slen = 0.25d0*min(dx,dy,dy1)! dy=dt for lattices
              facey(i,j,1)= slen
              facey(i,j,2)= 0.5d0*dy - slen
              facey(i,j,3)= slen - 0.5d0*dy1
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

      open(112,file='mesh.dat')
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
      real(8) ::rp,eta,rlg_x,rlg_y,disr,tmp

        ! Specify dimensions
        imax = 161   ! Width (number of points in x direction)
        jmax = 321   ! Height (number of points in y direction)

        pi  = 4.d0 * atan(1.d0)
        eta = 500000.0d0

        ! Two separate coefficients for x and y directions
        rlg_x = 80.0d0   ! Coefficient for x-direction
        rlg_y = 160.0d0   ! Coefficient for y-direction

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

    subroutine BoundaryCondition(ww)
      USE head
      implicit none

      real(8),intent(in):: ww(0:ntmax,0:ntmax,4)
      integer:: i,j,in,ip,jp,jn
      real(8):: dd,phi,pp,miu

      ! ini(1)=>phi , ini(2)=>p , ini(3)=>u , ini(4)=>v ,ini(5)=>rho, ini(6)=>miu
      do j = 1,jmax-1
          do i = 1,imax-1
              phi = ww(i,j,1)
              if (phi < 0.0d0 ) phi = 0.0d0
              if (phi > 1.0d0 ) phi = 1.0d0
              ini(i,j,1) = phi !order parameter
              dd = ddL*(1.0d0-phi)+phi*ddH !rho density
              ini(i,j,5) = dd ! store the density for gradient
              pp = ww(i,j,2)
              ini(i,j,2) = pp !pressure
              ini(i,j,3) = (3.0/dd)*ww(i,j,3) !u
              ini(i,j,4) = (3.0/dd)*ww(i,j,4) !v
          end do
      end do

      !boundaries the no slip boundaries
      ! 1=>phi , 2=>pressure , 3=>u , 4=>v
      i = 0
      do j = 1,jmax-1
          ini(i,j,1) = ini(i+1,j,1)
          ini(i,j,2) = ini(i+1,j,2)
          ini(i,j,3) = 0.d0
          ini(i,j,4) = -ini(i+1,j,4)
          dd = ddL*(1.d0-ini(i,j,1))+ini(i,j,1)*ddH
          ini(i,j,5) = dd
      end do

      i = imax
      do j = 1,jmax-1
          ini(i,j,1) = ini(i-1,j,1)
          ini(i,j,2) = ini(i-1,j,2)
          ini(i,j,3) = 0.d0
          ini(i,j,4) = -ini(i-1,j,4)
          dd = ddL*(1.d0-ini(i,j,1))+ini(i,j,1)*ddH
          ini(i,j,5) = dd
      end do

      j = 0
      do i = 1,imax-1
          ini(i,j,1) = ini(i,j+1,1)
          ini(i,j,2) = ini(i,j+1,2)
          ini(i,j,3) = -ini(i,j+1,3)
          ini(i,j,4) = 0.d0
          dd = ddL*(1.d0-ini(i,j,1))+ini(i,j,1)*ddH
          ini(i,j,5) = dd
      end do

      j = jmax
      do i = 1,imax-1
          ini(i,j,1) = ini(i,j-1,1)
          ini(i,j,2) = ini(i,j-1,2)
          ini(i,j,3) = -ini(i,j-1,3)
          ini(i,j,4) = 0.d0
          dd = ddL*(1.d0-ini(i,j,1))+ini(i,j,1)*ddH
          ini(i,j,5) = dd
      end do

      ! calculating the chemical potential function MiuC:
      ! miuc = 2A(PHI)*(PHI-1)*(2PHI-1)-kD2PHI
      ! call Laplacian

      do i=1,imax-1
        do j=1,jmax-1
            ip=i-1
            in=i+1
            jp=j-1
            jn=j+1
            ini(i,j,6) = 2.d0*beta*(ini(i,j,1))*(ini(i,j,1)-1.d0)*(2.d0*ini(i,j,1)-1.d0)
            ini(i,j,6) = ini(i,j,6)-&
            &kapa*(ini(in,jn,1)+ini(in,jp,1)+ini(ip,jn,1)+ini(ip,jp,1)+4.0*ini(ip,j,1)+&
                & 4.0*ini(in,j,1)+4.0*ini(i,jp,1)+4.0*ini(i,jn,1)-20.0*ini(i,j,1))/6.0
            if (abs(ini(i,j,6))<1e-10) ini(i,j,6) = 0.d0
        end do
       end do

      !boundary for chemical potential function
      i = 0
      do j = 1,jmax-1
          ini(i,j,6) = ini(i+1,j,6)
      end do

      i = imax
      do j = 1,jmax-1
          ini(i,j,6) = ini(i-1,j,6)
      end do

      j = 0
      do i = 1,imax-1
          ini(i,j,6) = ini(i,j+1,6)
      end do

      j = jmax
      do i = 1,imax-1
          ini(i,j,6) = ini(i,j-1,6)
      end do

      call gradients

  end subroutine BoundaryCondition
