module head
    implicit none

    !   6 3 5             8 4 7
    !    \|/               \|/
    !   2-9-1 streaming   1-9-2 collision
    !    /|\               /|\
    !   7 4 8             5 3 7


	integer,parameter::ntmax= 410
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
    integer::  imax,jmax,file_unit,iter,datatime,ii
	real(8)::  niuN_LX(-1:ntmax,-1:ntmax),niuN_LY(-1:ntmax,-1:ntmax),niuN_HX(-1:ntmax,-1:ntmax),niuN_HY(-1:ntmax,-1:ntmax)
	real(8)::  niuN_LXOLD(-1:ntmax,-1:ntmax),niuN_LYOLD(-1:ntmax,-1:ntmax),niuN_HXOLD(-1:ntmax,-1:ntmax),niuN_HYOLD(-1:ntmax,-1:ntmax)
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
	real(8)::  Fr,Re,We,Bo,Ca,U0,miu_ratio,dd_ratio,T_dim
    real(8)::  rr,dia,pi,dtt,M,beta,epsilon,sigma,etha,kapa,gy ! the parameters of chan-hilliard Eqn.
    real(8)::  Fsx(-1:ntmax,-1:ntmax),Fsy(-1:ntmax,-1:ntmax),Fgx,Fgy
    real(8)::  limx2(-1:ntmax,-1:ntmax),limy2(-1:ntmax,-1:ntmax) ! laplacian scheme
    real(8)::  MG0,thetaC,thetaS ! the contact angle and slope angle
    real(8)::  alpha1,A,Ki,ThetaR,omega ! SAW parameters
    real(8)::  Us(ntmax),Ys(ntmax),Ts(ntmax),dxx,dyy
    real(8)::  Fax(-1:ntmax,-1:ntmax),Fay(-1:ntmax,-1:ntmax)
    integer::  FluidType,flag_time
    real(8)::  D0,time,start_time,end_time,L,W
    real(8)::  Diameter, H_LU, L_LU, H_phy, L_phy, Hf, Lf, t_phy, t_Lu, V_Lu, V_phy

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
      real(8)::  feq0(0:8),ffi(0:8),ff(0:8),ggi(0:8),fff(0:8),df,pp,c11,c22
      real(8)::  tmp1,tmp2,tmp3,tmp4
end module

    program main
       USE head
       implicit none
       integer::i,j

       call cpu_time(start_time)
       call settings
       call gridgen
       call fluxMesh
       call initialize
       call LBFS
       call outPuts
       !call analytical_sol
       call cpu_time(end_time)

       open(unit=510,file='time_taken.txt',action='write')
       write(510,*) 'time taken : ' , end_time - start_time

    end program main

    subroutine settings
        USE head
        USE rho_0
        implicit none

        dt = 0.1d0
        ! Specify dimensions
        imax = 401   ! Width (number of points in x direction)
        jmax = 201   ! Height (number of points in y direction)

        H_lu = jmax - 1
        L_lu = imax - 1
        ! Lf = 4D0/400 = 0.02 mm / LU
        ! Hf = 2D0/200 = 0.02 mm / LU
        ! => D0(physical) / Lf = 100 LU for Diameter of droplet
        ! for translating X and Y :: Y(physical) = H * Hf = > 4  ,
        !                            X(physical) = L * Lf = > 8
        Diameter = 2 !=>2 mm 2 micro liter charactristic length => 50 LU for droplet

        ! V_phy (m/s) * (1 Lu / x m ) * ( y sec / 1 t_Lu ) = V_Lu now we need to find Y, We have X from Hf
        !
        ! to find length in physical space we need to L_lu* hf * (mm)
        !
        ! to find velocity in physical space we need to V_lu * hf/t_Lu (mm/s)

        H_phy = 2*Diameter
        L_phy = 4*Diameter

        Hf = H_phy/H_LU ! length factor
        Lf = L_phy/L_LU

        V_LU = 0.002 ! for stability purposes
        V_phy = 0.05 ! duo to Reynolds number

        t_Lu = V_LU * Lf * 1e-3 / V_phy ! time factor

        datatime = 1
        tmax = 80000
        FluidType = 1
        flag_time = 0

    end subroutine

    subroutine initialize
        USE head
        USE rho_0
        implicit none
        integer:: in,ip,jp,jn,ios
        real(8):: factor,cs2,VL,col1,col2,col3,col4,Vel_L,Vel_W


        wa(0)=4.d0/9.d0
        do i=1,4
            wa(i)=1.d0/9.d0
        end do
        do i=5,8
            wa(i)=1.d0/36.d0
        end do
        ex(:) = (/0.0,1.0,-1.0,0.0,0.0,1.0,-1.0,-1.0,0.0/)
        ey(:) = (/0.0,0.0,0.0,1.0,-1.0,1.0,1.0,-1.0,-1.0/)

        ! for normalizing in the lattices coordination :
        ! U_pys = 1 m/s and U_LBM = 0.1 for more stability
        ! important non-dimensionalised numbers : Re , We , Fr
        ! definition of Bond number : rho*g*R^2/sigma to find the g unit in Lattice Boltzmann coordination

        dxx = maxval(lenx(:))
        dyy = maxval(leny(:))

        thetaC = 100.d0
        thetaS = 0.d0
        D0 = 2.0e-3
        Re = 100.d0
        We = 0.27d0
        Ca = We/Re !=> in order of magnitude of 1*10e-3
        Bo = 3.0
        U0 = 0.002d0
        ddH= 10.d0
        dtt = maxval(lenx)
        etha= 10.d0
        rr=50.0d0*dtt  ! radius of the drop
        dia=rr*2.0
        Cs = 1.d0/sqrt(3.d0)
        pi = 4.0*atan(1.0)
        epsilon = 3.d0*dtt  ! interface thickness
        sigma = ddH*U0**2*dia/We  ! surface tension
        !miuH  = ddH*U0*dia/Re     ! miu for the denser flow
        if (FluidType == 1) then
            miuH = ddH*U0*dia/Re
        elseif(FluidType == 2) then
            miuH  = ddH*U0**(2-q)*D0**q/Re     ! miu for the denser flow
            p = miuH / ddH ! finding the niu for power-law
        end if
        kapa=1.5*sigma*epsilon
        beta=12.0*sigma/epsilon
        miu_ratio=51.0  ! the viscosity ratio
        dd_ratio=842.0  ! density ratio
        miuL = miuH/miu_ratio
        ddL  = ddH/dd_ratio
        niuH = miuH/ddH
        niuL = miuL/ddL ! lighter flow
        ! find denser flow
        gy   = Bo*sigma/(ddH*rr**2)
        Fgx  = gy*sin(thetaS*pi/180.d0); Fgy = gy*cos(thetaS*pi/180.d0)

        !==========================================================================================
        ! setup for SAW parameters                                                                |
        ! to find alpha :                                                                         |
        !                                                                                         |
        ! alpha1 = sqrt (1 - Vl/Vw) => Vl sound speed in liquid and Vw is a property of substrace |
        ! Vw varies with frequency
        ! =========================================================================================

        Vel_L = 1497 * t_Lu / (Lf*1.0e-3)
        Vel_W = 3990 * t_Lu / (Lf*1.0e-3)
        !alpha1 = 2.47d0
        alpha1 = sqrt (1- (Vel_L/Vel_W)**2)
        ThetaR = 23.d0
        A      = 1.e-9/(Lf*1e-3)
        omega  = 2.0*pi*20e6*t_Lu
        Ki     = -1340.d0 * Lf * 1e-3

        open(unit=610, file='init.txt', status='old', action='read' , iostat = ios)

        if (ios /= 0) then

            do i=0,imax
                do j=0,jmax
                    ! initial order parameter for "phi"
                    ini(i,j,1) = 0.0
                    ww0(i,j,1) = 0.0
                    factor = sqrt((x(i)-200.0d0)**2+(y(j))**2)
                    ini(i,j,1) = 0.5d0+0.5d0*tanh((dia-2.d0*factor)/epsilon)
                    ww0(i,j,1) = ini(i,j,1)

                    ini(i,j,3) = 0.0 ! initial velocity for u
                    ini(i,j,4) = 0.0 ! initial velocity for v
                    ww0(i,j,3) = ini(i,j,3)
                    ww0(i,j,4) = ini(i,j,4)

                    ini(i,j,5) = ddH*ini(i,j,1)+(1.d0-ini(i,j,1))*ddL ! initial density
                    ini(i,j,2) = 2.d0*sigma/dia*ini(i,j,1)            ! initial pressure
                    ww0(i,j,2) = ini(i,j,2)
                end do
            end do

        else

            do i=1,imax-1
                do j=1,jmax-1
                    read(610,*,iostat=ios) col1, col2, col3, col4
                    if (ios /= 0 ) exit
                    ini(i,j,1) = col1
                    ini(i,j,2) = col2
                    ini(i,j,3) = col3
                    ini(i,j,4) = col4
                    ini(i,j,5) = ddH*ini(i,j,1)+(1.d0-ini(i,j,1))*ddL
                    ww0(i,j,1) = ini(i,j,1)
                    ww0(i,j,2) = ini(i,j,2)
                    ww0(i,j,3) = ini(i,j,5)*ini(i,j,3)*Cs**2
                    ww0(i,j,4) = ini(i,j,5)*ini(i,j,4)*Cs**2
                end do
            end do

        end if

        ! calculating the MG0 at as a initial gas mass :
        VL = 0.0d0
        do i=1,imax-1
            do j=1,jmax-1
                if (ini(i,j,1) >= 0 .AND. ini(i,j,1) <= 0.05 ) then
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
                 & 4.0*ini(in,j,1)+4.0*ini(i,jp,1)+4.0*ini(i,jn,1)-20.0*ini(i,j,1))/(6.0*dxx**2)
            end do
        end do

        !boundary for chemical potential function
        !i = 0
        ini(0,1:jmax-1,6) = ini(1,1:jmax-1,6)

        !i = imax
        ini(imax,1:jmax-1,6) = ini(imax-1,1:jmax-1,6)

        !j = 0
        ini(1:imax-1,0,6) = ini(1:imax-1,1,6)

        !j = jmax
        ini(1:imax-1,jmax,6) = ini(1:imax-1,jmax-1,6)

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
        call Force

    end subroutine initialize
    ! main program loop
    subroutine LBFS
        USE head
        USE ieee_arithmetic

        implicit none
        integer:: i,j,imid,jmid
        real(8):: res0,res1,res,wwt(ntmax,ntmax,4)

        file_unit = 1

        do iter = 1,tmax

          ! storing old variables

          wwt(1:imax-1,1:jmax-1,1) = ww0(1:imax-1,1:jmax-1,1) ! phi
          wwt(1:imax-1,1:jmax-1,2) = ww0(1:imax-1,1:jmax-1,2) ! pressure
          wwt(1:imax-1,1:jmax-1,3) = ww0(1:imax-1,1:jmax-1,3) ! rho*Cs^2*u
          wwt(1:imax-1,1:jmax-1,4) = ww0(1:imax-1,1:jmax-1,4) ! rho*Cs^2*v

          if (FluidType == 2) then
            niuN_HXOLD = niuN_HX
            niuN_HYOLD = niuN_HY
          end if

          call BC
          call Force

          rho_facexO = rho_facex
          rho_faceyO = rho_facey

          call fluxff ! calculate flux at new step

          res0 = 0.d0
          res1 = 0.d0

          do j=1,jmax-1
              do i=1,imax-1

                ww0(i,j,1) = wwt(i,j,1)-flux0(i,j,1)*(1.0d0/vol(i,j))*dt
                ww0(i,j,2) = wwt(i,j,2)-flux0(i,j,2)*(1.0d0/vol(i,j))*dt&
                & +(ini(i,j,3)*limx(i,j,5)+ini(i,j,4)*limy(i,j,5))*dt*Cs**2
                ww0(i,j,3) = wwt(i,j,3)-flux0(i,j,3)*(1.0d0/vol(i,j))*dt+(Cs**2*Fsx(i,j)+Fax(i,j))*dt
                ww0(i,j,4) = wwt(i,j,4)-flux0(i,j,4)*(1.0d0/vol(i,j))*dt+(Cs**2*Fsy(i,j)+Fay(i,j))*dt

                res0 = res0+sqrt((wwt(i,j,3)-ww0(i,j,3))**2+(wwt(i,j,4)-ww0(i,j,4))**2)
                res1 = res1+sqrt(ww0(i,j,3)**2+ww0(i,j,4)**2)

              end do
          end do

          res = res0/res1

          if (ieee_is_nan(res)) then
            print *, 'exit at iter : ' , iter , char(7)
            exit
          end if

          if(mod(iter,datatime) .eq. 0) then
            call cpu_time(time)
            write(*,'(A12, 2X, I7, 2X , A11, E12.5E2, 2X, A13, F7.1 , A1)') &
            &    'iteration:: ', iter , 'residual:: ',res, 'time taken:: ', time, 's'
            write(*,*) '==================================================================='
            write(*,*)
            if (flag_time == 1) then
                call output_time
            end if
          end if
        end do

        call outPuts

    end subroutine

    subroutine output_time
        USE head
        USE rho_0
        implicit none
        real(8):: Umax,Vmax,V
        character(len=100)::filename,directory
        directory = 'time_steps/'

        Ts(ii) = iter*dt*t_Lu

        write(filename, '(A, "TimeStep_", I8.8, ".dat")') trim(directory), iter
        open(unit=file_unit, file=trim(filename))
        open(unit=610, file='init.txt')

        write(file_unit,*) 'TITLE = "Time-dependent"'
        write(file_unit,*) 'variables= "x(mm)","y(mm)","phi","P","u(mm/s)","v(mm/s)","rho","Fax","Fay",V(mm/s)'
        write(file_unit,*) 'zone t=" time = ' , Ts(ii) , ' " i=',imax-1 ,'j=',jmax-1 ,'F=POINT'
        do j=1,jmax-1
          do i=1,imax-1
            write(file_unit,*) cenx(i)*Lf,ceny(j)*Hf,ini(i,j,1),ini(i,j,2),&
            &                  ini(i,j,3)*(Lf/t_Lu),ini(i,j,4)*(Lf/t_Lu),ini(i,j,5),Fax(i,j),Fay(i,j),&
            &                  sqrt((ini(i,j,3)*(Lf/t_Lu))**2+(ini(i,j,4)*(Lf/t_Lu))**2)
            write (610,*) ini(i,j,1),ini(i,j,2),ini(i,j,3),ini(i,j,4)
          end do
        end do
        close(file_unit)
        ii = ii + 1

        Umax = maxval(ini(:,:,3)*(Lf/t_Lu))
        Vmax = maxval(ini(:,:,4)*(Lf/t_Lu))

        open(unit=112,file='vel_avg.txt')
        write(112,*) "max Xvelocity : " , maxval(ini(:,:,3)*(Lf/t_Lu)), 'mm/s'
        write(112,*) "min Xvelocity : " , minval(ini(:,:,3)*(Lf/t_Lu)), 'mm/s'
        write(112,*) "max Yvelocity : " , maxval(ini(:,:,4)*(Lf/t_Lu)), 'mm/s'
        write(112,*) "min Yvelocity : " , minval(ini(:,:,4)*(Lf/t_Lu)), 'mm/s'
        write(112,*) "abs velocity  : " , sqrt(Umax**2+Vmax**2), 'mm/s'
        write(112,*) "iteration : " , iter , 'time : ', Ts(ii-1)
        close(112)
        close(610)

    end subroutine output_time

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

        ! calculation of the volume of bubble inter phase area (for 2-D case the area of bubble)
        do i=1,imax-1
            do j=1,jmax-1
                if (ini(i,j,1)>=0.1 .AND. ini(i,j,1) <= 0.9)then
                    VL = VL + vol(i,j)
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
                if (ini(i,j,1) >= 0.1 .AND. ini(i,j,1) <= 0.9) then
                    temp = temp + vol(i,j)*M*(ini(in,jn,6)+ini(in,jp,1)+ini(ip,jn,6)+ini(ip,jp,6)+4.0*ini(ip,j,6)+&
                         & 4.0*ini(in,j,6)+4.0*ini(i,jp,6)+4.0*ini(i,jn,6)-20.0*ini(i,j,6))/6.0
                end if
            end do
        end do

        do i=1,imax-1
            do j=1,jmax-1
                if(ini(i,j,1) >= 0.1 .AND. ini(i,j,1) <= 0.9) then
                    qs(i,j) = 1.d0/(VL)*((mL-M0)/((ddH-ddL)*dt) - temp)
                else
                    qs(i,j) = 0.d0
                end if
            end do
        end do

    end subroutine mass

    subroutine analytical_sol

        implicit none
        real(8):: thetaC,Weqs(0:200),Heqs(0:200),factor
        integer::i
        real(8)::pi
        pi = 4.0*atan(1.0)

        do i=1,180
            thetaC = i*pi/180.d0
            factor = sqrt((pi) / (2 * thetaC - sin(2 * thetaC)))
            Weqs(i) = sin(thetaC) * factor
            Heqs(i) = ((1 - cos(thetaC)) / 2) * factor
        end do

        open(unit=10,file='Analytical.dat')
        write(10,*) 'TITLE = "DROPLET ANALYTICAL SOLUTION"'
        write(10,*) 'VARIABLES = "thetaC" , "Weqs" , "Heqs"'
        write(10,*) 'ZONE t="Analytical"'
        do i=1,180
            write(10,*) i, Weqs(i) , Heqs(i)
        end do

        close(10)

    end subroutine
    ! calculate the multiple relaxation time
    subroutine MRT
        USE head
        implicit none
        real(4)::a1
        integer::i,j,k

        ! matrix M
        tm(0,:) = (/0.0,0.0,0.0,0.0,0.0,1.0,-1.0,1.0,-1.0/)
        tm(1,:) = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
        tm(2,:) = (/-4.0,-1.0,-1.0,-1.0,-1.0,2.0,2.0,2.0,2.0/)
        tm(3,:) = (/4.0,-2.0,-2.0,-2.0,-2.0,1.0,1.0,1.0,1.0/)
        tm(4,:) = (/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
        tm(5,:) = (/0.0,-2.0,0.0,2.0,0.0,1.0,-1.0,-1.0,1.0/)
        tm(6,:) = (/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
        tm(7,:) = (/0.0,0.0,-2.0,0.0,2.0,1.0,1.0,-1.0,-1.0/)
        tm(8,:) = (/0.0,1.0,-1.0,1.0,-1.0,0.0,0.0,0.0,0.0/)

        a1 = 1.0/36.0

        ! inversion of M matrix
        tminv(0,:) = (/4.0*a1,2.0*a1,a1,6.0*a1,3.0*a1,-6.0*a1,-3.0*a1,0.0,-9.0*a1/)
        tminv(1,:) = (/4.0*a1,-4.0*a1,4.0*a1,0.0,0.0,0.0,0.0,0.0,0.0/)
        tminv(2,:) = (/4.0*a1,-a1,-2.0*a1,6.0*a1,-6.0*a1,0.0,0.0,9.0*a1,0.0/)
        tminv(3,:) = (/4.0*a1,-a1,-2.0*a1,0.0,0.0,6.0*a1,-6.0*a1,-9.0*a1,0.0/)
        tminv(4,:) = (/4.0*a1,-a1,-2.0*a1,-6.0*a1,6.0*a1,0.0,0.0,9.0*a1,0.0/)
        tminv(5,:) = (/4.0*a1,-a1,-2.0*a1,0.0,0.0,-6.0*a1,6.0*a1,-9.0*a1,0.0/)
        tminv(6,:) = (/4.0*a1,2.0*a1,a1,6.0*a1,3.0*a1,6.0*a1,3.0*a1,0.0,9.0*a1/)
        tminv(7,:) = (/4.0*a1,2.0*a1,a1,-6.0*a1,-3.0*a1,6.0*a1,3.0*a1,0.0,-9.0*a1/)
        tminv(8,:) = (/4.0*a1,2.0*a1,a1,-6.0*a1,-3.0*a1,-6.0*a1,-3.0*a1,0.0,9.0*a1/)

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
        sm(0) = s9

        ! calculate the updated M-1*S
        do i=1,9
            do j=1,9
                stmiv(i,j) = tminv(i,j)*sm(j)
            end do
        end do

    end subroutine
    ! the final outputs
    subroutine outPuts
        USE head
        implicit none
        integer::i,j,in,im,jn,jm
        integer::mid,itmp(ntmax),jtmp(ntmax)
        real(8):: Weqs,Heqs

        itmp = -1000.d0
        jtmp = -1000.d0

        mid = (imax-1)/2
        open(2,file='droplet.dat')
        write(2,*) 'variables= "x(mm)","y(mm)","phi","P","u(mm/s)","v(mm/s)","rho","miu","V(mm/s)"'
        write(2,*) 'zone i=',imax-1 ,'j=',jmax-1 ,'F=POINT'
        do j=1,jmax-1
            do i=1,imax-1
                write(2,*) cenx(i)*Lf,ceny(j)*Hf,ini(i,j,1),ini(i,j,2),&
                &          ini(i,j,3)*(Lf/t_Lu),ini(i,j,4)*(Lf/t_Lu),ini(i,j,5),ini(i,j,6),&
                &          sqrt((ini(i,j,3)*(Lf/t_Lu))**2+(ini(i,j,4)*(Lf/t_Lu))**2)
            end do
        end do
        close(2)

        ! finding Heq and Weq
        do i=1,imax-1
            if(ini(i,1,1) >= 0.5d0) then
                itmp(i) = i
            end if
        end do

        do j=1,jmax-1
            if(ini(mid,j,1) >= 0.5d0) then
                jtmp(j) = j
            end if
        end do

        im = int(maxval(itmp))
        in = int(minval(abs(itmp)))

        jm = int(maxval(jtmp))
        jn = int(minval(abs(jtmp)))

        Weqs = (cenx(im)-cenx(in))/dia
        Heqs = (ceny(jm)-ceny(jn))/dia

        open(3,file='present.dat')
        write(3,*) 'variables= "theta" , "Weqs" , "Heqs"'
        write(3,*) 'zone t="present"'
        write(3,*) thetaC , Weqs , Heqs

        close(3)

    end subroutine outPuts
    ! calculate the equilibrium function using the macro on the interfaces
    subroutine feq(dd,pp,uu,vv,ff)
        USE head
        implicit none
        integer:: n
        real(8):: dd,uu,vv,pp,ff(0:8)
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
        ff(0) = c1*(pp - dd*0.5d0*(uu*uu + vv*vv))

        return
    end subroutine
    ! calculate the equilibrium g function using the macro on the interfaces
    subroutine geq(phi,miu,uu,vv,gg)
        implicit none

        integer::n
        real(8)::phi,miu,uu,vv,gg(0:8)
        real(8)::a1,b1,c1,u2,uu2,vv2,etha

        etha = 10.0d0
        a1 = 1.d0/9.d0
        b1 = 1.d0/36.d0
        c1 = a1*4.d0

        gg(0) = phi + 3.0d0*miu*etha*(1.d0 - c1)
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
            case (0)
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
            case (0)
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
    ! calculating the strain rate then finding the s8 in x and y direction
    subroutine strain(feq,ff,niuIN,xx,rho,niuOUT)
        USE head
        implicit none
        integer::i,j,k,kk,nn
        real(8)::tmp1,tmp2,tmp3,tmppp,tao
        real(8)::M_bar(0:8)
        real(8),intent(IN)::feq(0:8),ff(0:8),niuIN,xx,rho !xx is facex(1) or facey(1) , rho is ini(i,j,1)
        real(8),intent(OUT)::niuOUT
        real(8)::  mom(0:8),SS(3),gama_dot,ftmp(0:8)
        ! we must use the previous time step Niu for reconstruction of gama-dot
        ! calculate Gama dot
        SS   = 0.d0
        tao  = niuIN/(Cs**2*xx)
        tmp2 = -0.5d0/(rho*Cs**2*xx)

        do i=0,8
            ftmp(i) = (feq(i)-ff(i))
        end do

        SS(1)=tmp2*(ex(0)*ex(0)*ftmp(0)+ex(1)*ex(1)*ftmp(1)+&
            &ex(2)*ex(2)*ftmp(2)+ex(3)*ex(3)*ftmp(3)+ex(4)*ex(4)*ftmp(4)+&
            &ex(5)*ex(5)*ftmp(5)+ex(6)*ex(6)*ftmp(6)+ex(7)*ex(7)*ftmp(7)+&
            &ex(8)*ex(8)*ftmp(8))!xx

        SS(2)=tmp2*(ex(0)*ey(0)*ftmp(0)+ex(1)*ey(1)*ftmp(1)+&
            &ex(2)*ey(2)*ftmp(2)+ex(3)*ey(3)*ftmp(3)+ex(4)*ey(4)*ftmp(4)+&
            &ex(5)*ey(5)*ftmp(5)+ex(6)*ey(6)*ftmp(6)+ex(7)*ey(7)*ftmp(7)+&
            &ex(8)*ey(8)*ftmp(8))!xy

        SS(3)=tmp2*(ey(0)*ey(0)*ftmp(0)+ey(1)*ey(1)*ftmp(1)+&
            &ey(2)*ey(2)*ftmp(2)+ey(3)*ey(3)*ftmp(3)+ey(4)*ey(4)*ftmp(4)+&
            &ey(5)*ey(5)*ftmp(5)+ey(6)*ey(6)*ftmp(6)+ey(7)*ey(7)*ftmp(7)+&
            &ey(8)*ey(8)*ftmp(8))!yy

        gama_dot = 2.d0*sqrt((SS(1)*SS(1)+2.0*SS(2)*SS(2)+SS(3)*SS(3)))
        if (gama_dot == 0.0) gama_dot = 1e-5
        !if (gama_dot >= 10)  gama_dot = 1.0
        niuOUT   = p*(abs(gama_dot))**(q-1.0)

        return
    end subroutine strain
    ! calculating the fluxes variables using interpolation
    subroutine fluxff

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
      real(8)::  feq0(0:8),geq0(0:8),ffi(0:8),ggi(0:8),fs(0:8),gs(0:8)
      real(8)::  tmp1,tmp2,tmp3,tmp4m,miufx(0:8),miufy(0:8)
      real(8)::  niuOUT


      flux0(0:imax,0:jmax,1) = 0.d0
      flux0(0:imax,0:jmax,2) = 0.d0
      flux0(0:imax,0:jmax,3) = 0.d0
      flux0(0:imax,0:jmax,4) = 0.d0


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

              nn = 0
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
              pp = ffi(1)+ffi(2)+ffi(3)+ffi(4)+ffi(5)+ffi(6)+ffi(7)+ffi(8)+ffi(0)
              phi= ggi(1)+ggi(2)+ggi(3)+ggi(4)+ggi(5)+ggi(6)+ggi(7)+ggi(8)+ggi(0)
              if(phi>1.d0) phi = 1.d0
              if(phi<0.d0) phi = 0.d0
              if(phi<1e-8) phi = 0.d0
              dd = ddL*(1.d0-phi)+phi*ddH
              uu = 1.d0/(dd*Cs**2)*(ffi(1)-ffi(2)+ffi(5)-ffi(6)-ffi(7)+ffi(8))
              vv = 1.d0/(dd*Cs**2)*(ffi(3)-ffi(4)+ffi(5)+ffi(6)-ffi(7)-ffi(8))
              miu = wa(1)*miufx(1)+wa(2)*miufx(2)+wa(3)*miufx(3)+wa(4)*miufx(4)+wa(5)*miufx(5)&
              &    +wa(6)*miufx(6)+wa(7)*miufx(7)+wa(8)*miufx(8)+wa(9)*miufx(0)
              if (abs(miu) < 1e-10 ) miu = 0.d0

              ! calculate the equilibrium functions on the center of interfaces:
              call feq(dd,pp,uu,vv,feq0)
              call geq(phi,miu,uu,vv,geq0)

              if(FluidType == 2) then ! for Non-Newtonian fluid:
                  !=====================================
                  ! for the NON Newtonian
                  !=====================================
                  ! update the relaxation time
                  call strain(feq0,ffi,niuN_HXOLD(i,j),facex(i,j,1),dd,niuOUT)
                  niuN_HX(i,j) = niuOUT
                  taoxH(i,j)   = niuN_HX(i,j)/(facex(i,j,1)*Cs**2)+0.5d0

                  taof = (taoxH(i,j)*taoxL(i,j))/(phi*taoxL(i,j)+(1.d0-phi)*taoxH(i,j))
                  taog = taogx(i,j)

              elseif(FluidType == 1) then ! for Newtonian fluid
                  !=====================================
                  ! for the Newtonian
                  !=====================================
                  ! update the relaxation time
                  taof = (taoxH(i,j)*taoxL(i,j))/(phi*taoxL(i,j)+(1.d0-phi)*taoxH(i,j))
                  taog = taogx(i,j)
              end if

              do k=0,8
                fs(k) = feq0(k)-(taof-0.5d0)*(feq0(k)-ffi(k))
                gs(k) = geq0(k)-(taog-0.5d0)*(geq0(k)-ggi(k))
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

              nn = 0
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
              pp = ffi(1)+ffi(2)+ffi(3)+ffi(4)+ffi(5)+ffi(6)+ffi(7)+ffi(8)+ffi(0)
              phi= ggi(1)+ggi(2)+ggi(3)+ggi(4)+ggi(5)+ggi(6)+ggi(7)+ggi(8)+ggi(0)
              if(phi>1.d0)  phi = 1.d0
              if(phi<0.d0)  phi = 0.d0
              if(phi<1e-8) phi = 0.d0
              dd = ddL*(1.d0-phi)+phi*ddH
              uu = 1.d0/(dd*Cs**2)*(ffi(1)-ffi(2)+ffi(5)-ffi(6)-ffi(7)+ffi(8))
              vv = 1.d0/(dd*Cs**2)*(ffi(3)-ffi(4)+ffi(5)+ffi(6)-ffi(7)-ffi(8))
              miu = wa(1)*miufy(1)+wa(2)*miufy(2)+wa(3)*miufy(3)+wa(4)*miufy(4)+wa(5)*miufy(5)&
              &    +wa(6)*miufy(6)+wa(7)*miufy(7)+wa(8)*miufy(8)+wa(9)*miufy(0)
              if (abs(miu) < 1e-10 ) miu = 0.d0

              call feq(dd,pp,uu,vv,feq0)
              call geq(phi,miu,uu,vv,geq0)

              if (FluidType == 2 ) then
                  !=====================================
                  ! for the Non-Newtonian
                  !=====================================
                  ! update the relaxation time
                  call strain(feq0,ffi,niuN_HYOLD(i,j),facey(i,j,1),dd,niuOUT)
                  niuN_HY(i,j) = niuOUT
                  taoyH(i,j)   = niuN_HY(i,j)/(facey(i,j,1)*Cs**2)+0.5d0

                  taof = (taoyH(i,j)*taoyL(i,j))/(phi*taoyL(i,j)+(1.d0-phi)*taoyH(i,j))
                  taog = taogy(i,j)

              elseif(FluidType == 1) then
                  !=====================================
                  ! for the Newtonian
                  !=====================================
                  ! update the relaxation time
                  taof = (taoyH(i,j)*taoyL(i,j))/(phi*taoyL(i,j)+(1.d0-phi)*taoyH(i,j))
                  taog = taogy(i,j)
              end if

              do k=0,8
                fs(k) = feq0(k)-(taof-0.5d0)*(feq0(k)-ffi(k))
                gs(k) = geq0(k)-(taog-0.5d0)*(geq0(k)-ggi(k))
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

    subroutine Force
    ! calculating the gradients 2nd order
      USE head
      implicit none
      integer::i,j
      real(8):: Xtmp,Ytmp

      ! F = -phi * grad miu
      do i=1,imax-1
        do j=1,jmax-1
            Fsx(i,j) = -ini(i,j,1)*limx(i,j,6)
            Fsy(i,j) = -ini(i,j,1)*limy(i,j,6)

            !==============================
            ! calculate the acoustic force:

            Xtmp = cenx(i)!*2.5e-5
            Ytmp = ceny(j)!*2.5e-5

            ! 1e4 , 1e5 , 1e6 , 1e7

            if(ini(i,j,1) >= 0.5) then
                Fax(i,j) = -1.d0*ini(i,j,5)*(1.d0+alpha1**2)*A**2*omega**2*Ki*&
                &exp(2.d0*(Ki*Xtmp+alpha1*Ki*Ytmp))*sin(ThetaR*pi/180.d0)*1e4

                Fay(i,j) = -1.d0*ini(i,j,5)*(1.d0+alpha1**2)*A**2*omega**2*Ki*&
                &exp(2.d0*(Ki*Xtmp+alpha1*Ki*Ytmp))*cos(ThetaR*pi/180.d0)*1e4

                !print *, Fay(i,j), Fax(i,j)
            else
                Fax(i,j) = 0.0
                Fay(i,j) = 0.0
            end if

            !==============================
        end do
      end do


        ! handling the boundary nodes :
        ! i=0 left B-C
        Fsx(0,1:jmax-1) = Fsx(1,1:jmax-1)
        Fsy(0,1:jmax-1) = Fsy(1,1:jmax-1)
        Fax(0,1:jmax-1) = Fax(1,1:jmax-1)
        Fay(0,1:jmax-1) = Fay(1,1:jmax-1)

        ! i = imax right B-C
        Fsx(imax,1:jmax-1) = Fsx(imax-1,1:jmax-1)
        Fsy(imax,1:jmax-1) = Fsy(imax-1,1:jmax-1)
        Fax(imax,1:jmax-1) = Fax(imax-1,1:jmax-1)
        Fay(imax,1:jmax-1) = Fay(imax-1,1:jmax-1)

        ! j=0 upper B-C
        Fsx(1:imax-1,0) = Fsx(1:imax-1,1)
        Fsy(1:imax-1,0) = Fsy(1:imax-1,1)
        Fax(1:imax-1,0) = Fax(1:imax-1,1)
        Fay(1:imax-1,0) = Fay(1:imax-1,1)

        ! j=jmax bottom B-C
        Fsx(1:imax-1,jmax) = Fsx(1:imax-1,jmax-1)
        Fsy(1:imax-1,jmax) = Fsy(1:imax-1,jmax-1)
        Fax(1:imax-1,jmax) = Fax(1:imax-1,jmax-1)
        Fay(1:imax-1,jmax) = Fay(1:imax-1,jmax-1)

    end subroutine Force

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

      !i = 0
      limx(0,1:jmax-1,1) = (ini(1,1:jmax-1,1)-ini(0,1:jmax-1,1))/lenx(0)
      limx(0,1:jmax-1,2) = (ini(1,1:jmax-1,2)-ini(0,1:jmax-1,2))/lenx(0)
      limx(0,1:jmax-1,3) = (ini(1,1:jmax-1,3)-ini(0,1:jmax-1,3))/lenx(0)
      limx(0,1:jmax-1,4) = (ini(1,1:jmax-1,4)-ini(0,1:jmax-1,4))/lenx(0)
      limx(0,1:jmax-1,5) = (ini(1,1:jmax-1,5)-ini(0,1:jmax-1,5))/lenx(0)
      limx(0,1:jmax-1,6) = (ini(1,1:jmax-1,6)-ini(0,1:jmax-1,6))/lenx(0)

      !i = imax
      limx(imax,1:jmax-1,1) = (ini(imax,1:jmax-1,1)-ini(imax-1,1:jmax-1,1))/lenx(imax)
      limx(imax,1:jmax-1,2) = (ini(imax,1:jmax-1,2)-ini(imax-1,1:jmax-1,2))/lenx(imax)
      limx(imax,1:jmax-1,3) = (ini(imax,1:jmax-1,3)-ini(imax-1,1:jmax-1,3))/lenx(imax)
      limx(imax,1:jmax-1,4) = (ini(imax,1:jmax-1,4)-ini(imax-1,1:jmax-1,4))/lenx(imax)
      limx(imax,1:jmax-1,5) = (ini(imax,1:jmax-1,5)-ini(imax-1,1:jmax-1,5))/lenx(imax)
      limx(imax,1:jmax-1,6) = (ini(imax,1:jmax-1,6)-ini(imax-1,1:jmax-1,6))/lenx(imax)


      !j = 0
      limy(1:imax-1,0,1) = (ini(1:imax-1,1,1)-ini(1:imax-1,0,1))/leny(0)
      limy(1:imax-1,0,2) = (ini(1:imax-1,1,2)-ini(1:imax-1,0,2))/leny(0)
      limy(1:imax-1,0,3) = (ini(1:imax-1,1,3)-ini(1:imax-1,0,3))/leny(0)
      limy(1:imax-1,0,4) = (ini(1:imax-1,1,4)-ini(1:imax-1,0,4))/leny(0)
      limy(1:imax-1,0,5) = (ini(1:imax-1,1,5)-ini(1:imax-1,0,5))/leny(0)
      limy(1:imax-1,0,6) = (ini(1:imax-1,1,6)-ini(1:imax-1,0,6))/leny(0)


      !j = jmax
      limy(1:imax-1,jmax,1) = (ini(1:imax-1,jmax,1)-ini(1:imax-1,jmax-1,1))/leny(jmax)
      limy(1:imax-1,jmax,2) = (ini(1:imax-1,jmax,2)-ini(1:imax-1,jmax-1,2))/leny(jmax)
      limy(1:imax-1,jmax,3) = (ini(1:imax-1,jmax,3)-ini(1:imax-1,jmax-1,3))/leny(jmax)
      limy(1:imax-1,jmax,4) = (ini(1:imax-1,jmax,3)-ini(1:imax-1,jmax-1,3))/leny(jmax)
      limy(1:imax-1,jmax,5) = (ini(1:imax-1,jmax,3)-ini(1:imax-1,jmax-1,3))/leny(jmax)
      limy(1:imax-1,jmax,6) = (ini(1:imax-1,jmax,6)-ini(1:imax-1,jmax-1,6))/leny(jmax)

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
                limx2(i,j) = (2.d0/(h(i)+h(i-1)))*((ini(ip,j,6)-ini(i,j,6))/h(i)-(ini(i,j,6)-ini(in,j,6))/h(i-1))
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
          dy = leny(j)
          do i = 1, imax
              taoxH(i,j) = niuH/(Cs**2*facex(i,j,1))+0.5d0
              taoxL(i,j) = niuL/(Cs**2*facex(i,j,1))+0.5d0
              taogx(i,j) = 1.2d0!M/(etha*slen)
          end do
      end do

      do i = 1,imax-1
          dx= lenx(i)
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
              slen = 0.25*min(dx,dy,dy1)! dy=dt for lattices
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
      real   ::rp,eta,rlg_x,rlg_y,disr,tmp

        pi  = 4.d0 * atan(1.d0)
        eta = 500000.0d0

        ! Two separate coefficients for x and y directions
        rlg_x = (real(imax) - 1.0) / 2.0 !L/2.0   ! Coefficient for x-direction
        rlg_y = (real(jmax) - 1.0) / 2.0 !w/2.0   ! Coefficient for y-direction

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

    subroutine BC
      USE head
      implicit none

      integer:: i,j,in,ip,jp,jn
      real(8):: dd,phi,pp,miu,thetaD
      real(8):: Vcl(-1:ntmax)


      ! ini(1)=>phi , ini(2)=>p , ini(3)=>u , ini(4)=>v ,ini(5)=>rho, ini(6)=>miu
      do j = 1,jmax-1
          do i = 1,imax-1
              phi = ww0(i,j,1)
              if (phi < 0.0d0 ) phi = 0.0d0
              if (phi > 1.0d0 ) phi = 1.0d0
              if (phi < 1e-8  ) phi = 0.d0
              ini(i,j,1) = phi !order parameter
              dd = ddL*(1.0d0-phi)+phi*ddH !rho density
              ini(i,j,5) = dd ! store the density for gradient
              pp = ww0(i,j,2)
              ini(i,j,2) = pp !pressure
              ini(i,j,3) = (3.0/dd)*ww0(i,j,3) !u
              ini(i,j,4) = (3.0/dd)*ww0(i,j,4) !v
          end do
      end do

      !boundaries the no slip boundaries
      ! 1=>phi , 2=>pressure , 3=>u , 4=>v
      ! no slip BC for the velocities
      ! n.grad(phi) = 4/eps*cos(thetaC)*phi*(1-phi)

      !i = 0
      ini(0,1:jmax-1,1) = ini(imax-1,1:jmax-1,1)
      ini(0,1:jmax-1,2) = ini(imax-1,1:jmax-1,2)
      ini(0,1:jmax-1,3) = ini(imax-1,1:jmax-1,3)
      ini(0,1:jmax-1,4) = ini(imax-1,1:jmax-1,4)
      ini(0,1:jmax-1,5) = ddL*(1.d0-ini(0,1:jmax-1,1))+ini(0,1:jmax-1,1)*ddH

      !i = imax
      ini(imax,1:jmax-1,1) = ini(1,1:jmax-1,1)
      ini(imax,1:jmax-1,2) = ini(1,1:jmax-1,2)
      ini(imax,1:jmax-1,3) = ini(1,1:jmax-1,3)
      ini(imax,1:jmax-1,4) = ini(1,1:jmax-1,4)
      ini(imax,1:jmax-1,5) = ddL*(1.d0-ini(imax,1:jmax-1,1))+ini(imax,1:jmax-1,1)*ddH

      !j = 0
      ! calculating dynamic contact angle modeling :
      ! first finding gradient of order parameter d(phi)/dx and dy :

      call V_surface(Vcl)

      do i=1,imax-1
        call Dynamic_CA(Vcl(i),thetaD)
        ini(i,0,1) = (4.d0/epsilon*cos(thetaD*pi/180)*ini(i,1,1)*(1.d0-ini(i,1,1)))*(leny(1))+ini(i,1,1)
      end do

      !ini(1:imax-1,0,1) = (4.d0/epsilon*cos(thetaC*pi/180)*ini(1:imax-1,1,1)*(1.d0-ini(1:imax-1,1,1)))*(leny(1))+ini(1:imax-1,1,1)
      ini(1:imax-1,0,2) = ini(1:imax-1,1,2)
      ini(1:imax-1,0,3) = -ini(1:imax-1,1,3)
      ini(1:imax-1,0,4) = 0.d0
      ini(1:imax-1,0,5) = ddL*(1.d0-ini(1:imax-1,0,1))+ini(1:imax-1,0,1)*ddH

      !j = jmax
      ini(1:imax-1,jmax,1) = ini(1:imax-1,jmax-1,1)
      ini(1:imax-1,jmax,2) = ini(1:imax-1,jmax-1,2)
      ini(1:imax-1,jmax,3) = ini(1:imax-1,jmax-1,3)
      ini(1:imax-1,jmax,4) = ini(1:imax-1,jmax-1,4)
      ini(1:imax-1,jmax,5) = ddL*(1.d0-ini(1:imax-1,jmax,1))+ini(1:imax-1,jmax,1)*ddH

      ! calculating the chemical potential function MiuC:
      ! miuc = 2A(PHI)*(PHI-1)*(2PHI-1)-kD2PHI
      ! call Laplacian
      ! n.grad(miu) = 0 for the BC

      do i=1,imax-1
        do j=1,jmax-1
            ip=i-1
            in=i+1
            jp=j-1
            jn=j+1
            ini(i,j,6) = 2.d0*beta*(ini(i,j,1))*(ini(i,j,1)-1.d0)*(2.d0*ini(i,j,1)-1.d0)
            ini(i,j,6) = ini(i,j,6)-&
            &kapa*(ini(in,jn,1)+ini(in,jp,1)+ini(ip,jn,1)+ini(ip,jp,1)+4.0*ini(ip,j,1)+&
                & 4.0*ini(in,j,1)+4.0*ini(i,jp,1)+4.0*ini(i,jn,1)-20.0*ini(i,j,1))/(6.0*dxx**2)
            if (abs(ini(i,j,6))<1e-10) ini(i,j,6) = 0.d0
        end do
       end do

      !boundary for chemical potential function
      !i = 0
      ini(0,1:jmax-1,6) = ini(1,1:jmax-1,6)

      !i = imax
      ini(imax,1:jmax-1,6) = ini(imax-1,1:jmax-1,6)

      !j = 0
      ini(1:imax-1,0,6) = ini(1:imax-1,1,6)

      !j = jmax
      ini(1:imax-1,jmax,6) = ini(1:imax-1,jmax-1,6)

      call gradients

    end subroutine BC

    subroutine V_surface(Vcl)
        USE head
        USE rho_0
        implicit none
        real(8),intent(out):: Vcl(-1:ntmax)
        real(8):: VclX(-1:ntmax),VclY(-1:ntmax)
        real(8):: NormVecX(-1:ntmax),NormVecY(-1:ntmax)
        real(8):: dphidx,dphidy,mag

        ! Vcl = n * u_cell * n
        ! Vcl = u - (u . n)n => (a*b)*c = b(a . c) - a(b . c)
        ! calculation on J=1 cells

        do i=1,imax-1
            dphidx = limx(i,1,1)
            dphidy = limy(i,1,1)
            mag    = sqrt(dphidx**2+dphidy**2)
            NormVecX(i) = 0.0
            NormVecY(i) = 0.0
            VclX(i)     = ini(i,1,3) - (ini(i,1,3)*NormVecX(i)-ini(i,1,4)*NormVecY(i))*NormVecX(i)
            VclY(i)     = ini(i,1,4) - (ini(i,1,3)*NormVecX(i)-ini(i,1,4)*NormVecY(i))*NormVecY(i)
            Vcl(i)      = VclX(i)!sqrt(VclX(i)**2 + VclY(i)**2)
        end do

    end subroutine

    subroutine Dynamic_CA(Vcl,thetaD)
        USE head
        USE rho_0

        implicit none
        real(8),intent(in):: Vcl
        real(8),intent(out):: thetaD
        integer:: func
        real   :: Ca_cl, kr, ka, theta_ad, theta_re

        kr = 9e-8
        ka = 9e-9

        thetaC = 90.0
        theta_ad = 114.0
        theta_re = 93.0

        func = 2 ! 1: Quasi-dynamic 2: Yokoi

        select case(func)
            case (1)
                if(Vcl < 0) then
                    thetaD = theta_re  ! theta re
                else
                    thetaD = theta_ad ! theta ad
                end if
            case (2)
                Ca_cl = miuH*Vcl/sigma
                if (Vcl <0) then
                    thetaD = max(thetaC + (Ca_cl/kr)**(1./3.) , theta_re)
                else
                    thetaD = min(thetaC + (Ca_cl/ka)**(1./3.) , theta_ad)
                end if
        end select

    end subroutine

