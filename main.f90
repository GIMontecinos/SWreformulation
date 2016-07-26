!!----------------------------------------------
!! Author  :  Gino I. Montecinos
!!
!! Date    :  08 - 07 -2016, CMM, Santiago.
!!
!! Purpose :  Solve the Shallow water equation
!!     The equation are based on a well balance
!!     set of equation.
!!
!! To compile
!!        gfortran main.f90 -o main
!!
!! To run
!!        ./main
!!----------------------------------------------
program  main

  implicit none 
  !! Definition of relevant global variables
  
  Integer :: imax   !! number of cells
  Integer :: nVar   !! unknowns
  Double precision :: x0,x1  !! domain [x0,x1] 
  Double precision :: tEnd !! output time
  Double precision :: cfl  !! CFL coefficient
  Double precision :: dt   !! Time step
  Double precision :: dx   !! Cell size

  Double precision :: ampl
  Double precision :: sigma
  Double precision :: epsilon
  !!  Equation parameters
  Double precision :: g = 9.81   !! gravity constant

  !! Arrays
  Double precision, allocatable :: u(:,:), q(:,:) , &
       &                     Flux(:,:),Source(:,:),x(:),&
       &                     Jump(:,:)
  Double precision, allocatable :: FL(:) 
  !!--------------------------------
  !! Auxiliary parameters
  Integer i,j,k,l,iter,i_file
  double precision time

  !!------------------------------
  !!  Parameters for ODE solver
  Integer ODE_It_Tot
  Double precision ODE_Tol

  
  !!---------------------------------
  !!  Get information from input.dat
  !!---------------------------------
  i_file = 1
  Open(i_file,File="input.dat")
  read(i_file,*)imax    !! imax
  read(i_file,*)x0      !! x0
  read(i_file,*)x1      !! x1
  read(i_file,*)nVar    !! nVar
  read(i_file,*)cfl     !! cfl
  read(i_file,*)tEnd    !! tEnd
  read(i_file,*)epsilon !! relaxation
  read(i_file,*)ampl    !! amplitude of bathymetry
  read(i_file,*)sigma   !! width of gaussian function
  read(i_file,*)ODE_It_Tot  !! Total iterations for ODE solver
  read(i_file,*)ODE_Tol  !! Total iterations for ODE solver
  Close(i_file)


  
  !!-------------------------------
  !! Allocate the memory
  !!-------------------------------

  allocate( x(imax ))
  allocate( u(nVar,imax ))
  allocate( q(nVar,imax ))
  allocate( Flux(nVar,imax ))
  allocate( Jump(nVar,imax ))
  allocate( Source(nVar,imax )) 

  allocate(FL(nVar))
  !!----------------------
  !!  Fill baricentres x.
  !!---------------------
  dx = (x1-x0) / real( imax)
  call fill_barycentre

  !!---------------
  !! fill the initial condition
  call fill_initial_condition

  time = 0.0

3 if(time < tEnd)then
     !!------------------------------
     !! Extract the time step from
     !! a CFL condition.
     call get_time_step( time, dt)

     !!------------------------------------------------
     !! Evolve the ODE
    do i = 1 , imax
       !        call solve_ODE(q(:,i),u(:,i) )
       call Num_ODE_solver(u(:,i),u(:,i),x(i),dt/2.0)
    end do

     
     !!----------------------------------
     !! Compute the numerical flux
    !! call fv_compute_rusanov_flux(Flux,FL)
     call fv_compute_muscl_flux(Flux,FL)


     
     !!     !!----------------------------
     !!     !! Compute the numerical flux
     !!    do i = 1 , imax
     !!      call PDESource( u(:,i), x(i) , Source(:,i))
     !!  end do

     !!----------------------
     !! Evolve the hyperbolic part
     do i = 1 , imax
        if(i==1)then
           u(:,i) = u(:,i)  - dt/dx*(Flux(:,i) - FL ) 
        else
           u(:,i) = u(:,i)  - dt/dx*(Flux(:,i) - Flux(:,i-1) ) 
        end if
     end do

     !!------------------------------------------------
     !! Evolve the ODE
     do i = 1 , imax
        !        call solve_ODE(q(:,i),u(:,i) )
        call Num_ODE_solver(u(:,i),u(:,i),x(i),dt/2.0 )
     end do
     
     call produce_output
     
     time = time + dt
     print*,"time= ",time," dt= ",dt
     goto 3
  end if


  deallocate(x)
  deallocate(u)
  deallocate(q)
  deallocate(Flux)
  deallocate(Source)
  deallocate(Jump)
  deallocate( FL )

  
  print*,"Program ends normally."


  stop

  !!------------------------------
  !!------------------------------
  !! Subroutines
contains

  !!-----------------------
  !! get time step
  subroutine get_time_step( time, dt)
    double precision, intent(in) :: time
    double precision, intent(out):: dt
    integer i,j
    double precision lam,lmax

    lmax = 0.0
    do i = 1 , imax
       call PDEEigenvalues(u(:,i) ,x(i),lam)
       lmax = max(lmax , lam)
    end do

    dt = cfl * dx / lmax

    if( time +dt > tEnd)then
       dt = tEnd - time
    end if
    
  end subroutine get_time_step
  
  !!----------------------
  !! Fill the barycentres
  subroutine  fill_barycentre
    !! Local variables
    integer i
    do i = 1 , imax
       x(i) = x0 + dx *(real( i) -0.5)
    end do
    
  end subroutine fill_barycentre


  !!------------------------
  !!  Minmod reconstruction
  subroutine minmod(QL,QC,QR,D)
    double precision, intent(in)  :: QL(:)
    double precision, intent(in)  :: QC(:)
    double precision, intent(in)  :: QR(:)
    double precision, intent(out) :: D(:)
    double precision DL(nVar), DR(nVar)
    integer i , j


    DL = QC - QL
    DR = QR - QC

    do j = 1, nVar
       if( DL(j) * DR(j) <=0 )then
          D(j) = 0.0
       else if( abs( DL(j) ) < abs( DR(j))  )then
          D(j) = DL(j)
       else
          D(j) = DR(j)
       end if
          
    end do
    
  end subroutine minmod

!> Purpose:
!>         mcb for reconstruction
!> Input: 
!>     QL,QC,QR : Data on the left, centre and right
!> Output:
!>       Delta : mcb limiter applied to slopes
!> 
 subroutine MCB(QL,QC,QR,Slope)
   double precision, intent(in)::QL(:),QC(:),QR(:)
   double precision, intent(out)::Slope(:)
   double precision DL(nVar),DR(nVar),beta
   integer j
   
   beta = 1.9
   DR = QR-QC
   DL = QC-QL
   
   do j=1,nVar
      if( DL(j)*DR(j) < 0.)then
         Slope(j) = 0.
         
         
      else if( 0.5*abs(DL(j) + DR(j))< beta* abs(DR(j) ) )then
         
         if( 0.5*abs(DL(j) + DR(j)) < beta*abs(DL(j) ) )then
            Slope(j) = .5*(DL(j) + DR(j))
         else
            Slope(j) = beta*DL(j)
         end if
      else
         if( abs(DR(j)) < abs(DL(j)) )then           
            Slope(j) = beta*DR(j)
         else
            Slope(j) = beta*DL(j)
         end if
      end if
   end do
   
   
 end subroutine MCB

  
  !!----------------------------
  !! Get numerical flux with
  !! Rusanov  method
  subroutine fv_compute_rusanov_flux(Flux,FluxL)
    double precision, intent(out) :: Flux(:,:),FluxL(:)
    integer i
    !!----------------------------
    !! Compute the numerical flux
    do i = 1 , imax
       if(i == imax)then
          call Rusanov( u(:,i) , u(:,i) , x(i)+dx/2.0 , Flux( :, i) )
       else
          call Rusanov( u(:,i) , u(:,i+1) , x(i)+dx/2.0 , Flux( :, i) )
       end if
       Flux(3,i) = 0.0
    end do
    
    call Rusanov( u(:,1) , u(:,1) , x(1)-dx/2.0 , FluxL  )
    FluxL(3) = 0.0
    
  end subroutine fv_compute_rusanov_flux

  !!---------------------------
  !!  Muscl-Hancock
  subroutine fv_compute_muscl_flux( Flux ,FluxL )
    double precision, intent(out) :: Flux(:,:),FluxL(:)
    double precision UL(nVar,imax) , UR(nVar,imax)
    double precision QL(nVar),QC(nVar),QR(nVar)
    double precision FL(nVar) , FR(nVar)
    integer i ,  j

    !!------------------------
    !!  Optain the increments by
    !!  using the minmod limiter
    do i = 1, imax

       if(i == 1)then

          QL = u(:,1)
          QC = u(:,1)
          QR = u(:,2)
          
       elseif( i == imax)then

          QL = u(:,imax-1)
          QC = u(:,imax)
          QR = u(:,imax)
       else
          QL = u(:,i-1)
          QC = u(:,i  )
          QR = u(:,i+1)

       end if
       
!!       call minmod(QL,QC,QR,Jump(:,i))
       call MCB(QL,QC,QR,Jump(:,i))
    end do


    !!------------------------------
    !! Get extrapolated values and
    !! evolve extrapolated values
    do i = 1, imax

       QL = u(:,i) - 0.5*Jump(:,i)
       QR = u(:,i) + 0.5*Jump(:,i)
       
       call PDEFlux( QL , x(i)-dx/2.0, FL )
       call PDEFlux( QR , x(i)+dx/2.0, FR )

       uL(:,i) = QL - dt/dx/2.0*(FR-FL)
       uR(:,i) = QR - dt/dx/2.0*(FR-FL) 
       
    end do

    !!------------------------------
    !!  Compute the numerical fluxes
    !!----------------------------
    !! Compute the numerical flux
    do i = 1 , imax
       if(i == imax)then
          call Rusanov( uR(:,i) , uR(:,i) , x(i)+dx/2.0 , Flux( :, i) )
       else
          call Rusanov( uR(:,i) , uL(:,i+1) , x(i)+dx/2.0 , Flux( :, i) )
       end if
       Flux(3,i) = 0.0
    end do

    call Rusanov( uL(:,1) , uL(:,1) , x(1)-dx/2.0 , FluxL  )
    FluxL(3) = 0.0
    
    

  end subroutine fv_compute_muscl_flux
  
  !!----------------------
  !! Bathymetry
  double precision function eta(x)
    double precision, intent(in) :: x
    double precision b
    

    b = 0.0 ;
    if( x > 0.2 .and. x < 0.35 )then
       b = ampl ;
    end if

    if( x > 0.45 .and. x < 0.65 )then
       b = ampl ;
    end if
    eta = b

    !eta = ampl * exp( - sigma * ( x - 0.5)**2 )
  end function eta
  
  !!------------------------
  !!  Definition of the initial
  !! condition
  subroutine initial_condition(x,Q)
    double precision, intent(in)  :: x
    double precision, intent(out) :: Q(:)
    double precision h, u
    double precision et

    et = 1e-6
    u = 1.0

    Q(1) = 1.0
    Q(3) = eta(x)
    h = Q(1) - Q(3)
    Q(2) =  u * h
    Q(4) = (eta(x+et)-eta(x-et) ) / (2.0*et)
    !(eta(x+dx/2.0)-eta(x-dx/2.0))/dx !0.0 !-2.0*sigma*(x-0.5)*Q(3)eppsilon
    
  end subroutine initial_condition
  
  !!-----------------------------------
  !! Call initial condition function
  !! and fill array u
  subroutine fill_initial_condition
    !! local variables
    integer i
    do i = 1 , imax
       call initial_condition(x(i),u(:,i) )
    end do
  end subroutine fill_initial_condition



  !!---------------
  !! Physical flux
  subroutine PDEFlux( Q,xi,F)
    double precision, intent(in) :: xi, Q(:)
    double precision, intent(out) :: F(:)
    double precision  b, h


    b  = Q(3)
    !b = eta(xi)
    h  = Q(1) - b
    
    F(1) = Q(2)
    F(2) = Q(2)**2 / h + 0.5*g*Q(1)**2 - g * Q(1) * b
    F(3) = 0.0
    F(4) = -b / epsilon
    
  end subroutine PDEFlux

  !!----------------------
  !! Egenvalues
  !!---------------------
  subroutine PDEEigenvalues(Q,xi,lam)
    double precision, intent(in)  :: xi, Q(:)
    double precision, intent(out) :: lam
    double precision h , b , u , c

    b = Q(3)
    h = Q(1) - b
    u = Q(2) / h

    c = sqrt( h * g )  
    lam = abs( u -c )
    lam = max( lam, abs(u) )
    lam = max( lam, abs(u+c) )
        
  end subroutine PDEEigenvalues

  
  !!--------------------
  !! Numerical source
  subroutine PDESource(Q,xi,S)
    double precision, intent(in)  :: xi,Q(:)
    double precision, intent(out) :: S(:)
    double precision bL, bR , h , bx

    !bL = eta(xi-dx/2.0)
    !bR = eta(xi+dx/2.0)
    
    !bx = (bR - bL) / dx 
  
    
    S(1) =   0.0
    S(2) = - g * Q(1) * Q(4)
    S(3) =   0.0
    S(4) = - Q(4) / epsilon
    
  end subroutine PDESource

  

  !!---------------------------
  !! Rusanov
  subroutine Rusanov(QL,QR,xi,F)
    double precision, intent(in)  :: xi,QL(:),QR(:)
    double precision, intent(out) :: F(:)
    double precision FL(nVar) , FR(nVar)
    double precision lamL,lamR,lmax

    call PDEEigenvalues(QL,xi,lamL)
    call PDEEigenvalues(QR,xi,lamR) 

    lmax = max( lamL , lamR )

    call PDEFlux(QL,xi,FL)
    call PDEFlux(QR,xi,FR)

    F = 0.5 * ( FL + FR) - 0.5 * lmax*(QR-QL)
    
  end subroutine Rusanov


  !!----------------------------
  !! Solve the ODE
  subroutine solve_ODE(U,Q)
    double precision, intent(in)   :: U(:)
    double precision, intent(out) :: Q(:)
    double precision Jac(nVar,nVar), H(nVar)
    double precision Sq(nVar), Su(nVar)
    integer it


    Q(1) = U(1);
    Q(2) = U(2) - g * U(1) * U(4) *epsilon*(1.0- exp(-dt/epsilon) );
    Q(3) = U(3) ;
    Q(4) = U(4) * exp( - dt /epsilon) ;    
    
    
  end subroutine solve_ODE

  !!----------------------------
  !! Solve the ODE numerically
  subroutine impl_ODE(U,xi,dtin,Q)
    double precision, intent(in)   :: U(:),xi,dtin
    double precision, intent(out) :: Q(:)
    double precision iJac(nVar,nVar), H(nVar)
    double precision Sq(nVar), Su(nVar) , Delta(nVar)
    double precision Tol,Err, dDt
    integer it, j , it_sol


    call PDESource(U,xi,Su)

    Q  = U;       
    do it = 1 , 160
       !!-------------------------------------
       !! Compose the jacobian matrix
       !!-------------------------------------
       iJac(1,1) = -1.0
       iJac(1,2) = 0.0
       iJac(1,3) = 0.0
       iJac(1,4) = 0.0
       iJac(2,1) = 0.5*Q(4)*dtin*g
       iJac(2,2) = -1.0
       iJac(2,3) = 0.0
       iJac(2,4) = (Q(1)*dtin*epsilon*g)/(2.0*epsilon+dtin)
       iJac(3,1) = 0.0
       iJac(3,2) = 0.0
       iJac(3,3) = -1.0
       iJac(3,4) = 0.0
       iJac(4,1) = 0.0
       iJac(4,2) = 0.0
       iJac(4,3) = 0.0
       iJac(4,4) = -(2.0*epsilon)/(2.0*epsilon+dtin)
       
       
       call PDESource(Q,xi,Sq)
       
       H = U + 0.5 * dtin *( Sq + Su) - Q
       
       Delta = matmul( iJac , H )
       
       Q = Q - Delta
       
       Err = 0.0
       do j = 1, nVar
          Err = Err + Delta(j)**2 ;
       end do
       
       if( sqrt( Err) < ODE_Tol )then

          goto 2
       end if
    end do
2      continue
       

  end subroutine impl_ODE

  !!---------------------------------------------------- 
  !! Solve the ODE numerically with an implicit method
  !! We solve Q'=S(Q), Q(0) = Q_i^n up to t = dt
  !!	  
  subroutine Num_ODE_solver(Qin,Qout,xi,ODE_tEnd)
   double precision, intent(in) :: Qin(:), xi, ODE_tEnd
    double precision, intent(out) :: Qout(:)
    double precision Q(nVar),dDt
    integer it_sol


    Q = Qin
    dDt = ODE_tEnd / real(ODE_It_Tot)
    do it_sol = 1, ODE_It_Tot
       
       call impl_ODE(Q,xi,dDt,Qout)
       Q = Qout
    end do
    
    
  end subroutine Num_ODE_solver
  
  !!----------------------------
  !! output result
  subroutine produce_output
    integer i, i_file

    i_file = 2
    Open(i_file, File = "approx.txt" )

    do i = 1 , imax
       write(i_file,*)x(i), u(3,i) , u(1,i) , u(2,i)/( u(1,i)-u(3,i)) ,u(4,i)
    end do
    
    Close(i_file)

   
  end subroutine produce_output
  
end program main

