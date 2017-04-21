! frequency unit is 100*10^14 [1/s]
! RWA is not made. Density matrix equation with damping is solved.


!This code is designed to calculate the siganl
!of FAST-CARS
!Energy Levels:
!=======================
! -------- a
!
!
!
!
!
!          ------- d
!         //
!        //
!   -------- b
!
! -------- c
!=======================



module para
implicit none

integer :: i , j, p, q

!!!!!           Unit of frequency is 10^14 HZ      !!!!!!!
! the atomic energy level frequencies
real :: w_a = 10.0, &
        w_b = 0.47,  &
        w_c = 0.0,   &
        w_d = 0.97  
! the laser frequeny
real :: nu_pu = 6.47, &      
           nu_st = 6.00, &   
           nu_pr = 7.53, &
           nu_dr = 0.50  

real :: om_pu = 0.1, & 
           om_st = 0.1, &
           om_pr = 0.1, &
           om_dr = 0.02
            

!!!!!!            Unit of time is    10^(-14) second !!!!!!
! the pulse parameters
real :: tau_pu=200.0, &
           tau_st=200.0, &
           tau_pr=500.0, &
           t0_pu=1000.0, &
           t0_st=1000.0, & 
           t0_pr=7000.0
!!    decay     
real(8) :: g_ab = 0.1, &
           g_ac = 0.1

integer, dimension(4,4) :: sigma_ab, sigma_ab_dag, sigma_ac, sigma_ac_dag 



end module para


program CARS_MATRIX
use para
implicit none
real(8), dimension(4,4) :: Hamiltonia
! additional parameters
real :: tmax = 10000.   ! the totle simulation time
real :: delta_t = 0.01
real(8) :: t = 0
integer :: N 
complex :: imag
real*8 pi

! atomic density operator
complex(8), dimension(4,4) :: rho, rho_0
complex(8), dimension(4,4) :: delta_rho_k1, &
                              delta_rho_k2, &
                              delta_rho_k3, &
                              delta_rho_k4

! the parameter for fft 
complex(8), dimension(:), allocatable :: ini, outo
integer ( kind = 4 ), parameter :: fftw_estimate = 128
integer ( kind = 4 ), parameter :: fftw_forward = -1
integer*8 plan




N = INT(tmax / delta_t) 
! unit complex number
imag=sqrt((-1.0,0.0))
! mathmatic number pi  
pi = 4.0*atan(1.0)    
allocate(ini(N))
allocate(outo(N))

!!!    Initiallized density matrix  and Hamiltonia        !!!!!!!!!!!!!!
 do i = 1 , 4 
    do j = 1 , 4
      rho_0(i,j) = (0, 0)
      Hamiltonia(i,j) = 0.0
    end do  
 end do 

! Initial state is in |c >  state --- ground state 
rho_0 (3,3) = 1
!!! ______________________________________________                    

write (*, *) "Number of step (N) is:", N
write (*, *) "Each step time (delta_t) is:", delta_t


!! loop the driving field intensity (0 -> 0.02)
do p = 0, 10, 1

om_dr = p * 0.002

write (*, *) "the driving field Rabi frequency:", om_dr

 t = 0.0
 rho = rho_0
 ! open (unit = 1, file = 'rho.txt', ACTION="write", STATUS="replace")

 do i = 1 , N 
     t = delta_t * (i - 1)  

 ! write (1 ,*) rho(3,1)  ! Hamiltonia
 ! write (1 ,*) "_______________________"

    call HD(t, Hamiltonia)   ! return the value of Hamiltonia(t)   
                             ! H(t, Hamiltonia) is the H(t_n) or H_n, corresponding to x_n
                             ! rho is the rho(t_n) or rho_n, corresponding to y_n
                             ! h is delta_t
                             ! k1 = h * f(x_n, y_n)

    delta_rho_k1 = - imag* (matmul(Hamiltonia, rho) - matmul(rho, Hamiltonia)) * delta_t 
    !!!!!!   _________________________________________________
                                    ! k2 = h* f(x_n + h / 2, y_n + k1/2) 
    call HD(t + delta_t/2, Hamiltonia)  ! return the value of Hamiltonia(t + delta_t/2) 
                                        ! 1): H(t_n + delta_t/2, Hamiltonia) corresponding to x_(n) + h/2
                                        ! 2): (H(t_n, Hamiltonia) + H(t_n, Hamiltonia))/2 corresponding to x_(n) + h/2
    
    delta_rho_k2 = - imag* (matmul(Hamiltonia, rho + delta_rho_k1/2) - matmul(rho + delta_rho_k1/2, Hamiltonia)) * delta_t 

    !!!!!!   _________________________________________________
                                    ! k3 = h* f(x_n + h / 2, y_n + k2/2) 
    delta_rho_k3 = - imag* (matmul(Hamiltonia, rho + delta_rho_k2/2) - matmul(rho + delta_rho_k2/2, Hamiltonia)) * delta_t 
    !!!!!!   _________________________________________________
                                    ! k4 = h* f(x_n + h, y_n + k3) 
    call HD(t + delta_t,  Hamiltonia)   ! return the value of Hamiltonia(t + delta_t) 
    delta_rho_k4 = - imag* (matmul(Hamiltonia, rho + delta_rho_k3)- matmul(rho + delta_rho_k3, Hamiltonia)) * delta_t 
 
 rho = rho + 1.0/6 * (delta_rho_k1 + 2 * delta_rho_k2 + 2 * delta_rho_k3 + delta_rho_k4)

 ini(i) = real(rho(3,1))

 end do 

 write (*, *) "the amplitude prob at ground state", rho(3,1)
 write (*, *) "the offdignal element of density matrix", rho(3,1)



 call dfftw_plan_dft_1d(plan,N,ini,outo,fftw_forward,fftw_estimate)
 call dfftw_execute_dft(plan,ini,outo)

 open (unit = 2, file = 'frequency.txt', ACTION="write", position = 'append')

!if (MOD(p , 1) == 0) then 
!  write (2, *) "______________________"
!end if

! this part is what we are concerning
 do i = 12653,12797,1  
  write(2,*) om_dr , 2.0 * pi * i/tmax, abs(outo(i)) !, i
 end do
 call dfftw_destroy_plan(plan)
 
 

end do 



!!!!!    the name "contains" must be here    !!!!!!!!
contains

 FUNCTION Omega_pu(t)
    real(8) :: Omega_pu, t
    Omega_pu = om_pu * exp(-(t - t0_pu)*(t - t0_pu)/(2*tau_pu*tau_pu)) * cos(nu_pu*t)
    return 
 END FUNCTION

 FUNCTION Omega_st(t)
    real(8) :: Omega_st, t
    Omega_st = om_st * exp(-(t -t0_st)*(t -t0_st)/(2*tau_st*tau_st)) * cos(nu_st*t)
    return
 END FUNCTION

  FUNCTION Omega_pr(t)
  real(8) :: Omega_pr, t
    Omega_pr = om_pr * exp(-(t -t0_pr)*(t -t0_pr)/(2*tau_pr*tau_pr)) * cos(nu_pr*t)
    return
 END FUNCTION


  FUNCTION Omega_dr(t)
  real (8) :: Omega_dr, t
    Omega_dr = om_dr * cos(nu_dr*t)
    return
 END FUNCTION

subroutine HD(t, ham)
  implicit none
  complex(8), dimension(4,4) :: rho_decay, rho
  real(8), dimension(4,4) :: ham
  real(8) :: t

    ham(1,1) = w_a
    ham(1,2) = -(Omega_st(t) + Omega_pr(t))
    ham(1,3) = -Omega_pu(t)
    ham(1,4) = 0

    ham(2,1) = ham(1,2)
    ham(2,2) = w_b
    ham(2,3) = 0
    ham(2,4) = -Omega_dr(t)

    ham(3,1) = ham(1,3)
    ham(3,2) = 0
    ham(3,3) = w_c
    ham(3,4) = 0

    ham(4,1) = 0
    ham(4,2) = ham(2,4)
    ham(4,3) = 0
    ham(4,4) = w_d

 end subroutine



end program CARS_MATRIX


