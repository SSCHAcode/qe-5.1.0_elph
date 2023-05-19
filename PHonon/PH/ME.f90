! Program originally written by Miguel Borinaga during his PhD thesis
! to solve the isotropic Migdal Eliashberg equations

program me

implicit none

!Integers
integer :: t_number ! Number of different temperatures
integer :: counter,reason
integer :: wmax_int ! Maximum freq integer (for array allocation and integration)
integer :: mats_max ! Maximum Matsubara frequency
integer :: nit_max  ! Maximum number of iterations for the self-consistent loop
integer :: i,j,k,l,t
integer :: nit

!Real scalars
double precision :: t_first,t_last,t_step      !Temperatures
double precision :: mu_star                    ! Coulomb pseudo mustar
double precision :: wc_cutoff                  ! Matsubara freq. cutoff
double precision :: realnumber                 ! Any real number (for data reading)
double precision :: w_step
double precision :: gap_guess
double precision :: s, ss, t1, t2
double precision :: mixing                     ! Mixing parameter between old an new gap 

!Real arrays

double precision , dimension(:,:) , allocatable :: w_a2f_array      ! Array for storing (w,a2F(w)) data of Eliashberg's function
double precision , dimension(:,:) , allocatable :: lambda_array     ! Array for storing (lambda(wi,wj)) 

double precision , dimension(:) , allocatable :: t_array            ! Temperature array
double precision , dimension(:) , allocatable :: delta_array, delta_old        ! Energy gap array in imaginary axis
double precision , dimension(:) , allocatable :: z_array            ! Z function array in imaginary axis
double precision , dimension(:) , allocatable :: phi_array          ! Phi function array in imaginary axis
double precision , dimension(:) , allocatable :: deltazero_array    ! Energy gap at first matsubara freq. vs temperature
double precision , dimension(:) , allocatable :: aux_array          ! Auxiliar array
double precision , dimension(:,:) , allocatable :: matsubara_freq   ! List with Matsubara frequencies at each temperature

!Characters

character(len=20) :: a2f_filename, lambda_filename,aux_string,x1
character(len=200) :: auxxx

! Parameters:

double precision :: pi=3.14159265358979323846
double precision :: KtoHartree=3.1668d-06
double precision :: meVtoHartree=3.674933d-05

!input data
NAMELIST / inputme /&

t_number, &     ! Number of different temperatures in which to do the calculation
t_first, &      ! Temperature starting point 
t_last, &       ! Temperature last point
mu_star, &      ! Value of the Coulomb pseudopotential.
wc_cutoff, &    ! Cutoff for the matsubara frequencies.
gap_guess, &    ! First guess for the energy gap (meV)
a2f_filename, & ! Filename for a2F(w) data
mixing, &       ! Mixing parameter for the self-consistent loop
nit_max, &      ! Maximum number of iterations for the self-consistent loop 
mats_max        ! Maximum number of Matsubara frequencies. The total number of Matsubara freqs. is (-mat_max,mats_max) 

!default values

t_number = 1
t_first = 10.0
t_last = 10.0
mu_star = 0.10
wc_cutoff = 0.02
gap_guess = 10.0d0
a2f_filename = 'a2F.dat'
mixing = 0.2
nit_max = 100000
mats_max = -1


print *, ''
print *, ' ****************************************'
print *, ' *                                      *'
print *, ' *         READ INPUT DATA              *'
print *, ' *                                      *'
print *, ' ****************************************'
print *, ''


read(*,inputme)
gap_guess=gap_guess*meVtoHartree

write (*,'(A24,f12.3,A2)') '  Starting temperature:  ',t_first,' K'
write (*,'(A24,f12.3,A2)') '  Last temperature:      ',t_last,' K'
write (*,'(A33,I3)')      '  Temperature steps:             ',t_number
write (*,'(A24,f12.3,A8)') '  Matsubara en. cutoff: ',wc_cutoff,' Hartree'
write (*,'(A24,f12.3)') '  Coulomb pseudo.(mu*): ',mu_star


! We buid temperature array

t_step=1
if (t_number.ne.1) t_step=((t_last-t_first)/(t_number-1))

allocate (t_array(t_number))
do i=1,t_number
 t_array(i)=(t_first+(i-1)*t_step)*KtoHartree
end do

print *, ''
print *, ' ****************************************'
print *, ' *                                      *'
print *, ' *             READ a2F(w)              *'
print *, ' *                                      *'
print *, ' ****************************************'
print *, ''

!We need to know how many frequencies do we have in real space, and the freq. step

open (unit=11, file=a2f_filename)
counter=0
read (11,*) auxxx
read (11,*) auxxx
do 
  read (11,*,iostat=reason) realnumber, realnumber, realnumber
  if (reason<0) exit
  counter=counter+1
end do
close (11)
wmax_int=counter

allocate (w_a2f_array(wmax_int,2))

open (unit=11, file=a2f_filename)
read (11,*) auxxx
read (11,*) auxxx
do i=1,wmax_int
  read (11,*) w_a2f_array(i,1),w_a2f_array(i,2),realnumber
end do
close (11)

w_a2f_array(:,1)=w_a2f_array(:,1)/2.0d0 !In hartree
w_step= w_a2f_array(2,1)-w_a2f_array(1,1) ! In Hartree


! Start loop on temperature

open (unit=12,file='t_gap.dat')
write (unit=12,fmt='(A)') '# T (K)      Gap (meV)' 

do i = 1, t_number

  write (*,*) ''
  write (*,'(A24,F10.2)') '  Temperature : ', t_array(i) / KtoHartree
  write (*,*) ''

  ! Calculate the number of Matsubara frequencies for this temperature considering the cutoff

!  if (mats_max .lt. 0) mats_max = int((wc_cutoff/(pi*t_array(i)) - 1.)/2.)
  mats_max = int((wc_cutoff/(pi*t_array(i)) - 1.)/2.)
  write (*,'(A28,2I10,A10,I10)') '  Number of Matsubara freq: ', & 
          -mats_max, mats_max,' Total :', 2*mats_max+1

  allocate(matsubara_freq(t_number,2*mats_max+1))

  k = 0
  do j = -mats_max, mats_max
    k = k + 1
    matsubara_freq(i,k) = pi * t_array(i) * dble(2*j+1)
  end do
  
  ! Allocate arrays

  allocate (lambda_array(2*mats_max+1,2*mats_max+1))
  allocate (aux_array(wmax_int)) 
  allocate (delta_array(2*mats_max+1))
  allocate (delta_old(2*mats_max+1))
  allocate (z_array(2*mats_max+1))

  ! Calculate lambda matrix and store it

  call cpu_time(t1)

  lambda_array = 0.0d0

  aux_array(:) = w_a2f_array(:,2)*2.0d0*w_a2f_array(:,1)/(w_a2f_array(:,1)**2.0d0&
                    +(matsubara_freq(i,1)-matsubara_freq(i,1))**2.0d0)
  !      call trapezoid_integral (aux_array(:),w_step,wmax_int,lambda_array(j,k))  
  call trapezoid (aux_array(:), w_a2f_array(:,1), lambda_array(1,1))

  do j = 2, 2*mats_max+1
    lambda_array(j,j) = lambda_array(1,1)
  end do
  
  do j = 1, 2*mats_max+1
    do k = j+1, 2*mats_max+1
      aux_array(:) = w_a2f_array(:,2)*2.0d0*w_a2f_array(:,1)/(w_a2f_array(:,1)**2.0d0&
                +(matsubara_freq(i,j)-matsubara_freq(i,k))**2.0d0)
!      call trapezoid_integral (aux_array(:),w_step,wmax_int,lambda_array(j,k))  
      call trapezoid (aux_array(:), w_a2f_array(:,1), lambda_array(j,k))  
!      print *, ' lambda ', j, k, lambda_array(j,k)
      lambda_array(k,j) = lambda_array(j,k)
    end do
  end do 

  call cpu_time(t2)

  write (*, '(A,F16.3,A)') ' Time needed to calculate lambda matrix = ', t2-t1, ' s.'  

  ! Assign starting gap array

  delta_array = 0.00000001d0
  delta_array(mats_max+1) = gap_guess
  delta_old = delta_array 

  ! Start self-consistent loop

  nit = 0

  call cpu_time(t1)

  do
    nit = nit + 1 
    if (nit .gt. nit_max) then
      write (*,*) ''
      write (*, '(A,I10,A)') ' Convergence not reached after ', nit_max, ' iterations'
      write (*, '(A)') ' Last value of gap written '
      write (*,'(F10.2,F16.6)') t_array(i)  / KtoHartree, delta_array(mats_max+1) / meVtoHartree
      write (*,*) ''
      write (unit=12,fmt='(F7.2,F16.6)') t_array(i)  / KtoHartree, delta_array(mats_max+1) / meVtoHartree              
      exit
    end if
    z_array = 0.0d0
    delta_array = 0.0d0
    do j = 1, 2*mats_max+1
      s = 0.0d0
      do k = 1, 2*mats_max+1
        s = s + lambda_array(j,k) * matsubara_freq(i,k) / &
                     sqrt( matsubara_freq(i,k)**2.0d0 + delta_old(k)**2.0d0 )
      end do
      z_array(j) = 1.0d0 + pi * t_array(i) * s / matsubara_freq(i,j) 
      s = 0.0d0
      do k = 1, 2*mats_max+1
        s = s + (lambda_array(j,k) - mu_star) * delta_old(k) / &
                     sqrt( matsubara_freq(i,k)**2.0d0 + delta_old(k)**2.0d0 )
      end do
      delta_array(j) = pi * t_array(i) * s / z_array(j)     
    end do 
    ! Update gap   
    delta_array = (1.0 - mixing) * delta_array + mixing * delta_old
    ! Check self-consistency condiction
    s = 100. * dot_product(delta_array-delta_old,delta_array-delta_old) / &
        dot_product(delta_old,delta_old)
    ss = abs(delta_array(mats_max+1)-delta_old(mats_max+1))
!    ss = sqrt(dot_product(delta_array,delta_array)) / meVtoHartree
!    if (s .lt. 0.00000001 .or. ss .lt. 0.001) then
    if (s .lt. 0.000000000001 .or. ss .lt. 0.0000001d0*meVtoHartree) then
!    if (abs(delta_old(mats_max+1)-delta_array(mats_max+1))/delta_old(mats_max+1)< 0.00001d0 &
!            .or. abs(delta_array(mats_max+1)-delta_old(mats_max+1))< 0.0001d0*meVtoHartree) then
!    if (abs(delta_array(mats_max+1)-delta_old(mats_max+1))< 0.0001d0*meVtoHartree) then
      write (*,*) ''
      write (*, '(A,I10,A)') ' Convergence reached in ', nit, ' iterations'
      write (*,'(F10.2,F16.6)') t_array(i)  / KtoHartree, delta_array(mats_max+1) / meVtoHartree    
      write (*,*) ''
      write (unit=12,fmt='(F7.2,F16.6)') t_array(i)  / KtoHartree, delta_array(mats_max+1) / meVtoHartree
      exit
    else
      delta_old = delta_array
    end if 
  end do

  call cpu_time(t2)

  write (*, '(A,F16.3,A)') ' Time spent in the self-consistent loop = ', t2-t1, ' s.'  

  ! Update gap for the initial value

  gap_guess = delta_array(mats_max+1) 

  ! If gap is lower than 10^-5 meV we assume it is already closed and we exit the temperature loop

  if (delta_array(mats_max+1) / meVtoHartree .lt. 0.00001) exit

  ! Deallocate staff

  deallocate (matsubara_freq)
  deallocate (lambda_array)
  deallocate (aux_array) 
  deallocate (delta_array)
  deallocate (delta_old)
  deallocate (z_array)

end do

close (unit=12)

contains

subroutine trapezoid_integral (f,h,lth,s)
      integer :: i
      double precision, intent(in) :: h
      integer, intent(in) :: lth
      double precision, dimension(:),intent(in) :: f
      double precision, intent(out) :: s
      s=0.0d0
      do i=1,lth-1
        s = s + (f(i)+f(i+1))/2
      end do
      s=s*h 
end subroutine trapezoid_integral

subroutine trapezoid(f,x,intf)

  double precision, dimension(:), intent(in) :: f, x
  double precision, intent(out) :: intf

  double precision :: a
  integer :: i

  a = 0.50D0 * (f(2) + f(1)) * (x(2) - x(1))

  do i = 2, size(x)-1
    a = a + 0.50D0 * (f(i+1) + f(i)) * (x(i+1) - x(i))
  end do

  intf = a

end subroutine trapezoid

end program me

