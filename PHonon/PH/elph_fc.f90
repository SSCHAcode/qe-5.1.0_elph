! This programs reads the force constants in real space
! and the electron-phonon matrices in q space. This allows
! to combine calculations of the electron-phonon coupling 
! and phonons in different q meshes.
!
! The program needs as input:
!   1) The list of all q points with its corresponding 
!      weight
!   2) File for force constants
!
! A typical input cal look like
!
!&input
!    asr          =   'crystal',  
!    amass(1)     =   195.078
!    amass(2)     =   1.00794
!    flfrc        =   '221.fc', 
!    fildyn       =   'pth.dyn'
!    nbroad       =   2
!    minbroad     =   5
! /
!0.10
!2
!   0.000000000000000E+00   0.000000000000000E+00   0.000000000000000E+00 1
!   0.000000000000000E+00  -0.577350269189585E+00   0.000000000000000E+00 3
!
! After the / we have the gaussian broadening chosen for the
! electron-phonon calculation, the number of q points, ant the
! list of q points with their weight.
!
! IMPORTANT:
! The files where the electron-phonon matrices are 
! stored are assumed to be called
!
! TRIM(fildyn)//TRIM(int_to_char(iq))//'.elph.'//TRIM(int_to_char(iq))
!
! where iq is the q phonon mode. 
!


Module ifconstants
  !
  ! All variables read from file that need dynamical allocation
  !
  USE kinds, ONLY: DP
  REAL(DP), ALLOCATABLE :: frc(:,:,:,:,:,:,:), tau_blk(:,:),  zeu(:,:,:), &
               m_loc(:,:)
  ! frc : interatomic force constants in real space
  ! tau_blk : atomic positions for the original cell
  ! zeu : effective charges for the original cell
  ! m_loc: the magnetic moments of each atom
  INTEGER, ALLOCATABLE  :: ityp_blk(:)
  ! ityp_blk : atomic types for each atom of the original cell
  !
  CHARACTER(LEN=3), ALLOCATABLE :: atm(:)
end Module ifconstants

program elph_fc

  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_bcast
!  USE mp_global,  ONLY : nproc, mpime, mp_startup, mp_global_end
  USE mp_world,   ONLY : world_comm
  USE mp_global,  ONLY : mp_startup, mp_global_end
  USE environment, ONLY : environment_start, environment_end
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_ifc_param, read_ifc
  USE cell_base,  ONLY : at, bg
  USE constants,  ONLY : RY_TO_THZ, RY_TO_CMM1 
  USE symm_base,  ONLY : set_sym
  USE rap_point_group,  ONLY : code_group

  USE ifconstants

implicit none 
  !
  INTEGER :: gid
  !
  ! variables *_blk refer to the original cell, other variables
  ! to the (super)cell (which may coincide with the original cell)
  !
  INTEGER:: nax, nax_blk
  INTEGER, PARAMETER:: ntypx=10, nrwsx=200
  REAL(DP), PARAMETER :: eps=1.0d-6
  INTEGER :: nr1, nr2, nr3, nsc, nk1, nk2, nk3, ntetra, ibrav, total_qpoints
  CHARACTER(LEN=256) :: flfrc, flfrq, flvec, fltau, fldos, filename
  CHARACTER(LEN=256) :: fldyn, fleig, fileq                                   ! << ION ERREA  
  CHARACTER(LEN=10)  :: asr
  CHARACTER(LEN=9)   :: symm_type
  LOGICAL :: dos, has_zstar
  COMPLEX(DP), ALLOCATABLE :: dyn(:,:,:,:), dyn_blk(:,:,:,:)
  double COMPLEX, dimension(:,:), ALLOCATABLE :: z
  double COMPLEX, dimension(:,:,:), ALLOCATABLE :: zq
  double complex, dimension(:,:), allocatable :: dyn_mat, dyn_mat_inv
  double complex, dimension(:,:,:), allocatable :: dyn_matq
  double precision, dimension(:,:), ALLOCATABLE :: tau, q, freq
  double precision, dimension(:), ALLOCATABLE :: w2
  double precision, dimension(:,:), allocatable :: wq, w2q 
  INTEGER, ALLOCATABLE:: tetra(:,:), ityp(:), itau_blk(:)
  REAL(DP) ::     omega,alat, &! cell parameters and volume
                  at_blk(3,3), bg_blk(3,3),  &! original cell
                  omega_blk,                 &! original cell volume
                  epsil(3,3),                &! dielectric tensor
                  amass(ntypx),              &! atomic masses
                  amass_blk(ntypx),          &! original atomic masses
                  atws(3,3),      &! lattice vector for WS initialization
                  rws(0:3,nrwsx)   ! nearest neighbor list, rws(0,*) = norm^2
  !
  INTEGER :: nat, nat_blk, ntyp, ntyp_blk, &
             l1, l2, l3,                   &! supercell dimensions
             nrws,                         &! number of nearest neighbor
             code_group_old

  INTEGER :: nspin_mag, nqs, ios
  !
  LOGICAL :: readtau, la2F, xmlifc, lo_to_split
  !
  REAL(DP) :: qhat(3), qh, DeltaE, Emin=0._dp, Emax, E, DOSofE(1)
  REAL(DP) :: celldm(6), delta, scra
  REAL(DP), ALLOCATABLE :: xqaux(:,:), weightq(:)
  INTEGER, ALLOCATABLE :: nqb(:)
  INTEGER :: n, i, iq, j, l, it, nq, nqx, na, nb, ndos, iout, nqtot, icar, jcar
  INTEGER :: nta, ntb
  INTEGER :: iout_dyn, iout_eig                                          ! << ION ERREA
  LOGICAL, EXTERNAL :: has_xml
  CHARACTER(LEN=15), ALLOCATABLE :: name_rap_mode(:)
  INTEGER, ALLOCATABLE :: num_rap_mode(:,:)
  LOGICAL, ALLOCATABLE :: high_sym(:)
  LOGICAL :: q_in_band_form
  double complex :: complex_number
  double precision :: broad 
  double precision :: rydthz, rydcm1, rydmev, rydk, dosef, lambda_tot, lambda_partial
  character (len=256) :: fileelph, fileconv, fildyn
  double precision, dimension(:,:), allocatable :: lambda, gamma_lw
  double precision, dimension(:), allocatable :: lambdav
  double complex, dimension(:,:,:), allocatable :: lambdav_mat
  double complex, dimension(:,:,:), allocatable :: elphmat
  double precision :: degauss
  CHARACTER(LEN=6) :: int_to_char
  integer :: nbroad
  double precision :: minbroad
  !
  NAMELIST /input/ flfrc, amass, asr, flfrq, flvec, at, dos,  &
       &           fldos, nk1, nk2, nk3, l1, l2, l3, ntyp, readtau, fltau, &
                   la2F, ndos, DeltaE, q_in_band_form, fldyn, fleig, fildyn, &
                   nbroad, minbroad      
  !
  !
  CALL mp_startup()
  CALL environment_start('MATDYN')
  !
  ! Useful constants
  
  rydthz = 13.6058d0*241.796d0
  rydcm1 = 13.6058d0*8065.5d0
  rydmev = 13.6058d0*1000.d0
  rydk   = 13.6058d0*11604.d0

  IF (ionode) CALL input_from_file ( )
     !
     ! ... all calculations are done by the first cpu
     !
     ! set namelist default
     !
     dos = .FALSE.
     deltaE = 1.0d0
     ndos = 1
     nk1 = 0
     nk2 = 0
     nk3 = 0
     asr  ='no'
     readtau=.FALSE.
     flfrc=' '
     fldos='matdyn.dos'
     flfrq='matdyn.freq'
     flvec='matdyn.modes'
     fldyn=' '                                                           ! << ION ERREA
     fleig=' '                                                           ! << ION ERREA
     fltau=' '
     amass(:) =0.d0
     amass_blk(:) =0.d0
     at(:,:) = 0.d0
     ntyp = 0
     l1=1
     l2=1
     l3=1
     la2F=.false.
     q_in_band_form=.FALSE.
     nbroad = 1 ! number of broagdenings for which the a2F(w) will be calculated
     minbroad = 5 ! minimum broadening in cm-1
     !
     !

     !
     !
     IF (ionode) READ (5,input,IOSTAT=ios)
     CALL mp_bcast(ios, ionode_id, world_comm)
     CALL errore('matdyn', 'reading input namelist', ABS(ios))
     CALL mp_bcast(dos,ionode_id, world_comm)
     CALL mp_bcast(deltae,ionode_id, world_comm)
     CALL mp_bcast(ndos,ionode_id, world_comm)
     CALL mp_bcast(nk1,ionode_id, world_comm)
     CALL mp_bcast(nk2,ionode_id, world_comm)
     CALL mp_bcast(nk3,ionode_id, world_comm)
     CALL mp_bcast(asr,ionode_id, world_comm)
     CALL mp_bcast(readtau,ionode_id, world_comm)
     CALL mp_bcast(flfrc,ionode_id, world_comm)
     CALL mp_bcast(fldos,ionode_id, world_comm)
     CALL mp_bcast(flfrq,ionode_id, world_comm)
     CALL mp_bcast(flvec,ionode_id, world_comm)
     CALL mp_bcast(fleig,ionode_id, world_comm)                          ! << ION ERREA
     CALL mp_bcast(fltau,ionode_id, world_comm)
     CALL mp_bcast(amass,ionode_id, world_comm)
     CALL mp_bcast(amass_blk,ionode_id, world_comm)
     CALL mp_bcast(at,ionode_id, world_comm)
     CALL mp_bcast(ntyp,ionode_id, world_comm)
     CALL mp_bcast(l1,ionode_id, world_comm)
     CALL mp_bcast(l2,ionode_id, world_comm)
     CALL mp_bcast(l3,ionode_id, world_comm)
     CALL mp_bcast(la2f,ionode_id, world_comm)
     CALL mp_bcast(q_in_band_form,ionode_id, world_comm)
     !
     ! convert masses to atomic units
     !
     amass(:) = amass(:) * 1.66042d-24/9.1095d-28*0.5d0
     !
     ! read force constants
     !
     ntyp_blk = ntypx ! avoids fake out-of-bound error
     xmlifc=has_xml(flfrc)
     IF (xmlifc) THEN
        CALL read_dyn_mat_param(flfrc,ntyp_blk,nat_blk)
        ALLOCATE (m_loc(3,nat_blk))
        ALLOCATE (tau_blk(3,nat_blk))
        ALLOCATE (ityp_blk(nat_blk))
        ALLOCATE (atm(ntyp_blk))
        ALLOCATE (zeu(3,3,nat_blk))
        CALL read_dyn_mat_header(ntyp_blk, nat_blk, ibrav, nspin_mag, &
                 celldm, at_blk, bg_blk, omega_blk, atm, amass_blk, &
                 tau_blk, ityp_blk,  m_loc, nqs, has_zstar, epsil, zeu )
        alat=celldm(1)
        call volume(alat,at_blk(1,1),at_blk(1,2),at_blk(1,3),omega_blk)
        CALL read_ifc_param(nr1,nr2,nr3)
        ALLOCATE(frc(nr1,nr2,nr3,3,3,nat_blk,nat_blk))
        CALL read_ifc(nr1,nr2,nr3,nat_blk,frc)
     ELSE
        CALL readfc ( flfrc, nr1, nr2, nr3, epsil, nat_blk, &
            ibrav, symm_type, alat, at_blk, ntyp_blk, &
            amass_blk, omega_blk, has_zstar)
     ENDIF
     !
     CALL recips ( at_blk(1,1),at_blk(1,2),at_blk(1,3),  &
          bg_blk(1,1),bg_blk(1,2),bg_blk(1,3) )
     !
     ! set up (super)cell
     !
     if (ntyp < 0) then
        call errore ('matdyn','wrong ntyp ', abs(ntyp))
     else if (ntyp == 0) then
        ntyp=ntyp_blk
     end if
     !
     ! masses (for mass approximation)
     !
     DO it=1,ntyp
        IF (amass(it) < 0.d0) THEN
           CALL errore ('matdyn','wrong mass in the namelist',it)
        ELSE IF (amass(it) == 0.d0) THEN
           IF (it.LE.ntyp_blk) THEN
              WRITE (stdout,'(a,i3,a,a)') ' mass for atomic type ',it,      &
                   &                     ' not given; uses mass from file ',flfrc
              amass(it) = amass_blk(it)
           ELSE
              CALL errore ('matdyn','missing mass in the namelist',it)
           END IF
        END IF
     END DO
     !
     ! lattice vectors
     !
     IF (SUM(ABS(at(:,:))) == 0.d0) THEN
        IF (l1.LE.0 .OR. l2.LE.0 .OR. l3.LE.0) CALL                    &
             &             errore ('matdyn',' wrong l1,l2 or l3',1)
        at(:,1) = at_blk(:,1)*DBLE(l1)
        at(:,2) = at_blk(:,2)*DBLE(l2)
        at(:,3) = at_blk(:,3)*DBLE(l3)
     END IF
     !
     CALL check_at(at,bg_blk,alat,omega)
     !
     ! the supercell contains "nsc" times the original unit cell
     !
     nsc = NINT(omega/omega_blk)
     IF (ABS(omega/omega_blk-nsc) > eps) &
          CALL errore ('matdyn', 'volume ratio not integer', 1)
     !
     ! read/generate atomic positions of the (super)cell
     !
     nat = nat_blk * nsc
     nax = nat
     !!!
     nax_blk = nat_blk
     !!!
     ALLOCATE ( tau (3, nat), ityp(nat), itau_blk(nat) )
     allocate ( dyn_mat (3*nat,3*nat))
     allocate ( dyn_mat_inv (3*nat,3*nat))
     !
     IF (readtau) THEN
        CALL read_tau &
             (nat, nat_blk, ntyp, bg_blk, tau, tau_blk, ityp, itau_blk)
     ELSE
        CALL set_tau  &
             (nat, nat_blk, at, at_blk, tau, tau_blk, ityp, ityp_blk, itau_blk)
     ENDIF
     !
     IF (fltau.NE.' ') CALL write_tau (fltau, nat, tau, ityp)
     !
     ! reciprocal lattice vectors
     !
     CALL recips (at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
     !
     ! build the WS cell corresponding to the force constant grid
     !
     atws(:,1) = at_blk(:,1)*DBLE(nr1)
     atws(:,2) = at_blk(:,2)*DBLE(nr2)
     atws(:,3) = at_blk(:,3)*DBLE(nr3)
     ! initialize WS r-vectors
     CALL wsinit(rws,nrwsx,nrws,atws)
     !
     ! end of (super)cell setup
     !
     IF (dos) THEN
        IF (nk1 < 1 .OR. nk2 < 1 .OR. nk3 < 1) &
             CALL errore  ('matdyn','specify correct q-point grid!',1)
        ntetra = 6 * nk1 * nk2 * nk3
        nqx = nk1*nk2*nk3
        ALLOCATE ( tetra(4,ntetra), q(3,nqx) )
        CALL gen_qpoints (ibrav, at, bg, nat, tau, ityp, nk1, nk2, nk3, &
             symm_type, ntetra, nqx, nq, q, tetra)
     ELSE
        !
        ! read q-point list
        !
        IF (ionode) READ (5,*) degauss
        IF (ionode) READ (5,*) nq
        CALL mp_bcast(nq, ionode_id, world_comm)
        allocate (weightq(nq))
        ALLOCATE ( q(3,nq) )
        ALLOCATE( tetra(1,1) )
        IF (.NOT.q_in_band_form) THEN
           DO n = 1,nq
              IF (ionode) READ (5,*) (q(i,n),i=1,3), weightq(n)
           END DO
           CALL mp_bcast(q, ionode_id, world_comm)
        ELSE
           ALLOCATE(nqb(nq))
           ALLOCATE(xqaux(3,nq))
           DO n = 1,nq
              IF (ionode) READ (5,*) (q(i,n),i=1,3), weightq(n)
           END DO
           CALL mp_bcast(q, ionode_id, world_comm)
           CALL mp_bcast(nqb, ionode_id, world_comm)
           nqtot=SUM(nqb(1:nq-1))+1
           xqaux(:,1:nq)=q(:,1:nq)
           DEALLOCATE(q)
           ALLOCATE(q(3,nqtot))
           nqtot=0
           DO i=1,nq-1
                 nqtot=nqtot+1
                 q(:,nqtot) = xqaux(:,i)
           ENDDO
           nqtot=nqtot+1
           q(:,nqtot)=xqaux(:,nq)
           nq=nqtot
           DEALLOCATE(xqaux)
           DEALLOCATE(nqb)
        END IF
     END IF
     !
     IF (asr /= 'no') THEN
        CALL set_asr (asr, nr1, nr2, nr3, frc, zeu, &
             nat_blk, ibrav, tau_blk)
     END IF
     !
     IF (flvec.EQ.' ') THEN
        iout=6
     ELSE
        iout=4
        IF (ionode) OPEN (unit=iout,file=flvec,status='unknown',form='formatted')
     END IF
     !
     IF (fldyn.EQ.' ') THEN                                              ! << ION ERREA
        iout_dyn=6                                                       ! << ION ERREA
     ELSE                                                                ! << ION ERREA
        iout_dyn=44                                                      ! << ION ERREA
        OPEN (unit=iout_dyn,file=fldyn,status='unknown',form='formatted')! << ION ERREA
     END IF                                                              ! << ION ERREA 
     !
     IF (fleig.EQ.' ') THEN                                                          ! << ION ERRE
        iout_eig=0                                                                   ! << ION ERRE
     ELSE                                                                            ! << ION ERRE
        iout_eig=313                                                                 ! << ION ERRE
        IF (ionode) OPEN (unit=iout_eig,file=fleig,status='unknown',form='formatted')! << ION ERRE
     END IF                                                                          ! << ION ERRE

     ALLOCATE ( dyn(3,3,nat,nat), dyn_blk(3,3,nat_blk,nat_blk) )
     ALLOCATE ( zq(nq,3*nat,3*nat), w2(3*nat), wq(nq,3*nat) )
     ALLOCATE ( z(3*nat,3*nat), w2q(nq,3*nat) )
     allocate ( dyn_matq (nq,3*nat,3*nat))
     allocate(gamma_lw(nq,3*nat))
     allocate(elphmat(nq,3*nat,3*nat))
     allocate(lambda(nq,3*nat))
     allocate(lambdav(nq))
     allocate(lambdav_mat(nq,3*nat,3*nat))

     if(la2F.and.ionode) open(300,file='dyna2F',status='unknown')
     IF (xmlifc) CALL set_sym(nat, tau, ityp, nspin_mag, m_loc, 6, 6, 6 )

     ALLOCATE(num_rap_mode(3*nat,nq))
     ALLOCATE(high_sym(nq))
     num_rap_mode=-1
     high_sym=.TRUE.

! After having read all the data start loop on q points

     DO n=1, nq

        print *,''
        print *,'***************************************'
        print *,'* Q POINT    ', n
        print *,'*'
        print *,'*       q -> ', q(:,n)
        print *,'*  weight -> ', weightq(n)
        print *,'***************************************'
        print *,''


        dyn(:,:,:,:) = (0.d0, 0.d0)

        lo_to_split=.FALSE.
        CALL setupmat (q(1,n), dyn, nat, at, bg, tau, itau_blk, nsc, alat, &
             dyn_blk, nat_blk, at_blk, bg_blk, tau_blk, omega_blk,  &
             epsil, zeu, frc, nr1,nr2,nr3, has_zstar, rws, nrws)

        IF (q(1,n)==0.d0 .AND. q(2,n)==0.d0 .AND. q(3,n)==0.d0) THEN
           !
           ! q = 0 : we need the direction q => 0 for the non-analytic part
           !
           IF ( n == 1 ) THEN
              ! if q is the first point in the list
              IF ( nq > 1 ) THEN
                 ! one more point
                 qhat(:) = q(:,n) - q(:,n+1)
              ELSE
                 ! no more points
                 qhat(:) = 0.d0
              END IF
           ELSE IF ( n > 1 ) THEN
              ! if q is not the first point in the list
              IF ( q(1,n-1)==0.d0 .AND. &
                   q(2,n-1)==0.d0 .AND. &
                   q(3,n-1)==0.d0 .AND. n < nq ) THEN
                 ! if the preceding q is also 0 :
                 qhat(:) = q(:,n) - q(:,n+1)
              ELSE
                 ! if the preceding q is npt 0 :
                 qhat(:) = q(:,n) - q(:,n-1)
              END IF
           END IF
           qh = SQRT(qhat(1)**2+qhat(2)**2+qhat(3)**2)
           ! write(*,*) ' qh,  has_zstar ',qh,  has_zstar
           IF (qh /= 0.d0) qhat(:) = qhat(:) / qh
           IF (qh /= 0.d0 .AND. .NOT. has_zstar) THEN
                CALL infomsg  &
                ('matdyn','Z* not found in file '//TRIM(flfrc)// &
                          ', TO-LO splitting at q=0 will be absent!')
           ELSE
              lo_to_split=.TRUE.
           ENDIF
           !
           CALL nonanal (nat, nat_blk, itau_blk, epsil, qhat, zeu, omega, dyn)
           !
        END IF

        ! Print dynamical matrix

        write (unit=*, fmt='(a)') ''
        write (unit=*, fmt='(a)') ' Printing dynamical matrix:'
        write (unit=*, fmt='(a)') ''
        do i = 1, nat
          do j = 1, nat
            write (unit=*, fmt='(2i5)') i, j
            do icar = 1, 3
              write (unit=*, fmt='(3(2f12.8,2x))') (dyn(icar,jcar,i,j), jcar=1,3)
            end do
          end do
        end do        

        !
        ! Write the dynamical matrix as a matrix and divide by the masses        
    
        do i=1,3*nat
          na= ( i-1 ) / 3 + 1
          icar=i-3*(na-1)
          do j=1,3*nat
            nb= ( j-1 ) / 3 + 1
            jcar=j-3*(nb-1)
            dyn_mat(i,j)=dyn(icar,jcar,na,nb) / sqrt(amass(ityp(na))*amass(ityp(nb)))
          enddo
        enddo
        
!        ! Divide by the masses
!        do i = 1,3*nat
!          na = (i-1)/3+1
!          nta = ityp(na)
!          do j = 1,3*nat
!            nb = (j-1)/3+1
!            ntb=ityp(nb)
!            dyn_mat(i,j)=dyn_mat(i,j)/sqrt(amass(nta)*amass(ntb))
!          end do
!        end do
    
        ! Check hermiticity and impose it
        do i = 1, 3*nat
          do j = 2, 3*nat
            scra = abs(dyn_mat(i,j) - conjg(dyn_mat(j,i)))
            if (scra .gt. 1.0d-8) then
              print *, '  WARNING: dynamical matrix non hermitian'
              print *, i, j, '->', dyn_mat(i,j)
              print *, j, i, '->', dyn_mat(j,i)
!              dyn_mat(j,i) = dyn_mat(i,j) 
!            else
!              dyn_mat(j,i) = dyn_mat(i,j) 
            end if
          end do
        end do        

        ! Diagonalize dynamical matrix

!        call diag_dynmat(dyn_mat,w2,z) 
        call cdiagh2(3*nat,dyn_mat,3*nat,w2,z)

        dyn_matq(n,:,:) = dyn_mat(:,:)  
        zq(n,:,:) = z(:,:)
        w2q(n,:) = w2(:)
          
        ! Print eigenvalues and eigenvectors
  
        print *, ''
        do i = 1, 3*nat
          if (w2(i) < 0) then
            wq(n,i) = - sqrt(-w2(i))
          else
            wq(n,i) = sqrt(w2(i))
          end if
        end do
        do i = 1, 3*nat
          print '(a,i3,a,f12.6)', '  w (cm-1)   ', i, ' = ', wq(n,i)*rydcm1
          do j = 1, nat
            print '(a,i3,i3,a,3(2f12.8,2x))', '  Eigenvector ', i, j,' = ', zq(n,3*(j-1)+1:3*j,i) 
          end do
          do j = 1, nat
            print '(a,i3,i3,a,3(2f12.8,2x))', '  Eigenvector / sqrt(M)', i, j,' = ', zq(n,3*(j-1)+1:3*j,i) / sqrt(amass(ityp(j))) 
          end do
        end do
        print *, ''

        ! Get the electron phonon coefficients

!        if ( n .lt. 10) then
!          write (fileelph, '(a,i1)') 'elph.d.mat.', n 
!        else if ( n .lt. 100) then
!          write (fileelph, '(a,i2)') 'elph.d.mat.', n 
!        else if ( n .lt. 1000) then
!          write (fileelph, '(a,i3)') 'elph.d.mat.', n 
!        end if
        fileelph=TRIM(fildyn)//TRIM(int_to_char(n))//'.elph.d.mat.'//TRIM(int_to_char(n))

        call read_elph(fileelph,degauss,q(:,n),w2q(n,:),zq(n,:,:),ityp,amass,gamma_lw(n,:), &
                       lambda(n,:),lambdav(n),elphmat(n,:,:),dosef)

        ! Print the electron-phonon data

        print *, ''
        do i = 1, nat*3
          print '(a,i3,a,f12.6,a,f12.6,a,f12.6)', ' Mode: ', i, ' Omega (cm-1) = ', wq(n,i)*rydcm1, &
                              ' Gamma (GHz) = ', gamma_lw(n,i)*rydthz*1000.0d0, ' lambda = ', lambda(n,i)
        end do

        print *, ''
        write(*,*) '  Sum of lambda_nu =', lambdav(n)

        ! Multiply the dynamical matrix by the masses in order to get 
        ! the force constants

        do i = 1,3*nat
          na = (i-1)/3+1
          nta = ityp(na)
          do j = 1,3*nat
            nb = (j-1)/3+1
            ntb=ityp(nb)
            dyn_mat(i,j)=dyn_mat(i,j)*sqrt(amass(nta)*amass(ntb))
          end do
        end do

        ! If this is the gamma point skip it since the determinant is equal to and
        ! the dynamical matrix has no inverse   
        if (q(1,n).eq.0.0d0 .and. q(2,n).eq.0.0d0 .and. q(3,n).eq.0.0d0) then
          print *, ' '     
          print *, ' It is the Gamma point: SKIPPING IT, electron-phonon contribution skipped!'     
          print *, ' '     
          cycle
        end if 

        ! Invert the dynamical matrix

        call invert_dynmat(dyn_mat,dyn_mat_inv)

        ! Calculate the value of lambda from the product
        ! of the matrices in coordinate space

        complex_number = 0.0d0
        do i = 1, 3*nat
          do j = 1, 3*nat
            lambdav_mat(n,i,j) = conjg(dyn_mat_inv(i,j))*elphmat(n,i,j) &
                                /2.0d0/dosef    
            complex_number = complex_number + lambdav_mat(n,i,j)  
          end do
        end do
        lambdav(n) = real(complex_number)

        write(*,*) ' Calculation with matrix product:'
        write(*,*) '  Sum of lambda_nu =', lambdav(n)  
   
     ENDDO

     ! Calculate lambda

     lambda_tot = 0.0d0
     total_qpoints = 0
 
     do i = 1, nq
       total_qpoints = total_qpoints + weightq(i)
       if (q(1,i).eq.0.0d0 .and. q(2,i).eq.0.0d0 .and. q(3,i).eq.0.0d0) then 
         total_qpoints = total_qpoints - 1  
         cycle
       end if
       lambda_tot = lambda_tot + lambdav(i)*weightq(i) 
     end do 

     print *, ''
     print *, '  Total number of q points -> ', total_qpoints

     print *, ''
     print *, '  lambda = ', lambda_tot / dble(total_qpoints)

     ! Calculate partial contributions of atoms to lambda using the matrix

     print *, ''
     print *, '  Partial contributions of atoms from matrix:'

     lambda_tot = 0.0d0
     do i = 1, nat 
       do j = 1, nat
         complex_number = 0.0d0
         do iq = 1, nq  
           if (q(1,iq).eq.0.0d0 .and. q(2,iq).eq.0.0d0 .and. q(3,iq).eq.0.0d0) cycle
           do icar = 1, 3
             do jcar = 1, 3     
               complex_number = complex_number + lambdav_mat(iq,(i-1)*3+icar,(j-1)*3+jcar) &
                                * weightq(iq) / total_qpoints 
             end do
           end do 
         end do 
         lambda_partial = real(complex_number)
         write(*,*) '    lambda contribution', i, j, ' = ', lambda_partial 
         lambda_tot = lambda_tot + lambda_partial  
       end do
     end do       

     lambda_tot = lambda_tot 
     
     print *, ''
     print *, '    Sum of partial contributions in matrix = ', lambda_tot

     ! Calculate the a2F(w), PDOS, lambda integral, partial contributions
     ! to a2F(w) and the value of lambda obtained from the Eliashberg
     ! function for different gaussian broadenings for the phonon delta.
     ! The broadening goes from 1 cm^-1 (files with 1)to 10 cm^-1 
     ! (files with 10) 

     do i = 1, nbroad

       print *, ''
       print *,'***************************************'
       print '(a,f6.3,a)','* Broadening ', dble(i) * minbroad, ' [cm-1]' 
       print *,'*'
       print *,'*   writing output files...'
       print *,'*  '
       print *,'***************************************'
       print *, ''

       broad = dble(i) * minbroad / rydcm1
 
       call eliashberg(wq,lambda,broad,q,weightq,zq,ityp,amass,elphmat,dosef)

     end do
 
     print *, ''  
     print *, ' DONE'  
     print *, ''  
     
!--------------------------------
!
! Here we include the subroutines 
!
!--------------------------------

contains

  ! Subroutine that calculates the Eliashberg function,
  ! its decomposition into different atoms, Tc using McMillan
  ! equation and the phonon density of states (PDOS). The Eliashberg
  ! function and the PDOS are calculated without any interpolation
  ! and using gaussians for the delta function with respect to the
  ! phonon frequencies.

  subroutine eliashberg(w,lambda,broad,q,weight,z,ityp,amass,elphmat,dosef)

    implicit none

    double precision, dimension(:,:), intent(in) :: w, lambda, q
    double precision, intent(in) :: broad, dosef
    double complex, dimension(:,:,:), intent(in) :: z, elphmat
    double precision, dimension(:), intent(in) :: amass, weight
    integer, dimension(:), intent(in) :: ityp
    
    double precision :: minfreq, maxfreq, w_aux, gaussian, twopi, rydcm1, rydmev, rydk 
    double precision :: lambda_tot, omega_log, mustar, pdos_tot, f1, f2, lambda1, lambda2, omega2, omega2bar
    double complex :: lmat
    double precision, dimension(10000) :: a2F, w_vector, lambda_vector, omega_log_vector, omega2_vector
    double precision, dimension(10000) :: pdos, lambda_int
    double complex, dimension(:,:,:), allocatable :: e
    double precision, dimension(:,:,:), allocatable :: a2Fpartial
    double precision, dimension(:,:), allocatable :: pdospartial
    integer :: i, j, l, na, total_qpoints, total_qpoints_a2f, iw, iq, inu, alpha, beta
    integer :: nat, nq, intbroad
    character (len=50) :: file_a2F, file_pdos, file_lambda

    ! Useful constants

    rydcm1 = 13.6058d0*8065.5d0
    rydmev = 13.6058d0*1000.d0
    rydk   = 13.6058d0*11604.d0
    twopi  = 6.283185307179586d0

    nat = size(w(1,:)) / 3
    nq  = size(w(:,1))

    intbroad = int(broad*rydcm1)

    allocate(a2Fpartial(10000,3*nat,3*nat))
    allocate(pdospartial(10000,nat))
    allocate(e(nq,3*nat,3*nat))
    e = z

    ! Total number of q points

    total_qpoints = 0

    do i = 1, nq
      total_qpoints = total_qpoints + weightq(i)
    end do    

    total_qpoints_a2f = total_qpoints 

    do i = 1, nq
      if (q(1,i).eq.0.0d0 .and. q(2,i).eq.0.0d0 .and. q(3,i).eq.0.0d0) then
        total_qpoints_a2f = total_qpoints - 1   
      end if
    end do

    ! Calculate the a2F(w) 

    minfreq = minval(w(:,:))
    maxfreq = maxval(w(:,:)) * 1.2d0

    file_a2F = 'a2F.'//trim(int_to_char(intbroad))//'.dat'
    file_lambda = 'lambda_int.'//trim(int_to_char(intbroad))//'.dat'
    file_pdos = 'pdos..'//trim(int_to_char(intbroad))//'.dat'
    
!    if (intbroad .lt. 10) then
!      write (file_a2F,'(a,i1,a)') 'a2F.', intbroad, '.dat'
!      write (file_lambda,'(a,i1,a)') 'lambda_int.', intbroad, '.dat'
!      write (file_pdos,'(a,i1,a)') 'pdos.', intbroad, '.dat'
!    else
!      write (file_a2F,'(a,i2,a)') 'a2F.', intbroad, '.dat'
!      write (file_lambda,'(a,i2,a)') 'lambda_int.', intbroad, '.dat'
!      write (file_pdos,'(a,i2,a)') 'pdos.', intbroad, '.dat'
!    end if
   
    open (unit=1,file=file_a2F)
    open (unit=2,file=file_pdos)

    write (unit=1, fmt=*) '# Broadening: ', intbroad, ' [cm-1]' 
    write (unit=1, fmt=*) '#    w [Ry]          a2F' 
    write (unit=2, fmt=*) '# Broadening: ', intbroad, ' [cm-1]' 
    write (unit=2, fmt=*) '#    w [Ry]          PDOS' 

    do i = 1, 10000
      w_aux = minfreq + dble(i-1)*(maxfreq-minfreq)/9999.0d0  
      a2F(i) = 0.0d0
      pdos(i) = 0.0d0
      w_vector(i) = w_aux
      do j = 1, nq
        do l = 1, 3*nat
          gaussian = exp(-(w_aux-w(j,l))**2.0d0/(2.0d0*broad**2.0d0)) &
                     / (broad*sqrt(twopi))
          pdos(i) = pdos(i) + gaussian * weightq(j) / dble(total_qpoints)
          if (q(1,j).eq.0.0d0 .and. q(2,j).eq.0.0d0 .and. q(3,j).eq.0.0d0) then
            cycle
          else     
            a2F(i) = a2F(i) + gaussian * lambda(j,l) * w(j,l) * weightq(j) &
                       / (2.0d0 * dble(total_qpoints_a2f))
          end if
        end do
      end do            
      lambda_vector(i) = 2.0d0 * a2F(i) / w_vector(i) 
      omega_log_vector(i) = log( w_vector(i) ) * a2F(i) / w_vector(i)   
      omega2_vector(i) =  a2F(i) * w_vector(i)
      write (unit=1, fmt=*) w_vector(i) , a2F(i), lambda_vector(i)  
      write (unit=2, fmt=*) w_vector(i) , pdos(i)
    end do

    close (unit=1)
    close (unit=2)

    ! Calculate lambda integrating the Eliashberg function

    na = 1   
    do i = 1, 10000
      if (w_vector(i) .lt. 0.0d0) then
        na = na + 1
      else
        exit
      end if          
    end do

    call trapezoid_int(lambda_vector(na:10000),w_vector(na:10000),lambda_int(na:10000),lambda_tot)
    call trapezoid(omega_log_vector(na:10000),w_vector(na:10000),omega_log)
    call trapezoid(pdos,w_vector,pdos_tot)

    ! Calculate other parameters for modified Allen-Dynes

    call trapezoid(omega2_vector(na:10000),w_vector(na:10000),omega2)

    omega2 = 2. * omega2 / lambda_tot
    omega2bar = omega2 ** 0.5

    print *, '  Integral of PDOS -> ', pdos_tot

    open (unit=9,file=file_lambda)

    write (unit=9, fmt=*) '# Broadening: ', intbroad, ' [cm-1]' 
    write (unit=9, fmt=*) '#    w [Ry]          lambda (w) '

    do i = na, 10000
      write (unit=9,fmt=*) w_vector(i), lambda_int(i)
    end do

    close (unit=9)

    print *, ''
    print *, '  Calculating lambda and w_log from Eliashberg function:'
    print *, '     lambda = ', lambda_tot 

    omega_log = exp(2.0d0*omega_log/lambda_tot)

    print *, '     w_log  = ', omega_log * rydcm1, ' [cm-1] ', omega_log * rydmev, ' [meV] '  
  
    print *, ''
    print *, '  Calculating Tc using McMillan equation: '
    do i = 1, 21
      mustar = 0.07d0 + dble(i-1) / (200.0d0) 
      print '(a,f6.3,a,f8.3)', '     mustar -> ', mustar, ', Tc [K] = ', omega_log * rydk * &
                       exp(-(1.04d0*(1.d0+lambda_tot))/(lambda_tot-mustar*(1.d0+0.62d0*lambda_tot))) &
                       / 1.2d0
    end do

    print *, ''
    print *, '  Calculating Tc using McMillan Allen-Dynes modified equation: '
    do i = 1, 21
      mustar = 0.07d0 + dble(i-1) / (200.0d0) 
      lambda1 = 2.46*(1+3.8*mustar)
      lambda2 = 1.82*(1+6.3*mustar)*omega2bar/omega_log
      f1 = (1.+(lambda_tot/lambda1)**1.5)**(0.33333333333)
      f2 = 1. + (omega2bar/omega_log - 1.) * lambda_tot**2. / (lambda_tot**2. + lambda2**2.)
      print '(a,f6.3,a,f8.3)', '     mustar -> ', mustar, ', Tc [K] = ', omega_log * rydk * f1 * f2 * &
                       exp(-(1.04d0*(1.d0+lambda_tot))/(lambda_tot-mustar*(1.d0+0.62d0*lambda_tot))) &
                       / 1.2d0
    end do

    ! Calculate partial contributions to Eliashberg function

    ! Divide polarization vectors by the mass
    do l = 1, nq
      do j = 1, 3 * nat
        do i = 1, 3 * nat
          na = (i - 1) / 3 + 1
          e (l,i, j) = e (l,i, j) / sqrt (amass (ityp (na) ) )
        enddo
      enddo      
    enddo

    ! Start loop to build the partial Eliashberg functions
    do i = 1, 3*nat
      do j = 1, 3*nat
        do iw = 1, 10000
          a2Fpartial(iw,i,j) = 0.0d0
          do iq = 1, nq 
            if (q(1,iq).eq.0.0d0 .and. q(2,iq).eq.0.0d0 .and. q(3,iq).eq.0.0d0) cycle
            lmat = (0.0d0,0.0d0)
            do inu = 1, 3*nat
              gaussian =  exp(-(w_vector(iw)-w(iq,inu))**2.0d0/(2.0d0*broad**2.0d0)) &
                          / (broad*sqrt(twopi))
              lmat = lmat + gaussian * conjg(e(iq,i,inu)) * e(iq,j,inu) / w(iq,inu)
            end do
            a2Fpartial(iw,i,j) = a2Fpartial(iw,i,j) + real(lmat * elphmat(iq,i,j)) &
                        * weight(iq) / (4.0d0*dosef*dble(total_qpoints))  
          end do            
        end do
      end do
    end do 

    ! Write partial contribution on file
    do i = 1, nat
      do j = 1, nat
        file_a2F ='a2Fpartial_'//trim(int_to_char(i))//'_'//trim(int_to_char(j))//'.'//trim(int_to_char(intbroad))//'.dat' 
!        if (intbroad .lt. 10) then
!          write (file_a2F, fmt='(a,i1,a,i1,a,i1,a)') 'a2Fpartial_', i, '_', j, '.', intbroad, '.dat'       
!        else
!          write (file_a2F, fmt='(a,i1,a,i1,a,i2,a)') 'a2Fpartial_', i, '_', j, '.', intbroad, '.dat'       
!        end if
        open (unit=1, file=file_a2F)
        write (unit=1, fmt=*) '# Broadening: ', intbroad, ' [cm-1]' 
        write (unit=1, fmt=*) '#    w [Ry]          a2F ' 
        a2F = 0.0d0
        do iw = 1, 10000 
          do alpha = 1, 3
            do beta = 1, 3
              a2F(iw) = a2F(iw) + a2Fpartial(iw,(i-1)*3+alpha,(j-1)*3+beta)
            end do
          end do
          write (unit=1,fmt=*)  w_vector(iw) , a2F(iw)   
        end do
        close (unit=1) 
      end do
    end do 

    ! Restore the polarization vectors 
    do l = 1, nq
      do j = 1, 3 * nat
        do i = 1, 3 * nat
          na = (i - 1) / 3 + 1
          e (l,i, j) = e (l,i, j) * sqrt (amass (ityp (na) ) )
        enddo
      enddo      
    enddo

    ! Calculate partial contributions to the PDOS 

    ! Start loop to build the partial PDOS                   
    do i = 1, nat
      do iw = 1, 10000
        pdospartial(iw,i) = 0.0d0
        do iq = 1, nq 
          lmat = (0.0d0,0.0d0)
          do inu = 1, 3*nat
            gaussian =  exp(-(w_vector(iw)-w(iq,inu))**2.0d0/(2.0d0*broad**2.0d0)) &
                        / (broad*sqrt(twopi))
            do alpha = 1, 3 
              lmat = lmat + gaussian * conjg(e(iq,3*(i-1)+alpha,inu)) * &
                            e(iq,3*(i-1)+alpha,inu) 
            end do
          end do
          pdospartial(iw,i) = pdospartial(iw,i) + real(lmat) &
                      * weight(iq) / (dble(total_qpoints))  
        end do
      end do
    end do 

    ! Write partial contribution on file
    do i = 1, nat
      file_a2F ='pdospartial_'//trim(int_to_char(i))//'.'//trim(int_to_char(intbroad))//'.dat'
!      if (intbroad .lt. 10) then 
!        write (file_a2F, fmt='(a,i1,a,i1,a)') 'pdospartial_', i, '.', intbroad, '.dat'     
!      else 
!        write (file_a2F, fmt='(a,i1,a,i2,a)') 'pdospartial_', i, '.', intbroad, '.dat'     
!      end if
      open (unit=1, file=file_a2F)
      write (unit=1, fmt=*) '# Broadening: ', intbroad, ' [cm-1]' 
      write (unit=1, fmt=*) '#    w [Ry]          PDOS ' 
      do iw = 1, 10000 
        write (unit=1,fmt=*)  w_vector(iw) , pdospartial(iw,i)   
      end do
      close (unit=1) 
    end do 

    deallocate(a2Fpartial)
    deallocate(pdospartial)
    deallocate(e)
 
  end subroutine eliashberg

  !This subroutine uses the trapezoidal rule for integration

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

  ! This subroutine uses the trapezoidal rule for integration
  ! but saves each step of the integral in a vector 

  subroutine trapezoid_int(f,x,v_int,intf)

    double precision, dimension(:), intent(in) :: f, x
    double precision, dimension(:), intent(out) :: v_int
    double precision, intent(out) :: intf

    double precision :: a
    integer :: i
 
    v_int(1) = 0.0d0

    a = 0.50D0 * (f(2) + f(1)) * (x(2) - x(1))

    v_int(2) = a

    do i = 2, size(x)-1
      a = a + 0.50D0 * (f(i+1) + f(i)) * (x(i+1) - x(i))
      v_int (i+1) = a
    end do

    intf = a

  end subroutine trapezoid_int

  ! Subroutine that inverts the dynamical matrix

  subroutine invert_dynmat(dyn,inv_dyn)

    double complex, dimension(:,:), intent(in) :: dyn
    double complex, dimension(:,:), intent(out) :: inv_dyn

    integer :: nat, info
    double precision :: scra
    integer, dimension(:), allocatable :: ipvt
    double precision, dimension(:), allocatable :: work
    double complex, dimension(:,:), allocatable :: check_inv

    nat = size(dyn(1,:)) / 3

    allocate(ipvt(3*nat))
    allocate(work(64*3*nat))
    allocate(check_inv(3*nat,3*nat))

    inv_dyn = dyn
    
    call zgetrf(3*nat,3*nat,inv_dyn,3*nat,ipvt,info)
    
    if(info.ne.0) then
       write(*,*) 'Error in zgetrf, info=',info
       stop
    endif

    call zgetri(3*nat,inv_dyn,3*nat,ipvt,work,size(work),info)

    if(info.ne.0) then
       write(*,*) 'Error in zgetri, info=',info
       stop
    endif 

    na = 0
    check_inv = matmul(dyn,inv_dyn)
    do i = 1, 3*nat
      do j = 1, 3*nat
        if (i.eq.j) then
          scra = abs(check_inv(i,j)) -1.0d0
        else
          scra = abs(check_inv(i,j))
        end if
        if (scra.lt.1.0d-8) na = na + 1
      end do
    end do
    print *, ''
    if (na.eq.9*nat**2) then
      print *, '  Inversion of matrix is correct! '
    else
      print *, '  Inversion of matrix is not correct! '
      print *, '  STOPPING...'
      stop
    end if
    print *, ''
    
    deallocate(work)
    deallocate(ipvt)
    deallocate(check_inv)
   
  end subroutine invert_dynmat


  ! Subroutine that diagonalizes the dynamical matrix
  ! and gives as an output the polarization vectors and
  ! the eigenvalues

  subroutine diag_dynmat(dynamical,w2,z)

    double complex, dimension(:,:), intent(in) :: dynamical
    double precision, dimension(:), intent(out) :: w2
    double complex, dimension(:,:), intent(out) :: z

    integer :: nat, info, lwork
    double complex, dimension(:,:), allocatable :: dyn_aux
    double precision, dimension(:), allocatable :: work, rwork

    nat = size(dynamical(1,:)) / 3

    lwork = 2*3*nat - 1

    allocate(dyn_aux(3*nat,3*nat))
    allocate(work(lwork))
    allocate(rwork(3*3*nat-2))

    dyn_aux = dynamical

    call zheev ('V','U',3*nat,dyn_aux,3*nat,w2,work,&
                lwork,rwork,info)
    
    z = dyn_aux  

    deallocate(dyn_aux)
    deallocate(work)
    deallocate(rwork)
   
  end subroutine diag_dynmat
 
  ! This subroutine reads the electron-phonon matrix elements
  ! and calculates the linewidth and lambda of each mode and 
  ! the sum of lambda with respect to the modes. Moreover, it
  ! also gives the electron-phonon matrices.

  subroutine read_elph(fileelph,degauss,q,w2,e,ityp,amass, &
                       gamma_lw,lambda,lambdav,elphmat,dosef)

  implicit none

  character (len=256), intent(in) :: fileelph 
  double precision, intent(in) :: degauss
  double precision, dimension(:), intent(in) :: w2, amass, q
  double complex, dimension(:,:), intent(in) :: e
  integer, dimension(:), intent(in) :: ityp
  double precision, dimension(:), intent(out) :: gamma_lw, lambda
  double precision, intent(out) :: lambdav
  double complex, dimension(:,:), intent(out) :: elphmat
  double precision, intent(out) :: dosef

  double complex, dimension(:,:), allocatable :: z   
  double precision, dimension(:), allocatable :: w, w2_, lambda_, gamma_lw_      
  double precision, dimension(3) :: qr
  double precision :: ef, degauss_
  double precision :: scra
  double complex :: complex_number
  integer :: nat, i, j, k, l, na, ngauss, nbr, nbr_, nsig

  nat = size(w2) / 3
  nbr = nat*3

  nsig = 30
  
  allocate(z(3*nat,3*nat))
  allocate(gamma_lw_(3*nat))
  allocate(lambda_(3*nat))
  allocate(w(3*nat))
  allocate(w2_(3*nat))

  do i = 1, 3*nat
    if (w2(i) < 0) then
      w(i) = - sqrt(-w2(i))
    else
      w(i) = sqrt(w2(i))
    end if
  end do

  z = e

  ! Divide polarization vectors by the mass
  do j = 1, 3 * nat
    do i = 1, 3 * nat
      na = (i - 1) / 3 + 1
      z (i, j) = z (i, j) / sqrt (amass (ityp (na) ) )
    enddo
  enddo

  ! Read the electron-phonon matrix from file
  open(unit=3,file=fileelph)

  do k=1,nsig
     read(unit=3,fmt=*) degauss_     
     read(unit=3,fmt=*) dosef        
     do i=1,3*nat
        do j=1,3*nat
           read(unit=3,fmt=*) elphmat(i,j)
        enddo
     enddo
     scra=dabs(degauss-degauss_)
     if(scra.lt.1.d-8) then     

       write (*,*) '  '
       write (*,*) '  Gaussian broadening coincides', degauss

       ! I calculate the electron-phonon coupling to see if
       !  it matches with the preceiding case.

       do j = 1, nbr
         complex_number = (0.0d0,0.0d0)
         do i = 1, 3 * nat
           do l = 1, 3 * nat
             complex_number = complex_number + conjg (z (i, j) ) &
                          * elphmat (i, l) * z (l, j) 
           enddo
         enddo
         print *,  complex_number
         gamma_lw(j) = 3.1415926 * real( complex_number ) / 2.d0
         lambda(j) = gamma_lw(j) / 3.1415926 / w2(j) / dosef
       enddo

       ! Calculate sum of lambda with respect to modes
       lambdav=0.d0
       do j=1,3*nat
         lambdav=lambdav+lambda(j)
       enddo
 
       exit

     end if
  enddo

  close(3)

  deallocate(z)
  deallocate(w)
  deallocate(w2_)
  deallocate(gamma_lw_)
  deallocate(lambda_)
  
  end subroutine read_elph
  

!
!-----------------------------------------------------------------------
SUBROUTINE readfc ( flfrc, nr1, nr2, nr3, epsil, nat,    &
                    ibrav, symm_type, alat, at, ntyp, amass, omega, has_zstar )
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE ifconstants,ONLY : tau => tau_blk, ityp => ityp_blk, frc, zeu
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE mp,         ONLY : mp_bcast 
  USE mp_world,   ONLY : world_comm
  !
  IMPLICIT NONE
  ! I/O variable
  CHARACTER(LEN=256) flfrc
  INTEGER ibrav, nr1,nr2,nr3,nat, ntyp
  REAL(DP) alat, at(3,3), epsil(3,3)
  LOGICAL has_zstar
  ! local variables
  INTEGER i, j, na, nb, m1,m2,m3
  INTEGER ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
  REAL(DP) amass(ntyp), amass_from_file, celldm(6), omega
  INTEGER nt
  CHARACTER(LEN=3) atm
  CHARACTER(LEN=9) symm_type
  !
  !
  IF (ionode) OPEN (unit=1,file=flfrc,status='old',form='formatted')
  !
  !  read cell data
  !
  IF (ionode)THEN
     READ(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
     if (ibrav==0) then
!        read(1,'(a)') symm_type
        read(1,*) ((at(i,j),i=1,3),j=1,3)
     end if
  ENDIF
  CALL mp_bcast(ntyp, ionode_id, world_comm)
  CALL mp_bcast(nat, ionode_id, world_comm)
  CALL mp_bcast(ibrav, ionode_id, world_comm)
  CALL mp_bcast(celldm, ionode_id, world_comm)
  IF (ibrav==0) THEN
     CALL mp_bcast(symm_type, ionode_id, world_comm)
     CALL mp_bcast(at, ionode_id, world_comm)
  ENDIF
  !
  CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
  alat = celldm(1)
  at = at / alat !  bring at in units of alat
  CALL volume(alat,at(1,1),at(1,2),at(1,3),omega)
  !
  !  read atomic types, positions and masses
  !
  DO nt = 1,ntyp
     IF (ionode) READ(1,*) i,atm,amass_from_file
     CALL mp_bcast(i,ionode_id,world_comm)
     CALL mp_bcast(atm,ionode_id,world_comm)
     CALL mp_bcast(amass_from_file,ionode_id,world_comm)
     IF (i.NE.nt) CALL errore ('readfc','wrong data read',nt)
     IF (amass(nt).EQ.0.d0) THEN
        amass(nt) = amass_from_file
     ELSE
        WRITE(stdout,*) 'for atomic type',nt,' mass from file not used'
     END IF
  END DO
  !
  ALLOCATE (tau(3,nat), ityp(nat), zeu(3,3,nat))
  !
  DO na=1,nat
     IF (ionode) READ(1,*) i,ityp(na),(tau(j,na),j=1,3)
     CALL mp_bcast(i,ionode_id,world_comm)
     IF (i.NE.na) CALL errore ('readfc','wrong data read',na)
  END DO
  CALL mp_bcast(ityp,ionode_id,world_comm)
  CALL mp_bcast(tau,ionode_id,world_comm)
  !
  !  read macroscopic variable
  !
  IF (ionode) READ (1,*) has_zstar
  CALL mp_bcast(has_zstar,ionode_id,world_comm)
  IF (has_zstar) THEN
     IF (ionode) READ(1,*) ((epsil(i,j),j=1,3),i=1,3)
     CALL mp_bcast(epsil,ionode_id,world_comm)
     IF (ionode) THEN
        DO na=1,nat
           READ(1,*)
           READ(1,*) ((zeu(i,j,na),j=1,3),i=1,3)
        END DO
     ENDIF
     CALL mp_bcast(zeu,ionode_id,world_comm)
  ELSE
     zeu  (:,:,:) = 0.d0
     epsil(:,:) = 0.d0
  END IF
  !
  IF (ionode) READ (1,*) nr1,nr2,nr3
  CALL mp_bcast(nr1,ionode_id,world_comm)
  CALL mp_bcast(nr2,ionode_id,world_comm)
  CALL mp_bcast(nr3,ionode_id,world_comm)
  !
  !  read real-space interatomic force constants
  !
  ALLOCATE ( frc(nr1,nr2,nr3,3,3,nat,nat) )
  frc(:,:,:,:,:,:,:) = 0.d0
  DO i=1,3
     DO j=1,3
        DO na=1,nat
           DO nb=1,nat
              IF (ionode) READ (1,*) ibid, jbid, nabid, nbbid
              CALL mp_bcast(ibid,ionode_id,world_comm)
              CALL mp_bcast(jbid,ionode_id,world_comm)
              CALL mp_bcast(nabid,ionode_id,world_comm)
              CALL mp_bcast(nbbid,ionode_id,world_comm)
              IF(i .NE.ibid  .OR. j .NE.jbid .OR.                   &
                 na.NE.nabid .OR. nb.NE.nbbid)                      &
                 CALL errore  ('readfc','error in reading',1)
              IF (ionode) READ (1,*) (((m1bid, m2bid, m3bid,        &
                          frc(m1,m2,m3,i,j,na,nb),                  &
                           m1=1,nr1),m2=1,nr2),m3=1,nr3)
               
              CALL mp_bcast(frc(:,:,:,i,j,na,nb),ionode_id,world_comm)
           END DO
        END DO
     END DO
  END DO
  !
  IF (ionode) CLOSE(unit=1)
  !
  RETURN
END SUBROUTINE readfc
!
!-----------------------------------------------------------------------
SUBROUTINE frc_blk(dyn,q,tau,nat,nr1,nr2,nr3,frc,at,bg,rws,nrws)
  !-----------------------------------------------------------------------
  ! calculates the dynamical matrix at q from the (short-range part of the)
  ! force constants
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  INTEGER nr1, nr2, nr3, nat, n1, n2, n3, &
          ipol, jpol, na, nb, m1, m2, m3, nint, i,j, nrws
  COMPLEX(DP) dyn(3,3,nat,nat)
  REAL(DP) frc(nr1,nr2,nr3,3,3,nat,nat), tau(3,nat), q(3), arg, &
               at(3,3), bg(3,3), r(3), weight, r_ws(3),  &
               total_weight, rws(0:3,nrws)
  REAL(DP), EXTERNAL :: wsweight
  !
  DO na=1, nat
     DO nb=1, nat
        total_weight=0.0d0
        DO n1=-3*nr1,3*nr1
           DO n2=-3*nr2,3*nr2
              DO n3=-3*nr3,3*nr3
                 !
                 ! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY SAFE RANGE!
                 !
                 DO i=1, 3
                    r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                    r_ws(i) = r(i) + tau(i,na)-tau(i,nb)
                 END DO
                 weight = wsweight(r_ws,rws,nrws)
                 IF (weight .GT. 0.0d0) THEN
                    !
                    ! FIND THE VECTOR CORRESPONDING TO R IN THE ORIGINAL CELL
                    !
                    m1 = MOD(n1+1,nr1)
                    IF(m1.LE.0) m1=m1+nr1
                    m2 = MOD(n2+1,nr2)
                    IF(m2.LE.0) m2=m2+nr2
                    m3 = MOD(n3+1,nr3)
                    IF(m3.LE.0) m3=m3+nr3
                    !
                    ! FOURIER TRANSFORM
                    !
                    arg = tpi*(q(1)*r(1) + q(2)*r(2) + q(3)*r(3))
                    DO ipol=1, 3
                       DO jpol=1, 3
                          dyn(ipol,jpol,na,nb) =                 &
                               dyn(ipol,jpol,na,nb) +            &
                               frc(m1,m2,m3,ipol,jpol,na,nb)     &
                               *CMPLX(COS(arg),-SIN(arg),kind=DP)*weight
                       END DO
                    END DO
                 END IF
                 total_weight=total_weight + weight
              END DO
           END DO
        END DO
        IF (ABS(total_weight-nr1*nr2*nr3).GT.1.0d-8) THEN
           WRITE(stdout,*) total_weight
           CALL errore ('frc_blk','wrong total_weight',1)
        END IF
     END DO
  END DO
  !
  RETURN
END SUBROUTINE frc_blk
!
!-----------------------------------------------------------------------
SUBROUTINE setupmat (q,dyn,nat,at,bg,tau,itau_blk,nsc,alat, &
     &         dyn_blk,nat_blk,at_blk,bg_blk,tau_blk,omega_blk, &
     &                 epsil,zeu,frc,nr1,nr2,nr3,has_zstar,rws,nrws)
  !-----------------------------------------------------------------------
  ! compute the dynamical matrix (the analytic part only)
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER:: nr1, nr2, nr3, nat, nat_blk, nsc, nrws, itau_blk(nat)
  REAL(DP) :: q(3), tau(3,nat), at(3,3), bg(3,3), alat,      &
                  epsil(3,3), zeu(3,3,nat_blk), rws(0:3,nrws),   &
                  frc(nr1,nr2,nr3,3,3,nat_blk,nat_blk)
  REAL(DP) :: tau_blk(3,nat_blk), at_blk(3,3), bg_blk(3,3), omega_blk
  COMPLEX(DP) dyn_blk(3,3,nat_blk,nat_blk)
  COMPLEX(DP) ::  dyn(3,3,nat,nat)
  LOGICAL has_zstar
  !
  ! local variables
  !
  REAL(DP) :: arg
  COMPLEX(DP) :: cfac(nat)
  INTEGER :: i,j,k, na,nb, na_blk, nb_blk, iq
  REAL(DP) qp(3), qbid(3,nsc) ! automatic array
  !
  !
  CALL q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !
  DO iq=1,nsc
     !
     DO k=1,3
        qp(k)= q(k) + qbid(k,iq)
     END DO
     !
     dyn_blk(:,:,:,:) = (0.d0,0.d0)
     CALL frc_blk (dyn_blk,qp,tau_blk,nat_blk,              &
          &              nr1,nr2,nr3,frc,at_blk,bg_blk,rws,nrws)
     IF (has_zstar) &
          CALL rgd_blk(nr1,nr2,nr3,nat_blk,dyn_blk,qp,tau_blk,   &
                       epsil,zeu,bg_blk,omega_blk,+1.d0)
     !
     DO na=1,nat
        na_blk = itau_blk(na)
        DO nb=1,nat
           nb_blk = itau_blk(nb)
           !
           arg=tpi* ( qp(1) * ( (tau(1,na)-tau_blk(1,na_blk)) -   &
                                (tau(1,nb)-tau_blk(1,nb_blk)) ) + &
                      qp(2) * ( (tau(2,na)-tau_blk(2,na_blk)) -   &
                                (tau(2,nb)-tau_blk(2,nb_blk)) ) + &
                      qp(3) * ( (tau(3,na)-tau_blk(3,na_blk)) -   &
                                (tau(3,nb)-tau_blk(3,nb_blk)) ) )
           !
           cfac(nb) = CMPLX(COS(arg),SIN(arg),kind=DP)/nsc
           !
        END DO ! nb
        !
        DO i=1,3
           DO j=1,3
              !
              DO nb=1,nat
                 nb_blk = itau_blk(nb)
                 dyn(i,j,na,nb) = dyn(i,j,na,nb) + cfac(nb) * &
                      dyn_blk(i,j,na_blk,nb_blk)
              END DO ! nb
              !
           END DO ! j
        END DO ! i
     END DO ! na
     !
  END DO ! iq
  !
  RETURN
END SUBROUTINE setupmat
!
!
!----------------------------------------------------------------------
SUBROUTINE set_asr (asr, nr1, nr2, nr3, frc, zeu, nat, ibrav, tau)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  CHARACTER (LEN=10), intent(in) :: asr
  INTEGER, intent(in) :: nr1, nr2, nr3, nat, ibrav
  REAL(DP), intent(in) :: tau(3,nat)
  REAL(DP), intent(inout) :: frc(nr1,nr2,nr3,3,3,nat,nat), zeu(3,3,nat)
  !
  INTEGER :: axis, n, i, j, na, nb, n1,n2,n3, m,p,k,l,q,r, i1,j1,na1
  REAL(DP) :: zeu_new(3,3,nat)
  REAL(DP), ALLOCATABLE :: frc_new(:,:,:,:,:,:,:)
  type vector
     real(DP),pointer :: vec(:,:,:,:,:,:,:)
  end type vector
  !
  type (vector) u(6*3*nat)
  ! These are the "vectors" associated with the sum rules on force-constants
  !
  integer :: u_less(6*3*nat),n_less,i_less
  ! indices of the vectors u that are not independent to the preceding ones,
  ! n_less = number of such vectors, i_less = temporary parameter
  !
  integer, allocatable :: ind_v(:,:,:)
  real(DP), allocatable :: v(:,:)
  ! These are the "vectors" associated with symmetry conditions, coded by
  ! indicating the positions (i.e. the seven indices) of the non-zero elements (there
  ! should be only 2 of them) and the value of that element. We do so in order
  ! to limit the amount of memory used.
  !
  real(DP), allocatable :: w(:,:,:,:,:,:,:), x(:,:,:,:,:,:,:)
  ! temporary vectors and parameters
  real(DP) :: scal,norm2, sum
  !
  real(DP) :: zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer :: zeu_less(6*3),nzeu_less,izeu_less
  ! indices of the vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(DP) :: zeu_w(3,3,nat), zeu_x(3,3,nat)
  ! temporary vectors

  ! Initialization. n is the number of sum rules to be considered (if asr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if axis='3')
  !
  if((asr.ne.'simple').and.(asr.ne.'crystal').and.(asr.ne.'one-dim') &
                      .and.(asr.ne.'zero-dim')) then
     call errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
  endif
  !
  if(asr.eq.'simple') then
     !
     ! Simple Acoustic Sum Rule on effective charges
     !
     do i=1,3
        do j=1,3
           sum=0.0d0
           do na=1,nat
              sum = sum + zeu(i,j,na)
           end do
           do na=1,nat
              zeu(i,j,na) = zeu(i,j,na) - sum/nat
           end do
        end do
     end do
     !
     ! Simple Acoustic Sum Rule on force constants in real space
     !
     do i=1,3
        do j=1,3
           do na=1,nat
              sum=0.0d0
               do nb=1,nat
                  do n1=1,nr1
                     do n2=1,nr2
                        do n3=1,nr3
                           sum=sum+frc(n1,n2,n3,i,j,na,nb)
                        end do
                     end do
                  end do
               end do
               frc(1,1,1,i,j,na,na) = frc(1,1,1,i,j,na,na) - sum
               !               write(6,*) ' na, i, j, sum = ',na,i,j,sum
            end do
         end do
      end do
      !
      return
      !
   end if

  if(asr.eq.'crystal') n=3
  if(asr.eq.'one-dim') then
     ! the direction of periodicity is the rotation axis
     ! It will work only if the crystal axis considered is one of
     ! the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
     ! z-direction)
     if (nr1*nr2*nr3.eq.1) axis=3
     if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis=1
     if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
     if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
     if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
        call errore('set_asr','too many directions of &
             & periodicity in 1D system',axis)
     endif
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        write(stdout,*) 'asr: rotational axis may be wrong'
     endif
     write(stdout,'("asr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(asr.eq.'zero-dim') n=6
  !
  ! Acoustic Sum Rule on effective charges
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the effective charges matrix on
  !
  zeu_u(:,:,:,:)=0.0d0
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu_new(i,j,na)=zeu(i,j,na)
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        ! These are the 3*3 vectors associated with the
        ! translational acoustic sum rules
        p=p+1
        zeu_u(p,i,j,:)=1.0d0
        !
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        ! These are the 3 vectors associated with the
        ! single rotational sum rule (1D system)
        p=p+1
        do na=1,nat
           zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
           zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
        enddo
        !
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           ! These are the 3*3 vectors associated with the
           ! three rotational sum rules (0D system - typ. molecule)
           p=p+1
           do na=1,nat
              zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
              zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
           enddo
           !
        enddo
     enddo
  endif
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  !
  nzeu_less=0
  do k=1,p
     zeu_w(:,:,:)=zeu_u(k,:,:,:)
     zeu_x(:,:,:)=zeu_u(k,:,:,:)
     do q=1,k-1
        r=1
        do izeu_less=1,nzeu_less
           if (zeu_less(izeu_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
           zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
        endif
     enddo
     call sp_zeu(zeu_w,zeu_w,nat,norm2)
     if (norm2.gt.1.0d-16) then
        zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
     else
        nzeu_less=nzeu_less+1
        zeu_less(nzeu_less)=k
     endif
  enddo
  !
  ! Projection of the effective charge "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules
  !
  zeu_w(:,:,:)=0.0d0
  do k=1,p
     r=1
     do izeu_less=1,nzeu_less
        if (zeu_less(izeu_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        zeu_x(:,:,:)=zeu_u(k,:,:,:)
        call sp_zeu(zeu_x,zeu_new,nat,scal)
        zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
     endif
  enddo
  !
  ! Final substraction of the former projection to the initial zeu, to get
  ! the new "projected" zeu
  !
  zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
  call sp_zeu(zeu_w,zeu_w,nat,norm2)
  write(stdout,'("Norm of the difference between old and new effective ", &
       & "charges: ",F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection of zeu")')
  !do k=1,p
  !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
  !  call sp_zeu(zeu_x,zeu_new,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu(i,j,na)=zeu_new(i,j,na)
        enddo
     enddo
  enddo
  !
  ! Acoustic Sum Rule on force constants
  !
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the force-constants matrix on
  !
  do k=1,18*nat
     allocate(u(k) % vec(nr1,nr2,nr3,3,3,nat,nat))
     u(k) % vec (:,:,:,:,:,:,:)=0.0d0
  enddo
  ALLOCATE (frc_new(nr1,nr2,nr3,3,3,nat,nat))
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       frc_new(n1,n2,n3,i,j,na,nb)=frc(n1,n2,n3,i,j,na,nb)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        do na=1,nat
           ! These are the 3*3*nat vectors associated with the
           ! translational acoustic sum rules
           p=p+1
           u(p) % vec (:,:,:,i,j,na,:)=1.0d0
           !
        enddo
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        do na=1,nat
           ! These are the 3*nat vectors associated with the
           ! single rotational sum rule (1D system)
           p=p+1
           do nb=1,nat
              u(p) % vec (:,:,:,i,MOD(axis,3)+1,na,nb)=-tau(MOD(axis+1,3)+1,nb)
              u(p) % vec (:,:,:,i,MOD(axis+1,3)+1,na,nb)=tau(MOD(axis,3)+1,nb)
           enddo
           !
        enddo
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           do na=1,nat
              ! These are the 3*3*nat vectors associated with the
              ! three rotational sum rules (0D system - typ. molecule)
              p=p+1
              do nb=1,nat
                 u(p) % vec (:,:,:,i,MOD(j,3)+1,na,nb)=-tau(MOD(j+1,3)+1,nb)
                 u(p) % vec (:,:,:,i,MOD(j+1,3)+1,na,nb)=tau(MOD(j,3)+1,nb)
              enddo
              !
           enddo
        enddo
     enddo
  endif
  !
  allocate (ind_v(9*nat*nat*nr1*nr2*nr3,2,7), v(9*nat*nat*nr1*nr2*nr3,2) )
  m=0
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       ! These are the vectors associated with the symmetry constraints
                       q=1
                       l=1
                       do while((l.le.m).and.(q.ne.0))
                          if ((ind_v(l,1,1).eq.n1).and.(ind_v(l,1,2).eq.n2).and. &
                               (ind_v(l,1,3).eq.n3).and.(ind_v(l,1,4).eq.i).and. &
                               (ind_v(l,1,5).eq.j).and.(ind_v(l,1,6).eq.na).and. &
                               (ind_v(l,1,7).eq.nb)) q=0
                          if ((ind_v(l,2,1).eq.n1).and.(ind_v(l,2,2).eq.n2).and. &
                               (ind_v(l,2,3).eq.n3).and.(ind_v(l,2,4).eq.i).and. &
                               (ind_v(l,2,5).eq.j).and.(ind_v(l,2,6).eq.na).and. &
                               (ind_v(l,2,7).eq.nb)) q=0
                          l=l+1
                       enddo
                       if ((n1.eq.MOD(nr1+1-n1,nr1)+1).and.(n2.eq.MOD(nr2+1-n2,nr2)+1) &
                            .and.(n3.eq.MOD(nr3+1-n3,nr3)+1).and.(i.eq.j).and.(na.eq.nb)) q=0
                       if (q.ne.0) then
                          m=m+1
                          ind_v(m,1,1)=n1
                          ind_v(m,1,2)=n2
                          ind_v(m,1,3)=n3
                          ind_v(m,1,4)=i
                          ind_v(m,1,5)=j
                          ind_v(m,1,6)=na
                          ind_v(m,1,7)=nb
                          v(m,1)=1.0d0/DSQRT(2.0d0)
                          ind_v(m,2,1)=MOD(nr1+1-n1,nr1)+1
                          ind_v(m,2,2)=MOD(nr2+1-n2,nr2)+1
                          ind_v(m,2,3)=MOD(nr3+1-n3,nr3)+1
                          ind_v(m,2,4)=j
                          ind_v(m,2,5)=i
                          ind_v(m,2,6)=nb
                          ind_v(m,2,7)=na
                          v(m,2)=-1.0d0/DSQRT(2.0d0)
                       endif
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  ! Note that the vectors corresponding to symmetry constraints are already
  ! orthonormalized by construction.
  !
  n_less=0
  allocate (w(nr1,nr2,nr3,3,3,nat,nat), x(nr1,nr2,nr3,3,3,nat,nat))
  do k=1,p
     w(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     do l=1,m
        !
        call sp2(x,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
        do r=1,2
           n1=ind_v(l,r,1)
           n2=ind_v(l,r,2)
           n3=ind_v(l,r,3)
           i=ind_v(l,r,4)
           j=ind_v(l,r,5)
           na=ind_v(l,r,6)
           nb=ind_v(l,r,7)
           w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)-scal*v(l,r)
        enddo
     enddo
     if (k.le.(9*nat)) then
        na1=MOD(k,nat)
        if (na1.eq.0) na1=nat
        j1=MOD((k-na1)/nat,3)+1
        i1=MOD((((k-na1)/nat)-j1+1)/3,3)+1
     else
        q=k-9*nat
        if (n.eq.4) then
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           i1=MOD((q-na1)/nat,3)+1
        else
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           j1=MOD((q-na1)/nat,3)+1
           i1=MOD((((q-na1)/nat)-j1+1)/3,3)+1
        endif
     endif
     do q=1,k-1
        r=1
        do i_less=1,n_less
           if (u_less(i_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp3(x,u(q) % vec (:,:,:,:,:,:,:), i1,na1,nr1,nr2,nr3,nat,scal)
           w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) - scal* u(q) % vec (:,:,:,:,:,:,:)
        endif
     enddo
     call sp1(w,w,nr1,nr2,nr3,nat,norm2)
     if (norm2.gt.1.0d-16) then
        u(k) % vec (:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) / DSQRT(norm2)
     else
        n_less=n_less+1
        u_less(n_less)=k
     endif
  enddo
  !
  ! Projection of the force-constants "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules and symmetry contraints
  !
  w(:,:,:,:,:,:,:)=0.0d0
  do l=1,m
     call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
     do r=1,2
        n1=ind_v(l,r,1)
        n2=ind_v(l,r,2)
        n3=ind_v(l,r,3)
        i=ind_v(l,r,4)
        j=ind_v(l,r,5)
        na=ind_v(l,r,6)
        nb=ind_v(l,r,7)
        w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)+scal*v(l,r)
     enddo
  enddo
  do k=1,p
     r=1
     do i_less=1,n_less
        if (u_less(i_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
        call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
        w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) + scal*u(k)%vec(:,:,:,:,:,:,:)
     endif
     deallocate(u(k) % vec)
  enddo
  !
  ! Final substraction of the former projection to the initial frc, to get
  ! the new "projected" frc
  !
  frc_new(:,:,:,:,:,:,:)=frc_new(:,:,:,:,:,:,:) - w(:,:,:,:,:,:,:)
  call sp1(w,w,nr1,nr2,nr3,nat,norm2)
  write(stdout,'("Norm of the difference between old and new force-constants:",&
       &     F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection IFC")')
  !do l=1,m
  !  call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("l= ",I8," frc_new|v(l)= ",F15.10)') l,scal
  !enddo
  !do k=1,p
  !  x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
  !  call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," frc_new|u(k)= ",F15.10)') k,scal
  !  deallocate(u(k) % vec)
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       frc(n1,n2,n3,i,j,na,nb)=frc_new(n1,n2,n3,i,j,na,nb)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  deallocate (x, w)
  deallocate (v, ind_v)
  deallocate (frc_new)
  !
  return
end subroutine set_asr
!
!----------------------------------------------------------------------
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  integer i,j,na,nat
  real(DP) zeu_u(3,3,nat)
  real(DP) zeu_v(3,3,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp_zeu
!
!
!----------------------------------------------------------------------
subroutine sp1(u,v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        do nb=1,nat
          do n1=1,nr1
            do n2=1,nr2
              do n3=1,nr3
                scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp1
!
!----------------------------------------------------------------------
subroutine sp2(u,v,ind_v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space). u is coded in the usual way
  ! but v is coded as explained when defining the vectors corresponding to the
  ! symmetry constraints
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  integer ind_v(2,7)
  real(DP) v(2)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,2
    scal=scal+u(ind_v(i,1),ind_v(i,2),ind_v(i,3),ind_v(i,4),ind_v(i,5),ind_v(i,6), &
         ind_v(i,7))*v(i)
  enddo
  !
  return
  !
end subroutine sp2
!
!----------------------------------------------------------------------
subroutine sp3(u,v,i,na,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! like sp1, but in the particular case when u is one of the u(k)%vec
  ! defined in set_asr (before orthonormalization). In this case most of the
  ! terms are zero (the ones that are not are characterized by i and na), so
  ! that a lot of computer time can be saved (during Gram-Schmidt).
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do j=1,3
    do nb=1,nat
      do n1=1,nr1
        do n2=1,nr2
          do n3=1,nr3
            scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp3
!
!-----------------------------------------------------------------------
SUBROUTINE q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !-----------------------------------------------------------------------
  ! generate list of q (qbid) that are G-vectors of the supercell
  ! but not of the bulk
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER :: nsc
  REAL(DP) qbid(3,nsc), at_blk(3,3), bg_blk(3,3), at(3,3), bg(3,3)
  !
  INTEGER, PARAMETER:: nr1=4, nr2=4, nr3=4, &
                       nrm=(2*nr1+1)*(2*nr2+1)*(2*nr3+1)
  REAL(DP), PARAMETER:: eps=1.0d-7
  INTEGER :: i, j, k,i1, i2, i3, idum(nrm), iq
  REAL(DP) :: qnorm(nrm), qbd(3,nrm) ,qwork(3), delta
  LOGICAL lbho
  !
  i = 0
  DO i1=-nr1,nr1
     DO i2=-nr2,nr2
        DO i3=-nr3,nr3
           i = i + 1
           DO j=1,3
              qwork(j) = i1*bg(j,1) + i2*bg(j,2) + i3*bg(j,3)
           END DO ! j
           !
           qnorm(i)  = qwork(1)**2 + qwork(2)**2 + qwork(3)**2
           !
           DO j=1,3
              !
              qbd(j,i) = at_blk(1,j)*qwork(1) + &
                         at_blk(2,j)*qwork(2) + &
                         at_blk(3,j)*qwork(3)
           END DO ! j
           !
           idum(i) = 1
           !
        END DO ! i3
     END DO ! i2
  END DO ! i1
  !
  DO i=1,nrm-1
     IF (idum(i).EQ.1) THEN
        DO j=i+1,nrm
           IF (idum(j).EQ.1) THEN
              lbho=.TRUE.
              DO k=1,3
                 delta = qbd(k,i)-qbd(k,j)
                 lbho = lbho.AND. (ABS(NINT(delta)-delta).LT.eps)
              END DO ! k
              IF (lbho) THEN
                 IF(qnorm(i).GT.qnorm(j)) THEN
                    qbd(1,i) = qbd(1,j)
                    qbd(2,i) = qbd(2,j)
                    qbd(3,i) = qbd(3,j)
                    qnorm(i) = qnorm(j)
                 END IF
                 idum(j) = 0
              END IF
           END IF
        END DO ! j
     END IF
  END DO ! i
  !
  iq = 0
  DO i=1,nrm
     IF (idum(i).EQ.1) THEN
        iq=iq+1
        qbid(1,iq)= bg_blk(1,1)*qbd(1,i) +  &
                    bg_blk(1,2)*qbd(2,i) +  &
                    bg_blk(1,3)*qbd(3,i)
        qbid(2,iq)= bg_blk(2,1)*qbd(1,i) +  &
                    bg_blk(2,2)*qbd(2,i) +  &
                    bg_blk(2,3)*qbd(3,i)
        qbid(3,iq)= bg_blk(3,1)*qbd(1,i) +  &
                    bg_blk(3,2)*qbd(2,i) +  &
                    bg_blk(3,3)*qbd(3,i)
     END IF
  END DO ! i
  !
  IF (iq.NE.nsc) CALL errore('q_gen',' probably nr1,nr2,nr3 too small ', iq)
  RETURN
END SUBROUTINE q_gen
!
!-----------------------------------------------------------------------
SUBROUTINE check_at(at,bg_blk,alat,omega)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP) :: at(3,3), bg_blk(3,3), alat, omega
  REAL(DP) :: work(3,3)
  INTEGER :: i,j
  REAL(DP), PARAMETER :: small=1.d-6
  !
  work(:,:) = at(:,:)
  CALL cryst_to_cart(3,work,bg_blk,-1)
  !
  DO j=1,3
     DO i =1,3
        IF ( ABS(work(i,j)-NINT(work(i,j))) > small) THEN
           WRITE (stdout,'(3f9.4)') work(:,:)
           CALL errore ('check_at','at not multiple of at_blk',1)
        END IF
     END DO
  END DO
  !
  omega =alat**3 * ABS(at(1,1)*(at(2,2)*at(3,3)-at(3,2)*at(2,3))- &
                       at(1,2)*(at(2,1)*at(3,3)-at(2,3)*at(3,1))+ &
                       at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))
  !
  RETURN
END SUBROUTINE check_at
!
!-----------------------------------------------------------------------
SUBROUTINE set_tau (nat, nat_blk, at, at_blk, tau, tau_blk, &
     ityp, ityp_blk, itau_blk)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER nat, nat_blk,ityp(nat),ityp_blk(nat_blk), itau_blk(nat)
  REAL(DP) at(3,3),at_blk(3,3),tau(3,nat),tau_blk(3,nat_blk)
  !
  REAL(DP) bg(3,3), r(3) ! work vectors
  INTEGER i,i1,i2,i3,na,na_blk
  REAL(DP) small
  INTEGER NN1,NN2,NN3
  PARAMETER (NN1=8, NN2=8, NN3=8, small=1.d-8)
  !
  CALL recips (at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
  !
  na = 0
  !
  DO i1 = -NN1,NN1
     DO i2 = -NN2,NN2
        DO i3 = -NN3,NN3
           r(1) = i1*at_blk(1,1) + i2*at_blk(1,2) + i3*at_blk(1,3)
           r(2) = i1*at_blk(2,1) + i2*at_blk(2,2) + i3*at_blk(2,3)
           r(3) = i1*at_blk(3,1) + i2*at_blk(3,2) + i3*at_blk(3,3)
           CALL cryst_to_cart(1,r,bg,-1)
           !
           IF ( r(1).GT.-small .AND. r(1).LT.1.d0-small .AND.          &
                r(2).GT.-small .AND. r(2).LT.1.d0-small .AND.          &
                r(3).GT.-small .AND. r(3).LT.1.d0-small ) THEN
              CALL cryst_to_cart(1,r,at,+1)
              !
              DO na_blk=1, nat_blk
                 na = na + 1
                 IF (na.GT.nat) CALL errore('set_tau','too many atoms',na)
                 tau(1,na)    = tau_blk(1,na_blk) + r(1)
                 tau(2,na)    = tau_blk(2,na_blk) + r(2)
                 tau(3,na)    = tau_blk(3,na_blk) + r(3)
                 ityp(na)     = ityp_blk(na_blk)
                 itau_blk(na) = na_blk
              END DO
              !
           END IF
           !
        END DO
     END DO
  END DO
  !
  IF (na.NE.nat) CALL errore('set_tau','too few atoms: increase NNs',na)
  !
  RETURN
END SUBROUTINE set_tau
!
!-----------------------------------------------------------------------
SUBROUTINE read_tau &
     (nat, nat_blk, ntyp, bg_blk, tau, tau_blk, ityp, itau_blk)
  !---------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode_id, ionode
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  !
  IMPLICIT NONE
  !
  INTEGER nat, nat_blk, ntyp, ityp(nat),itau_blk(nat)
  REAL(DP) bg_blk(3,3),tau(3,nat),tau_blk(3,nat_blk)
  !
  REAL(DP) r(3) ! work vectors
  INTEGER i,na,na_blk
  !
  REAL(DP) small
  PARAMETER ( small = 1.d-6 )
  !
  DO na=1,nat
     IF (ionode) READ(5,*) (tau(i,na),i=1,3), ityp(na)
     CALL mp_bcast(tau(:,na),ionode_id,world_comm)
     CALL mp_bcast(ityp(na),ionode_id,world_comm)
     IF (ityp(na).LE.0 .OR. ityp(na) .GT. ntyp) &
          CALL errore('read_tau',' wrong atomic type', na)
     DO na_blk=1,nat_blk
        r(1) = tau(1,na) - tau_blk(1,na_blk)
        r(2) = tau(2,na) - tau_blk(2,na_blk)
        r(3) = tau(3,na) - tau_blk(3,na_blk)
        CALL cryst_to_cart(1,r,bg_blk,-1)
        IF (ABS( r(1)-NINT(r(1)) ) .LT. small .AND.                 &
            ABS( r(2)-NINT(r(2)) ) .LT. small .AND.                 &
            ABS( r(3)-NINT(r(3)) ) .LT. small ) THEN
           itau_blk(na) = na_blk
           go to 999
        END IF
     END DO
     CALL errore ('read_tau',' wrong atomic position ', na)
999  CONTINUE
  END DO
  !
  RETURN
END SUBROUTINE read_tau
!
!-----------------------------------------------------------------------
SUBROUTINE write_tau(fltau,nat,tau,ityp)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,   ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER nat, ityp(nat)
  REAL(DP) tau(3,nat)
  CHARACTER(LEN=*) fltau
  !
  INTEGER i,na
  !
  IF (.NOT.ionode) RETURN
  OPEN (unit=4,file=fltau, status='new')
  DO na=1,nat
     WRITE(4,'(3(f12.6),i3)') (tau(i,na),i=1,3), ityp(na)
  END DO
  CLOSE (4)
  !
  RETURN
END SUBROUTINE write_tau
!
!-----------------------------------------------------------------------
SUBROUTINE gen_qpoints (ibrav, at_, bg_, nat, tau, ityp, nk1, nk2, nk3, &
     symm_type, ntetra, nqx, nq, q, tetra)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : set_sym_bl, find_sym, s, ftau, irt, nsym, &
                         nrot, t_rev, time_reversal,  sname
  !
  IMPLICIT NONE
  ! input
  INTEGER :: ibrav, nat, nk1, nk2, nk3, ntetra, ityp(*)
  REAL(DP) :: at_(3,3), bg_(3,3), tau(3,nat)
  character(LEN=9)    :: symm_type
  ! output
  INTEGER :: nqx, nq, tetra(4,ntetra)
  REAL(DP) :: q(3,nqx)
  ! local
  LOGICAL :: nosym_evc=.false., nofrac=.false.
  REAL(DP) :: xqq(3), wk(nqx), mdum(3,nat)
  !
  time_reversal = .true.
  t_rev(:) = 0
  xqq (:) =0.d0
  at = at_
  bg = bg_
  CALL set_sym_bl()
  !
  CALL kpoint_grid ( nrot, time_reversal, s, t_rev, bg, nqx, &
                           0,0,0, nk1,nk2,nk3, nq, q, wk)
  !
  CALL find_sym ( nat, tau, ityp, 6, 6, 6, .not.time_reversal, mdum )
  !
  CALL irreducible_BZ (nrot, s, nsym, time_reversal, at, bg, nqx, nq, q, wk, &
                       t_rev)
  !
  IF (ntetra /= 6 * nk1 * nk2 * nk3) &
       CALL errore ('gen_qpoints','inconsistent ntetra',1)
  !
  CALL tetrahedra (nsym, s, time_reversal, at, bg, nqx, 0, 0, 0, &
       nk1, nk2, nk3, nq, q, wk, ntetra, tetra)
  !
  RETURN
END SUBROUTINE gen_qpoints
!
!
!
!-----------------------------------------------------------------------
subroutine readfg ( ifn, nr1, nr2, nr3, nat, frcg )
  !-----------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : ionode, ionode_id, stdout
  USE mp,          ONLY : mp_bcast
  USE mp_world,    ONLY : world_comm
  implicit none
  ! I/O variable
  integer, intent(in) ::  nr1,nr2,nr3, nat
  real(DP), intent(out) :: frcg(nr1,nr2,nr3,3,3,nat,nat)
  ! local variables
  integer i, j, na, nb, m1,m2,m3, ifn
  integer ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
  !
  !
  IF (ionode) READ (ifn,*) m1, m2, m3
  CALL mp_bcast(m1, ionode_id, world_comm)
  CALL mp_bcast(m2, ionode_id, world_comm)
  CALL mp_bcast(m3, ionode_id, world_comm)
  if ( m1 /= nr1 .or. m2 /= nr2 .or. m3 /= nr3) &
       call errore('readfg','inconsistent nr1, nr2, nr3 read',1)
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              IF (ionode) read (ifn,*) ibid, jbid, nabid, nbbid
              CALL mp_bcast(ibid, ionode_id, world_comm)
              CALL mp_bcast(jbid, ionode_id, world_comm)
              CALL mp_bcast(nabid, ionode_id, world_comm)
              CALL mp_bcast(nbbid, ionode_id, world_comm)
              
              if(i.ne.ibid.or.j.ne.jbid.or.na.ne.nabid.or.nb.ne.nbbid)  then
                  write(stdout,*) i,j,na,nb,'  <>  ', ibid, jbid, nabid, nbbid
                  call errore  ('readfG','error in reading',1)
              else
                  IF (ionode) read (ifn,*) (((m1bid, m2bid, m3bid,     &
                                 frcg(m1,m2,m3,i,j,na,nb), &
                                 m1=1,nr1),m2=1,nr2),m3=1,nr3)
              endif
              CALL mp_bcast(frcg(:,:,:,i,j,na,nb), ionode_id, world_comm)
           end do
        end do
     end do
  end do
  !
  IF (ionode) close(ifn)
  !
  return
end subroutine readfg
!
!

end program elph_fc

