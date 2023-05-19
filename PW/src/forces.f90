!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE forces()
  !----------------------------------------------------------------------------
  !
  ! ... This routine is a driver routine which computes the forces
  ! ... acting on the atoms. The complete expression of the forces
  ! ... contains four parts which are computed by different routines:
  !
  ! ...  a)  force_lc,    local contribution to the forces
  ! ...  b)  force_cc,    contribution due to NLCC
  ! ...  c)  force_ew,    contribution due to the electrostatic ewald term
  ! ...  d)  force_us,    contribution due to the non-local potential
  ! ...  e)  force_corr,  correction term for incomplete self-consistency
  ! ...  f)  force_hub,   contribution due to the Hubbard term
  ! ...  g)  force_london, semi-empirical correction for dispersion forces
  !
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE cell_base,     ONLY : at, bg, alat, omega  
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, zv, amass, extfor
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : ngm, gstart, ngl, nl, igtongl, g, gg, gcutm
  USE lsda_mod,      ONLY : nspin
  USE symme,         ONLY : symvector
  USE vlocal,        ONLY : strf, vloc
  USE force_mod,     ONLY : force, lforce
  USE scf,           ONLY : rho
  USE ions_base,     ONLY : if_pos
  USE ldaU,          ONLY : lda_plus_u, U_projection
  USE extfield,      ONLY : tefield, forcefield
  USE control_flags, ONLY : gamma_only, remove_rigid_rot, textfor, &
                            iverbosity, llondon, lxdm, ts_vdw
#ifdef __ENVIRON
  USE environ_base,  ONLY : do_environ, env_static_permittivity, rhopol
  USE environ_base,  ONLY : env_extcharge_n, rhoexternal
  USE fft_interfaces,  ONLY : fwfft
  USE environ_main,  ONLY : calc_fenviron
#endif
  USE bp,            ONLY : lelfield, gdir, l3dstring, efield_cart, &
                            efield_cry,efield
  USE uspp,          ONLY : okvan
  USE martyna_tuckerman, ONLY: do_comp_mt, wg_corr_force
  USE london_module, ONLY : force_london
  USE xdm_module,    ONLY : force_xdm
  USE tsvdw_module,  ONLY : FtsvdW
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: forcenl(:,:), &
                           forcelc(:,:), &
                           forcecc(:,:), &
                           forceion(:,:), &
                           force_disp(:,:),&
                           force_disp_xdm(:,:),&
                           force_mt(:,:), &
                           forcescc(:,:), &
                           forces_bp_efield(:,:), &
                           forceh(:,:)
    ! nonlocal, local, core-correction, ewald, scf correction terms, and hubbard
#ifdef __ENVIRON
  REAL(DP), ALLOCATABLE :: force_environ(:,:), force_tmp(:,:)
#endif
!
! aux is used to store a possible additional density
! now defined in real space
!
  COMPLEX(DP), ALLOCATABLE :: auxg(:), auxr(:)
!
  REAL(DP) :: sumfor, sumscf, sum_mm
  REAL(DP),PARAMETER :: eps = 1.e-12_dp
  INTEGER  :: ipol, na
    ! counter on polarization
    ! counter on atoms
  !
  !
  CALL start_clock( 'forces' )
  !
  ALLOCATE( forcenl( 3, nat ), forcelc( 3, nat ), forcecc( 3, nat ), &
            forceh( 3, nat ), forceion( 3, nat ), forcescc( 3, nat ) )
  !    
  forcescc(:,:) = 0.D0
  forceh(:,:)   = 0.D0
  force (:,:)   = 0.D0
  !
  WRITE( stdout, '(/,5x,"Forces acting on atoms (Ry/au):", / )')
  !
  ! ... The nonlocal contribution is computed here
  !
  CALL force_us( forcenl )
  !
  ! ... The local contribution
  !
  CALL force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
                 g, rho%of_r, nl, nspin, gstart, gamma_only, vloc, &
                 forcelc )
  !
  ! ... The NLCC contribution
  !
  CALL force_cc( forcecc )
  !
  ! ... The Hubbard contribution
  !     (included by force_us if using beta as local projectors)
  !
  IF ( lda_plus_u .AND. U_projection.NE.'pseudo' ) CALL force_hub( forceh )
  !
  ! ... The ionic contribution is computed here
  !
  CALL force_ew( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
                 gg, ngm, gstart, gamma_only, gcutm, strf, forceion )
  !
  ! ... the semi-empirical dispersion correction
  !
  IF ( llondon ) THEN
    !
    ALLOCATE ( force_disp ( 3 , nat ) )
    force_disp ( : , : ) = 0.0_DP
    force_disp = force_london( alat , nat , ityp , at , bg , tau )
    !
  END IF
  IF (lxdm) THEN
     ALLOCATE (force_disp_xdm(3,nat))
     force_disp_xdm = 0._dp
     force_disp_xdm = force_xdm(nat)
  end if
     
  !
  ! ... The SCF contribution
  !
  CALL force_corr( forcescc )
  !
  IF (do_comp_mt) THEN
    !
    ALLOCATE ( force_mt ( 3 , nat ) )
    CALL wg_corr_force( .true.,omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, &
                        nspin, rho%of_g, force_mt )
  END IF
  !
#ifdef __ENVIRON
  IF (do_environ) THEN
    !
    ALLOCATE(force_tmp(3,nat))
    ALLOCATE(force_environ(3,nat))
    !
    force_environ=0.0_dp
    !
    if(do_comp_mt) then
      force_tmp=0.0_dp
      ALLOCATE(auxr(dfftp%nnr))
      ALLOCATE(auxg(ngm))
      auxg = CMPLX(0.0_dp,0.0_dp)
      auxr = CMPLX(0.0_dp,0.0_dp)
      if(env_static_permittivity .GT. 1.D0) auxr = CMPLX(rhopol(:),0.0, kind=DP)
      if(env_extcharge_n .GT. 0) auxr = auxr + CMPLX(rhoexternal(:),0.0, kind=DP) 
      CALL fwfft ('Dense', auxr, dfftp)
      auxg(:)=auxr(nl(:))
      CALL wg_corr_force(.false.,omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, &
                        1, auxg, force_tmp)
      force_environ = force_environ + force_tmp
      DEALLOCATE(auxr,auxg)
      WRITE( stdout, '(5x,"Tmp environment correction to forces")')
      DO na = 1, nat
         WRITE( stdout, 9035) na, ityp(na), ( force_tmp(ipol,na), ipol = 1, 3 )
      END DO
    endif 
    ! ... Computes here the solvent contribution
    !
    IF ( env_static_permittivity .GT. 1.D0 ) THEN
      force_tmp=0.0_dp 
      CALL force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
                     g, rhopol, nl, 1, gstart, gamma_only, vloc, &
                     force_tmp )
      force_environ = force_environ + force_tmp
    ENDIF
    ! 
    ! ... Computes here the external charges contribution
    !
    IF ( env_extcharge_n .GT. 0 ) THEN 
      force_tmp = 0.0_DP
      CALL force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
                   g, rhoexternal, nl, 1, gstart, gamma_only, vloc, &
                   force_tmp )
      force_environ = force_environ + force_tmp
    ENDIF
    !
    ! ... Add the other environment contributions
    !
    CALL calc_fenviron( dfftp%nnr, nspin, nat, rho%of_r, force_environ )
    !
    DEALLOCATE(force_tmp)
    !
    WRITE( stdout, '(5x,"The external environment correction to forces")')
    DO na = 1, nat
       WRITE( stdout, 9035) na, ityp(na), ( force_environ(ipol,na), ipol = 1, 3 )
    END DO
    !
    force = force_environ
    !
    DEALLOCATE(force_environ)
  END IF
  !
#endif
  !
  ! Berry's phase electric field terms
  !
  if(lelfield) then
     ALLOCATE ( forces_bp_efield (3,nat) )
     forces_bp_efield(:,:)=0.d0
     if(.not.l3dstring) then
        if(okvan) call  forces_us_efield(forces_bp_efield,gdir,efield)
        call forces_ion_efield(forces_bp_efield,gdir,efield)
     else
        if(okvan)then
           do ipol=1,3
              call  forces_us_efield(forces_bp_efield,ipol,efield_cry(ipol))
           enddo
        endif
        do ipol=1,3
           call  forces_ion_efield(forces_bp_efield,ipol,efield_cart(ipol))
        enddo
     endif
  endif
  !
  ! ... here we sum all the contributions and compute the total force acting
  ! ... on the crystal
  !
  DO ipol = 1, 3
     !
     sumfor = 0.D0
     !
     DO na = 1, nat
        !
        force(ipol,na) = force(ipol,na)    + &
                         forcenl(ipol,na)  + &
                         forceion(ipol,na) + &
                         forcelc(ipol,na)  + &
                         forcecc(ipol,na)  + &
                         forceh(ipol,na)   + &
                         forcescc(ipol,na)
        !
        IF ( llondon ) force(ipol,na) = force(ipol,na) + force_disp(ipol,na)
        IF ( lxdm )    force(ipol,na) = force(ipol,na) + force_disp_xdm(ipol,na)
        ! factor 2 converts from Ha to Ry a.u.
        IF ( ts_vdw )  force(ipol,na) = force(ipol,na) + 2.0_dp*FtsvdW(ipol,na)
        IF ( tefield ) force(ipol,na) = force(ipol,na) + forcefield(ipol,na)
        IF (lelfield)  force(ipol,na) = force(ipol,na) + forces_bp_efield(ipol,na)
        IF (do_comp_mt)force(ipol,na) = force(ipol,na) + force_mt(ipol,na) 

        sumfor = sumfor + force(ipol,na)
        !
     END DO
     !
     ! ... impose total force = 0
     !
     DO na = 1, nat
        !
        force(ipol,na) = force(ipol,na) - sumfor / DBLE( nat ) 
        !
     END DO
     !
#ifdef __MS2
     !
     ! ... impose total force of the quantum subsystem /= 0
     !
     DO na = 1, nat
        !
        force(ipol,na) = force(ipol,na) + sumfor / DBLE( nat )
        !
     END DO
     !
#endif
     !
  END DO
  !
  ! ... resymmetrize (should not be needed, but ...)
  !
  CALL symvector ( nat, force )
  !
  IF ( remove_rigid_rot ) &
     CALL remove_tot_torque( nat, tau, amass(ityp(:)), force  )
  !
  IF( textfor ) force(:,:) = force(:,:) + extfor(:,:)
  !
  ! ... call void routine for user define/ plugin patches on forces
  !
  CALL plugin_forces()
  !
  ! ... write on output the forces
  !
  DO na = 1, nat
     !
     WRITE( stdout, 9035) na, ityp(na), force(:,na)
     !
  END DO
  !
  ! ... forces on fixed coordinates are set to zero ( C.S. 15/10/2003 )
  !
  force(:,:)    = force(:,:)    * DBLE( if_pos )
  forcescc(:,:) = forcescc(:,:) * DBLE( if_pos )
  !
  IF ( iverbosity > 0 ) THEN
     IF ( do_comp_mt ) THEN
        WRITE( stdout, '(5x,"The Martyna-Tuckerman correction term to forces")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), ( force_mt(ipol,na), ipol = 1, 3 )
        END DO
     END IF
     !
     WRITE( stdout, '(5x,"The non-local contrib.  to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcenl(ipol,na), ipol = 1, 3 )
     END DO
     WRITE( stdout, '(5x,"The ionic contribution  to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forceion(ipol,na), ipol = 1, 3 )
     END DO
     WRITE( stdout, '(5x,"The local contribution  to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcelc(ipol,na), ipol = 1, 3 )
     END DO
     WRITE( stdout, '(5x,"The core correction contribution to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcecc(ipol,na), ipol = 1, 3 )
     END DO
     WRITE( stdout, '(5x,"The Hubbard contrib.    to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forceh(ipol,na), ipol = 1, 3 )
     END DO
     WRITE( stdout, '(5x,"The SCF correction term to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcescc(ipol,na), ipol = 1, 3 )
     END DO
     !
     IF ( llondon) THEN
        WRITE( stdout, '(/,5x,"Dispersion contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (force_disp(ipol,na), ipol = 1, 3)
        END DO
     END IF
     !
     IF (lxdm) THEN
        WRITE( stdout, '(/,5x,"XDM contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (force_disp_xdm(ipol,na), ipol = 1, 3)
        END DO
     END IF
     !
     IF ( ts_vdw) THEN
        WRITE( stdout, '(/,5x,"TS-VDW contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (2.0d0*FtsvdW(ipol,na), ipol=1,3)
        END DO
     END IF
     !
  END IF
  !
  sumfor = 0.D0
  sumscf = 0.D0
  !
  DO na = 1, nat
     !
     sumfor = sumfor + force(1,na)**2 + force(2,na)**2 + force(3,na)**2
     sumscf = sumscf + forcescc(1,na)**2 + forcescc(2,na)**2+ forcescc(3,na)**2
     !
  END DO
  !
  sumfor = SQRT( sumfor )
  sumscf = SQRT( sumscf )
  !
  WRITE( stdout, '(/5x,"Total force = ",F12.6,5X, &
              &  "Total SCF correction = ",F12.6)') sumfor, sumscf
  !
  IF ( llondon .AND. iverbosity > 0 ) THEN
     !
     sum_mm = 0.D0
     DO na = 1, nat
        sum_mm = sum_mm + &
                 force_disp(1,na)**2 + force_disp(2,na)**2 + force_disp(3,na)**2
     END DO
     sum_mm = SQRT( sum_mm )
     WRITE ( stdout, '(/,5x, "Total Dispersion Force = ",F12.6)') sum_mm
     !
  END IF
  !
  IF ( lxdm .AND. iverbosity > 0 ) THEN
     !
     sum_mm = 0.D0
     DO na = 1, nat
        sum_mm = sum_mm + &
                 force_disp_xdm(1,na)**2 + force_disp_xdm(2,na)**2 + force_disp_xdm(3,na)**2
     END DO
     sum_mm = SQRT( sum_mm )
     WRITE ( stdout, '(/,5x, "Total XDM Force = ",F12.6)') sum_mm
     !
  END IF
  !
  DEALLOCATE( forcenl, forcelc, forcecc, forceh, forceion, forcescc )
  IF ( llondon )  DEALLOCATE ( force_disp )
  IF ( lxdm ) DEALLOCATE( force_disp_xdm ) 
  IF ( lelfield ) DEALLOCATE ( forces_bp_efield )
  !
  lforce = .TRUE.
  !
  CALL stop_clock( 'forces' )
  !
  IF ( ( sumfor < 10.D0*sumscf ) .AND. ( sumfor > eps ) ) &
  WRITE( stdout,'(5x,"SCF correction compared to forces is large: ", &
                   &  "reduce conv_thr to get better values")')
  !
  IF(ALLOCATED(force_mt))   DEALLOCATE( force_mt )

  RETURN
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE forces