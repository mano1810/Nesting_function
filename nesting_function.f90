!! 2020.12 MANO Poobodin
!===============================================
! Fermi surface calculating program
! & Nesting function calculating program
! Read eigen values from Quantum Espresso / VASP output
! READED files [ Quantum Espresso ]
! 1. prefix.band       < from bands.x   >
! 2. prefix.band.proj  < from projwfc.x >
! READED files [ VASP ]
! 1. PROCAR
!===============================================
!
! &input_plot
! file_read     : prefix of files to be read
! file_k        : name of Fermi surface output file
! mode          : 'qe' , 'vasp'
! ldecom        : 0 for non-decommposed,
!                 1 for decomposed Fermi surface
! lsuscep       : 0 for no lindhard or nesting,
!                 1 for including lindhard or nesting
! ef            : Fermi energy
! degauss       : degauss for constaining calculation to Fermi level
! thresh        : exclude small value from output-weight
! n_con         : number of orbital to be considered in decomposed bands
! /
! LIST OF DECOMPOSED BANDS
! &suscep_q
! mode          : 'lindhard' or 'nesting'
! k_resoved     : 0 : NO , 1 : YES
! ibnd_start    : start band in band loop
! ibnd_final    : end band in band loop (This makes calculation faster)
! ndegauss      : number of degauss used in calculation
! degauss_init  : initial value for degauss
! delta_gauss   : interval for each degauss
! window        : energy window at Fermi level
! divide_k      : use subset of total k-point
! divide_q      : use subset of total q-point. Useful in convergence test
!                 and in speeding up to calculation
! file_q        : output file for lindhard/nesting function
! file_conv     : output file for convergence test
! omega         : coefficient for imaginary part
! /
! ------------------------------------------------ MAIN PROGRAM
!
program read_band
  !
  implicit none
  !
  integer :: nbnd, nks, ios, ink, ibnd, i,                     &
             idum, nat, ntyp, fdum, natomwfc,                  &
             nks_rap, nbnd_rap, iwfc, start_band ,             &
             end_band, ldecom, lsuscep, n_con, jnk, inat
  integer, allocatable :: orbital_limit(:)
  real(kind=8) , allocatable :: e(:,:) , k(:,:), sumproj(:,:), &
             selected(:) , tmp(:), proj_vasp(:)
  real(kind=8) :: proj, degauss, ef, weight, thresh
  real(kind=8), external :: wgauss, w0gauss
  NAMELIST /plot/ nks, nbnd
  character(len=80) :: file_read, filename2, file_k
  character(len=30) temp_1, temp_2, temp_3, dft
  NAMELIST /input_plot/ file_read, file_k, ldecom, lsuscep,    &
                        ef, degauss, thresh, n_con, dft
  !
  !--------- defalut
  file_k       = 'FERMI.dat'
  ldecom       = 0
  lsuscep      = 1
  degauss      = 0.02
  thresh       = 0.00
  dft          = 'qe'
  !------------
  !
  ! open standard input
  open(unit=111, file = 'input' , action = 'read' )
  read(111, input_plot)
  !
  allocate( orbital_limit(n_con) )
  read(111, * ) (orbital_limit(i) , i = 1, n_con )
  close(111)
  !
  ! open band structure file

  !------------ QE

  if (dft .eq. 'qe' ) then

    write(*,*) 'DFT mode              : Quantum Espresso'

    open( unit = 112, file = file_read , action = 'read' , iostat=ios)
    if (ios/=0) STOP 'Error opening file : .band'
    read(112,  plot, iostat=ios)
    !
    allocate( k(3,nks) )
    allocate( e(nbnd,nks) )
    !
    do ink  = 1, nks
      read( 112, * ) ( k(i,ink), i=1,3 )
      read( 112, * ) ( e(i,ink), i=1,nbnd )
    enddo
    close(112)
    !
    write(*,*) 'Read from             :     ', trim( file_read )
    !
    ! open projected file
    filename2=trim(file_read)//".proj"
    open(unit=113 , file = filename2 , action = 'read' , iostat=ios)
    if (ios/=0) STOP 'Error opening file : .proj'

    !-------- start of reading unnecessary information
    read(113, *, iostat=ios) ! empty line
    read(113, '(8i8)', iostat=ios) idum, idum, idum, idum, idum, &
             idum, nat, ntyp ! FFT grid (ignored), nat, ntyp
    read(113, '(i6,6f12.8)', iostat=ios) idum, fdum, fdum, fdum, fdum, &
           fdum, fdum ! ibrav, celldm(1-6)

    if (idum == 0) then ! read and discard the three cell vectors
      read(113, *, iostat=ios)
      read(113, *, iostat=ios)
      read(113, *, iostat=ios)
    endif

    do i=1,1+nat+ntyp
      read(113, *, iostat=ios) ! discard atomic positions
    enddo
    read(113, *, iostat=ios)
    !-------- end of reading unnecessary information

    read(113, '(3i8)', iostat=ios) natomwfc, nks_rap, nbnd_rap
    read(113, *,  iostat=ios) ! discard another line

    if (nks_rap/=nks.or.nbnd_rap/=nbnd.or.ios/=0) then
      write(*,'("file with projections not compatible with bands")')
      STOP
    else
      allocate( sumproj(nbnd,nks) )
      write(*,*) 'Number of wfc         : ', natomwfc
      write(*,*) 'Number of nk          : ', nks_rap
      write(*,*) 'Number of band        : ', nbnd_rap
    endif

    ! -------- read weight for each orbital from standard input
    allocate( selected(natomwfc ) )
    selected( : )  = 0.0
    do i = 1, n_con
      idum = orbital_limit(i)
      selected( idum ) = 1.0
    enddo
    write(*,*) 'Number of decomposed  : ', n_con
    !
  elseif (dft .eq. 'vasp')then
    !
    write(*,*) 'DFT mode              :         VASP '
    !
    file_read = 'PROCAR'
    !
    open( unit = 99 , file = file_read , action = 'read' )
    !
    read( 99 , '()')
    !
    read( 99 , *) temp_1 ,temp_2, temp_3, nks, &
                  temp_1 ,temp_2, temp_3, nbnd, &
                  temp_1 ,temp_2, temp_3, nat
    !
    write(*,*) 'Read from             :       ', trim( file_read )
    write(*,*) 'Number of wfc         : ', 9*nat
    write(*,*) 'Number of nk          : ', nks
    write(*,*) 'Number of band        : ', nbnd
    !
    ! allocatee vector with nk, nbnd, nat
    !
    natomwfc = 9 * nat
    !
    allocate( k(3,nks) )
    allocate( e(nbnd, nks) )
    allocate( sumproj(nbnd,nks) )
    allocate( selected(natomwfc ) )
    allocate( proj_vasp(natomwfc ) )
    !
    sumproj(:,:)  = 0.0
    selected( : )  = 0.0
    !
    do i = 1, n_con
      idum = orbital_limit(i)
      selected( idum ) = 1.0
    enddo
    !
    write(*,*) 'Number of decomposed  : ', n_con
    !
    do ink = 1, nks
      read(99, '()')
      read(99,'(a18, 3f11.8)') temp_1, k(1,ink), k(2,ink), k(3,ink)
      read(99,'()')
      !
      do ibnd = 1, nbnd
        read(99,*) temp_1, idum, temp_2, temp_3, &
                    e(ibnd,ink)
        read(99,'()')
        read(99,'()')
        !
        proj_vasp(:) = 0.0
        !
        do inat = 1, nat
          !
          read(99,*) idum, ( proj_vasp( (inat-1)*nat + i ), i =1,9 )
          !
          proj_vasp = dot_product(proj_vasp, selected)
          !
        enddo  ! end inat
        !
        sumproj(ibnd,ink) = sumproj(ibnd,ink) + sum(proj_vasp, dim=1)
        !
        read(99,'()')
        read(99,'()')
        !
      enddo  !  end ibnd
      !
    enddo  !  end ink
    !
    close(99)
    !
  else
    write(*,*) '**** Works only "qe" or "vasp" mode'
    STOP
  endif
  !
  !--------- end of reading standard input

  allocate( tmp(nks) )

  !--------- DO you need to calculate orbital decomposed band ?
  !------------- ldecom : 0 = NO
  !------------- ldecom : 1 = YES
  !--------- DO NOT NEED
  if (ldecom .eq. 0) then
    !
    sumproj(:,:) = 1.0
    write(*,*) 'WITHOUT DECOMPOSED FERMI SURFACE'
    ! ---------- sum over weight
    call sum_on_weight( e, sumproj, tmp, ef, nks, nbnd, thresh, degauss )

  !--------- NEED
  elseif (ldecom .eq. 1) then
    !
    write(*,*) 'ORBITAL DECOMPOSED FERMI SURFACE'
    !
    if (dft .eq. 'qe') then
      call fermi_surface( orbital_limit, selected, sumproj , n_con, nbnd,nks, natomwfc)
      ! ---------- sum over weight
      call sum_on_weight( e, sumproj, tmp, ef, nks, nbnd, thresh, degauss )
    endif
    if (dft .eq. 'vasp') then
      call sum_on_weight( e, sumproj, tmp, ef, nks, nbnd, thresh, degauss )
    endif
    !
  endif ! ldecom

  !
  ! ---------- writing FERMI surface output file
  write(*,*) 'writing Fermi surface output file to ' , trim(file_k)
  open(unit = 114  , file = file_k , action = 'write' )
  write(114, * ) '# kx  ky  kz  weight_on_fermi'
  do ink = 1, nks
    !
    write(114,'(3f12.6, f14.8)') (k(i , ink) , i=1,3) , tmp(ink)
    !
  enddo
  !
  !--------- DO you need to calculate SUSCEPTIBILITY ?
  !------------- lsuscep : 0 = NO
  !------------- lsuscep : 1 = YES
  !--------- DO NOT NEED
  if (lsuscep .eq. 0 ) then
    !
    write(*,*) 'Lindhard or nesting function will not be calculated'
  !--------- NEED
  elseif (lsuscep .eq. 1 ) then
    !
    call susceptibility( k, e, nbnd, nks,ef , sumproj)
    !
  endif ! lsuscep

end program


!------------------------------------------------ FUNCTIONS

!! Distribution function -- Fermi-deirac distribution

function wgauss (x, degauss)
  implicit none
  real(kind=8) :: x, degauss, wgauss
  real(kind=8), parameter :: maxarg = 200.d0
  !
  if (x.lt. - maxarg) then
     wgauss = 0.d0
  elseif (x.gt.maxarg) then
     wgauss = 1.d0
  else
     wgauss = 1.0d0 / (1.0d0 + exp ( - x/ degauss) )
  endif
  return
  !
end function wgauss

!! derivative of Fermi-dirac distribution function = delta function

function w0gauss (x, degauss)
  implicit none
  real(kind=8) :: w0gauss, x, degauss
  if (abs (x) .le. degauss*5 ) then
    w0gauss = 1.0d0 / (2.0d0 + exp ( - x/ degauss ) + exp ( + x/degauss) ) /degauss
  else
    w0gauss = 0.d0
  endif
  return
end function w0gauss


!------------------------------------------------ CALCULATE FERMI SURFACE


subroutine fermi_surface( orbital_limit, selected, sumproj, n_con, nbnd,nks, natomwfc )

  implicit none

  integer :: iwfc, ios, ink, nks, ibnd, nbnd,       &
             idum,fdum , n_con, natomwfc
  integer :: orbital_limit(n_con)
  real(kind=8) :: proj
  real(kind=8) :: sumproj(nbnd,nks) , selected(natomwfc )

  do iwfc = 1, maxval( orbital_limit(:) )
    !
    read(113, *,  IOSTAT=ios)
    !
    if (selected( iwfc ) .lt. 0.1) then
      !
      write(*,*) '   pass  reading iwfc : ' , iwfc
      !
      do ink = 1 , nks
        !
        do ibnd = 1, nbnd
          !
          read(113, *,  IOSTAT=ios)
          !
        enddo ! ibnd
        !
      enddo ! ink
      !
    else
      !
      write(*,*) '         READING iwfc : ' , iwfc , ' << CONSIDERED'
      !
      do ink = 1 , nks
        !
        do ibnd = 1, nbnd
          !
          read(113, '(2i8,f20.10)', IOSTAT=ios) idum,fdum,proj
          !
          if (idum/=ink.or.fdum/=ibnd.or.ios/=0) then
            !
            write(*,'("file with projections not compatible with bands")')
            STOP
            !
          else
            !
            sumproj(ibnd,ink) = sumproj(ibnd,ink) + proj * selected( iwfc )
            !
          endif
          !
        enddo ! ibnd
      !
      enddo ! ink

    endif
    !
  enddo ! iwfc

end subroutine fermi_surface


!------------------------------------------------ SUMMATION ON WEIGHT of ORBITALS


subroutine sum_on_weight( e, sumproj, tmp, ef, nks, nbnd, thresh, degauss )

  implicit none

  integer :: ink, nks, ibnd, nbnd
  real(kind=8)  :: thresh , ef, degauss
  real(kind=8)  :: tmp( nks ) , e(nbnd , nks ) , sumproj(nbnd,nks)
  real(kind=8), external :: wgauss, w0gauss

  tmp(:) = 0.0
  !
  do ink = 1, nks
    !
    do ibnd = 1, nbnd
      !
      tmp(ink) = tmp(ink) + &
                 w0gauss(e(ibnd , ink ) - ef , degauss) * sumproj(ibnd,ink)
      !
    enddo ! ibnd
    !
    if ( tmp(ink) .lt. thresh ) then
       tmp(ink) = 0.0
    endif
    !
  enddo  ! ink

end subroutine sum_on_weight


!------------------------------------------------ CALCULATE NESTING FUNCTION


subroutine susceptibility( k, e, nbnd, nks, ef , sumproj)

  implicit none

  integer :: ink, jnk , nks , count_tmp , nkx, nky , ibnd, nbnd,     &
         iinqx, inqx, iinqy, inqy, divide_q, divide_k, ibnd_start,   &
         ibnd_final , ndegauss, nq, inkx_tmp, inky_tmp,              &
         iinkx, inkx, iinky, inky, idegauss, jbnd, k_resolve
  real(kind=8) :: k( 3, nks) , e( nbnd, nks ), window, degauss_init, &
         delta_gauss, degauss, eigen_tmp_kk, eigen_tmp_kq, omega,    &
         diff_e, del_kk, del_kq, fd_kk, fd_kq, ef, sumproj(nbnd,nks)
  real(kind=8), allocatable ::  k_init( :,:,: ) , k_final(:,:,:) ,   &
                                eigen( :,:,: ), qpt(:, :, :),        &
                                suscep(:,:), weight_k(:,:),          &
                                int_suscep(:), sumproj_re(:,:,:),    &
                                k_resolved(:,:,:), suscep_k(:,:)
  real(kind=8), external :: wgauss, w0gauss
  character(len=80) :: file_q, file_conv, mode, file_real, file_imag
  complex(kind=8) :: ci, suscep_tmp
  NAMELIST / suscep_q / mode, ibnd_start, ibnd_final, ndegauss,      &
         degauss_init, delta_gauss, window, divide_k, divide_q,      &
         file_q, file_conv, omega, k_resolve

  !-------- default
  omega        = 0.001
  mode         = 'lindhard'
  k_resolve    = 0
  divide_q     = 1
  divide_k     = 1
  window       = 1
  ibnd_start   = 1
  ibnd_final   = nbnd
  ndegauss     = 10
  degauss_init = 0.005
  delta_gauss  = 0.005
  file_q       = 'SUSCEP.dat'
  file_conv    = 'CONVERGE.dat'
  !-----------------------
  !
  open(unit=111, file = 'input' , action = 'read' )
  read(111, suscep_q)
  close(111)
  !
  nkx = int(sqrt(real(nks)))
  nky = int(sqrt(real(nks)))

  allocate( k_init( nkx , nky , 3) )
  allocate( qpt( nkx , nky , 3) )
  allocate( eigen( nbnd , nkx , nky ) )
  allocate( sumproj_re(nbnd , nkx , nky) )

! reshape k_vector

  do ink = 1, nkx
    !
    do jnk = 1, nky
      !
      count_tmp = (ink -1 )*nky + jnk
      k_init( ink, jnk , : ) = k( : , count_tmp )
      do ibnd = 1, nbnd
        eigen(ibnd, ink, jnk )     = e(ibnd , count_tmp)
        sumproj_re(ibnd, ink, jnk) = sumproj(ibnd,count_tmp)
      enddo
      !
    enddo
    !
  enddo

  if (trim(mode) .eq. 'lindhard') then
    !
    write(*,*) "START calculating Lindhard function"
    !
    file_real = ""//trim(file_q)
    ! file_imag = "Im_"//trim(file_q)
    !
    open(unit=120 , file = file_real , action = 'write' )
    write(120,*) '# qx  qy  qz ', trim(mode), '_on_each_degauss[Re]'
    !
    ! open(unit=121 , file = file_imag , action = 'write' )
    ! write(121,*) '# qx  qy  qz ', trim(mode), '_on_each_degauss[Im]'
    !
  elseif (trim(mode) .eq. 'nesting') then
    !
    write(*,*) "START calculating Nesting function"
    !
    open(unit=120 , file = file_q , action = 'write' )
    write(120,*) '# qx  qy  qz ', trim(mode), '_on_each_degauss'
    !
  else
    !
    write(*,*) "only works for 'lindhard' and 'nesting' "
    STOP
    !
  endif

  allocate( suscep(2,ndegauss) )
  allocate( weight_k(nkx, nky) )
  allocate( int_suscep(ndegauss) )

  if (k_resolve .eq. 1 ) then
    !
    allocate( k_resolved(ndegauss, int(nkx / divide_k) ,int(nky / divide_k)) )
    allocate( suscep_k(2,ndegauss) )
    !
    k_resolved(:,:,:) = 0.0
    !
    if (trim(mode) .eq. 'lindhard') then
      write(*,*) '-- calculate k resolved Lindhard function '
    elseif (trim(mode) .eq. 'nesting') then
      write(*,*) '-- calculate k resolved Nesting function '
    endif
    !
  endif

  ci = (0.0, 1.0)
  int_suscep(:) = 0.0
  weight_k(:,:) = 1.0 / (nkx * nky / divide_k / divide_q )

  do iinqx = 1, nkx / divide_q
  !
  inqx = iinqx * divide_q
  do iinqy = 1, nky / divide_q
   !
   inqy = iinqy * divide_q
   !
   nq = (inqx - 1)*nky + inqy
   qpt(inqx, inqy, :) = k_init(inqx, inqy, :)
   suscep(:,:) = 0.0
   !
   do iinkx = 1, nkx / divide_k
     !
     inkx = iinkx * divide_k
     do iinky = 1, nky / divide_k
       !
       inky = iinky * divide_k
       !
       inkx_tmp = inkx + inqx - 1
       inky_tmp = inky + inqy - 1
       !
       ! periodic boundary condition
       !
       if (inkx_tmp .gt. nkx) then
          inkx_tmp = inkx_tmp - nkx
       endif
       if (inky_tmp .gt. nky) then
          inky_tmp = inky_tmp - nky
       endif
       !
       ! save information for k resolved result
       if (k_resolve .eq. 1 ) then
         suscep_k(:,:) = 0.0
       endif
       !
       if ( ( minval(abs((eigen(:,inkx,inky) - ef))) .lt. window ) .and. &
            ( minval(abs((eigen(:,inkx_tmp,inky_tmp) - ef))) .lt. window  )  ) then
         !
         do ibnd = ibnd_start, ibnd_final
           !
           do jbnd = ibnd_start, ibnd_final
            !
             do idegauss = 1, ndegauss
               !
               degauss = degauss_init + (idegauss-1)*delta_gauss
               !
               eigen_tmp_kk = eigen(ibnd, inkx, inky) - ef
               eigen_tmp_kq = eigen(jbnd, inkx_tmp, inky_tmp) - ef
               diff_e = eigen(ibnd, inkx, inky) -                  &
                        eigen(jbnd, inkx_tmp, inky_tmp)
               !
               ! delta function weights
               !
               del_kk = w0gauss( eigen_tmp_kk , degauss )
               del_kq = w0gauss( eigen_tmp_kq , degauss )
               !
               ! Fermi-dirac distribution weight
               !
               fd_kk = wgauss( eigen_tmp_kk , degauss )
               fd_kq = wgauss( eigen_tmp_kq , degauss )
               !
               ! susceptibility
               !
               if (trim(mode) .eq. 'lindhard') then
                 suscep_tmp = (fd_kk - fd_kq) / ( diff_e + ci*omega ) *  &
                               del_kk * del_kq *                         &
                               weight_k(inkx, inky)

                 suscep(1,idegauss) = suscep(1,idegauss)  +  real( suscep_tmp )
                 ! suscep(2,idegauss) = suscep(2,idegauss)  +  aimag( suscep_tmp )
                 if (k_resolve .eq. 1 ) then
                   suscep_k(1,idegauss) = suscep_k(1,idegauss)  +  real( suscep_tmp )
                 endif
                 !
                 !  nesting function
                 !
               elseif (trim(mode) .eq. 'nesting') then
                 suscep_tmp = del_kk * del_kq *                      &
                              weight_k(inkx, inky) *                 &
                              sumproj_re(ibnd, inkx, inky) *         &
                              sumproj_re(jbnd, inkx_tmp, inky_tmp)

                 suscep(1,idegauss) = suscep(1,idegauss)  +  real( suscep_tmp )
                 !
                 if (k_resolve .eq. 1 ) then
                   suscep_k(1,idegauss) = suscep_k(1,idegauss)  +  real( suscep_tmp )
                 endif
                 !
               endif
                !
               int_suscep(idegauss) = int_suscep(idegauss)    &
                            + suscep(1,idegauss) * weight_k(inqx,inqy)
               !
               if (k_resolve .eq. 1 ) then
                 k_resolved(idegauss , iinkx, iinky) = k_resolved(idegauss , iinkx, iinky) + &
                                       suscep_k(1,idegauss) * weight_k(inkx,inky)
               endif
               !
             enddo  ! idegauss
             !
           enddo  ! jbnd
           !
         enddo  ! ibnd
         !
       endif
       !
     enddo  ! inky
     !
   enddo  ! inkx
   !
   ! writing REAL part of Lindhard function
   !
   write(120, 8000) qpt(inqx, inqy ,1), qpt(inqx, inqy,2),qpt(inqx, inqy,3), &
                  (suscep(1,idegauss), idegauss = 1, ndegauss)
   !
   ! writing Imaginary part of Lindhard function
   !
   ! write(121, 8000) qpt(inqx, inqy ,1), qpt(inqx, inqy,2),qpt(inqx, inqy,3), &
                  ! (suscep(2,idegauss), idegauss = 1, ndegauss)
   !
  enddo  ! inqy
  !
  write(*,7000) iinqx, nky/divide_q
  !
  enddo  ! inqx
  !
  int_suscep(:) =  int_suscep(:)  / (nkx*nky /divide_q /divide_k)
  !
  ! writing file for convergence test
  !
  open(unit = 130 , file = file_conv , action = 'write')
  !
  write(130,*)  ' # degauss     SUM_on_q_sapce '
  !
  do idegauss = 1, ndegauss
    !
    degauss = degauss_init + (idegauss-1)*delta_gauss
    write(130,9060) degauss ,  int_suscep(idegauss)
    !
  enddo
  !
  ! writing file for k resolved result
  !
  if (k_resolve .eq. 1 ) then
    !
    open(unit = 222 , file = 'K_RESOLVED.dat' , action = 'write' )
    !
    write(222,*) '#  kx  ky  kz  k_resolved_on_each_degauss'
    !
    do iinkx = 1, nkx/divide_k
      !
      inkx = iinkx * divide_k
      !
      do iinky = 1, nky/divide_k
        !
        inky = iinky * divide_k
        !
        write(222, 8000) k_init(inkx, inky ,1), k_init(inkx, inky,2),k_init(inkx, inky,3), &
                       (k_resolved(idegauss , iinkx, iinky), idegauss = 1, ndegauss)
        !
      enddo
      !
    enddo
    !
  endif

  7000 format (" -- -- PROCESS :  q points number ", i5, "  from " , i8)
  8000 format (3f12.6 , 100f12.6)
  9000 format ( 100f12.6)
  9050 format ( 100i10)
  9060 format ( f9.6, "    " , f15.6 )

end subroutine susceptibility
