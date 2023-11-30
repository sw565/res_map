program res_map

  use bmad
  use sim_utils
  use m_cardano
  
  implicit none
  
  type(lat_struct), target :: lat
  type (ele_pointer_struct), allocatable :: eles(:)
  type (coord_struct), allocatable :: orbit(:)
  type (coord_struct) orb0
  type (normal_modes_struct) mode
  type (rad_int_all_ele_struct) rad_int
  
  real(rp) betax, alphax, gammax
  real(rp) betay, alphay, gammay
  real(rp) o_sfp(10,6), o_ufp(10,6)
  real(rp) j_sfp(10,2), j_ufp(10,2)
  real(rp) Qx_hz/0.0/, Qy_hz/0.0/, Qz_hz/0.0/, frev/390.136/, target_tunes(3)
  real(rp), allocatable :: dk1(:)
  real(rp) closed(6)
  
  character(255) file_name, lat_file
  character(40) :: regex_mask = ''

  logical set_tunes/.false./, ok, error,  use_nonlinear_twiss/.false./

  integer n_arg, ix_cache/0/, j, map_order, n_fixpoint

  real(rp) bx(10), ax(10), gx(10), by(10), ay(10), gy(10)
  real(rp) bx_u(10), ax_u(10), gx_u(10), by_u(10), ay_u(10), gy_u(10)
  real(rp) symp_err
  logical err_flag

 
  character(len=*), parameter :: fmt = '(3(1x,f8.5,sp,f8.5,"i"))'
  type(t_cubic_solution) :: p
  real(rp) jx_fp(2)

  
  lat_file = 'lat1.lat'
  map_order = 5
  n_fixpoint = 3
  
! input list
  namelist /input/ lat_file, set_tunes, Qx_hz, Qy_hz, Qz_hz, regex_mask, map_order, &
       use_nonlinear_twiss,  n_fixpoint

! read in the paramter file
  n_arg = cesr_iargc()
  if (n_arg > 1) then
    print *, 'Usage: trib_scan <input_file>'
    print *, 'Default: <input_file> = input.in'
    stop
  endif

  if (n_arg == 1) call cesr_getarg(1, file_name)
  if (n_arg == 0) then
    print '(a, $)', ' Input command file <CR=input.in>: '
    read (*, '(a)') file_name
  endif

  if (file_name == '')   file_name = 'input.in'
  print *, 'Opening: ', trim(file_name)
  open (unit= 1, file = file_name, status = 'old')
  read (1, nml = input)
  close (unit = 1)
  write(6, nml = input)

 
  call bmad_parser (lat_file,lat)
  print *,' bmad_parser done'

  frev = c_light/lat%param%total_length/1E3 ! in kHz
  call set_on_off(rfcavity$, lat, on$)
  call reallocate_coord (orbit, lat%n_ele_track)
  call closed_orbit_calc(lat,orbit,6)
  call track_all(lat,orbit)
  call lat_make_mat6(lat, -1, ref_orb = orbit) 
  call twiss_at_start(lat)
  call twiss_propagate_all(lat)
  call radiation_integrals (lat, orbit, mode, ix_cache, 0, rad_int)
  call calc_z_tune(lat%branch(0))

  orb0 = orbit(0)
  
  if (Qx_hz .eq. 0.0) then
     Qx_hz = lat%a%tune/twopi*frev
  endif

  if (Qy_hz .eq. 0.0) then
     Qy_hz = lat%b%tune/twopi*frev
  end if

  if (Qz_hz .eq. 0.0) then
     Qz_hz = lat%z%tune/twopi*frev
  end if

! set tunes
  if (set_tunes) then
     write (*,'(a21,f12.6,a9,f12.6,a9,f12.6)') 'Current tunes : Qx = ',lat%a%tune/twopi*frev, &
          '    Qy = ',lat%b%tune/twopi*frev,'    Qz = ', lat%z%tune/twopi*frev
     write (*,'(a21,f12.6,a9,f12.6,a9,f12.6)') 'Target tunes :  Qx = ',Qx_hz, '    Qy = ',Qy_hz,'    Qz = ', Qz_hz

     target_tunes(1)=(Qx_hz/frev+int(lat%ele(lat%n_ele_track)%a%phi/twopi))*twopi ! integer+fractional tunes
     target_tunes(2)=(Qy_hz/frev+int(lat%ele(lat%n_ele_track)%b%phi/twopi))*twopi 
     target_tunes(3)=Qz_hz/frev*twopi
       
     call choose_quads_for_set_tune(lat%branch(0), dk1, eles, regex_mask)
     ok = set_tune(target_tunes(1), target_tunes(2), dk1, eles, lat%branch(0), orbit)
     call set_z_tune(lat%branch(0), target_tunes(3))
  else
     print *, 'Use lattice design tunes!'
     print *, lat%a%tune/twopi, lat%b%tune/twopi,lat%z%tune/twopi
  end if

  call twiss_at_start(lat)
  call twiss_propagate_all(lat)
  call set_on_off(rfcavity$, lat, off$)

  write (*,'(a21,f12.6,a9,f12.6,a9,f12.6, a4)') 'Correct tunes :  Qx = ',lat%a%tune/twopi*frev, &
       '    Qy = ',lat%b%tune/twopi*frev,'    Qz = ', lat%z%tune/twopi*frev, ' kHz'   
  print '(2(a23,es12.4))', 'horizontal emittance = ',mode%a%emittance, &
       'Vertical emittance = ', mode%b%emittance

  betax=lat%ele(0)%a%beta
  gammax=lat%ele(0)%a%gamma
  alphax=lat%ele(0)%a%alpha
  betay=lat%ele(0)%b%beta
  gammay=lat%ele(0)%b%gamma
  alphay=lat%ele(0)%b%alpha
  write(6,*) " "
  write(6,*) "beta_x: ", betax, "beta_y:  ", betay
  write(6,*) "gamma_x: ", gammax, "gamma_y: ", gammay
  write(6,*) "alpha_x: ", alphax, "alpha_y: ", alphay
  write(6,*) " "
  
  select case (n_fixpoint)
  case (3)
     call ptc_calc_fps (lat, map_order, jx_fp, o_sfp, o_ufp, closed)
  case (4)
     call ptc_calc_fps (lat, map_order, jx_fp, o_sfp, o_ufp, closed)
  case (5)
     call ptc_calc_fps_5nux (lat, map_order, jx_fp, o_sfp, o_ufp, closed)
  end select
 
  write(6,*) "Stable fixed points: "     
  do j=1, n_fixpoint*2
     write(6,'(4e16.6)') o_sfp(j,1:4)
  end do

  write(6,*) "Unstable fixed points: "     
  do j=1, n_fixpoint*2
     write(6,'(4e16.6)') o_ufp(j,1:4)
  end do
  
  if (use_nonlinear_twiss) then
     ! Find the Twiss parameters at different fixed points
     do j=1, n_fixpoint*2
        orb0%vec=o_sfp(j,:)
        call twiss_from_tracking(lat, orb0, symp_err, err_flag)
        bx(j)=lat%ele(0)%a%beta
        gx(j)=lat%ele(0)%a%gamma
        ax(j)=lat%ele(0)%a%alpha
        by(j)=lat%ele(0)%b%beta
        gy(j)=lat%ele(0)%b%gamma
        ay(j)=lat%ele(0)%b%alpha
        print *,"symp_err:", symp_err
     end do

     do j=1, n_fixpoint*2
        orb0%vec=o_ufp(j,:)
        call twiss_from_tracking(lat, orb0, symp_err, err_flag)
        bx_u(j)=lat%ele(0)%a%beta
        gx_u(j)=lat%ele(0)%a%gamma
        ax_u(j)=lat%ele(0)%a%alpha
        by_u(j)=lat%ele(0)%b%beta
        gy_u(j)=lat%ele(0)%b%gamma
        ay_u(j)=lat%ele(0)%b%alpha
        print *,"symp_err:", symp_err
     end do

     write(6,*) " "
     do j=1, n_fixpoint*2
        write(6,*) j, "betax: ", bx(j), "gammax: ", gx(j), "alphax: ", ax(j)
     end do
     write(6,*) " "

     j_sfp(:,1) = (gx*(o_sfp(:,1)-closed(1))**2 + 2*ax*(o_sfp(:,1)-closed(1))*(o_sfp(:,2)-closed(2)) + bx*(o_sfp(:,2)-closed(2))**2)/2
     j_sfp(:,2) = (gy*(o_sfp(:,3)-closed(3))**2 + 2*ay*(o_sfp(:,3)-closed(3))*(o_sfp(:,4)-closed(4)) + by*(o_sfp(:,4)-closed(4))**2)/2

     j_ufp(:,1) = (gx_u*(o_ufp(:,1)-closed(1))**2 + 2*ax_u*(o_ufp(:,1)-closed(1))*(o_ufp(:,2)-closed(2)) + bx_u*(o_ufp(:,2)-closed(2))**2)/2
     j_ufp(:,2) = (gy_u*(o_ufp(:,3)-closed(3))**2 + 2*ay_u*(o_ufp(:,3)-closed(3))*(o_ufp(:,4)-closed(4)) + by_u*(o_ufp(:,4)-closed(4))**2)/2
  else
     
     ! Twiss parameters are different at different fixed points
     
     j_sfp(:,1) = (gammax*(o_sfp(:,1)-closed(1))**2 + 2*alphax*(o_sfp(:,1)-closed(1))*(o_sfp(:,2)-closed(2)) + betax*(o_sfp(:,2)-closed(2))**2)/2
     j_sfp(:,2) = (gammay*(o_sfp(:,3)-closed(3))**2 + 2*alphay*(o_sfp(:,3)-closed(3))*(o_sfp(:,4)-closed(4)) + betay*(o_sfp(:,4)-closed(4))**2)/2

     j_ufp(:,1) = (gammax*(o_ufp(:,1)-closed(1))**2 + 2*alphax*(o_ufp(:,1)-closed(1))*(o_ufp(:,2)-closed(2)) + betax*(o_ufp(:,2)-closed(2))**2)/2
     j_ufp(:,2) = (gammay*(o_ufp(:,3)-closed(3))**2 + 2*alphay*(o_ufp(:,3)-closed(3))*(o_ufp(:,4)-closed(4)) + betay*(o_ufp(:,4)-closed(4))**2)/2
  end if
  
  write(6,*) " "
  write(6,*) "Stable fixed point actions: "
  write(6,*) j_sfp(1:n_fixpoint*2,1)
  write(6,*) "Unstable fixed point actions: "
  write(6,*) j_ufp(1:n_fixpoint*2,1)
  write(6,*) " "

  
  ! first three are accurate fix points and last three are estimated 
  open(10, file="j_vs_tune.txt",action='write',position='append')
  write(10,'(14es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, j_sfp(1:n_fixpoint*2,1), j_ufp(1:n_fixpoint*2,1)
  write(6,'(14es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, j_sfp(1:n_fixpoint*2,1), j_ufp(1:n_fixpoint*2,1)
  close(10)

  ! PTC estimated actions of fixed points
  open(10, file="jx_est_map.txt",action='write',position='append')
  write(10,'(4es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, jx_fp
  write(6,'(4es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, jx_fp
  close(10)
     
  select case (n_fixpoint)
  case (3)  
     open(10, file="sfp_vs_tune.txt",action='write',position='append')
     write(10,'(30es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, closed(1:4), &
          o_sfp(1,1:4),o_sfp(2,1:4),o_sfp(3,1:4),o_sfp(4,1:4),o_sfp(5,1:4),o_sfp(6,1:4)
     close(10)

     open(10, file="ufp_vs_tune.txt",action='write',position='append')
     write(10,'(30es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, closed(1:4), &
          o_ufp(1,1:4),o_ufp(2,1:4),o_ufp(3,1:4),o_ufp(4,1:4),o_ufp(5,1:4),o_ufp(6,1:4)
     close(10)

  case (4)    
     open(10, file="sfp_vs_tune.txt",action='write',position='append')
     write(10,'(38es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, closed(1:4), &
          o_sfp(1,1:4), o_sfp(2,1:4), o_sfp(3,1:4), o_sfp(4,1:4), &
          o_sfp(5,1:4), o_sfp(6,1:4), o_sfp(7,1:4), o_sfp(8,1:4)
     close(10)

     open(10, file="ufp_vs_tune.txt",action='write',position='append')
     write(10,'(38es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, closed(1:4), &
          o_ufp(1,1:4), o_ufp(2,1:4), o_ufp(3,1:4), o_ufp(4,1:4), &
          o_ufp(5,1:4), o_ufp(6,1:4), o_ufp(7,1:4), o_ufp(8,1:4)
     close(10)

  case (5)
     open(10, file="sfp_vs_tune.txt",action='write',position='append')
     write(10,'(46es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, closed(1:4), &
          o_sfp(1,1:4), o_sfp(2,1:4), o_sfp(3,1:4), o_sfp(4,1:4), o_sfp(5,1:4), &
          o_sfp(6,1:4), o_sfp(7,1:4), o_sfp(8,1:4), o_sfp(9,1:4), o_sfp(10,1:4)
     close(10)

     open(10, file="ufp_vs_tune.txt",action='write',position='append')
     write(10,'(46es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, closed(1:4), &
          o_ufp(1,1:4), o_ufp(2,1:4), o_ufp(3,1:4), o_ufp(4,1:4), o_ufp(5,1:4), &
          o_ufp(6,1:4), o_ufp(7,1:4), o_ufp(8,1:4), o_ufp(9,1:4), o_ufp(10,1:4)
     
  end select
  
  open(10, file="twiss_vs_tune.txt",action='write',position='append')
  write(10,'(8es16.8)') lat%a%tune/twopi*frev,lat%b%tune/twopi*frev, alphax, betax, gammax, alphay, betay, gammay
  close(10)
  
end program res_map




