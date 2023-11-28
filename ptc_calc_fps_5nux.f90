subroutine ptc_calc_fps_5nux (lat, map_order, jx_fp, o_sfp, o_ufp, closed)

  use bmad
  use ptc_interface_mod
  use ptc_layout_mod
  use pointer_lattice, latptc => lat, dummy => pi, dummy1 =>twopi
  use m_cardano
  
  implicit none
  
  type(lat_struct), target :: lat
  integer map_order
!  complex(rp), optional :: alpha_xx, g_fpp
  
  type(probe) xs0 
  type(probe_8) xs
  type(layout), pointer :: ring
  type(fibre),pointer:: p, ptc_fibre
  type(internal_state), target :: state

  type(c_damap) one_turn_map, id
  type(c_normal_form) normal_form

  real(rp) closed(6)

  !resonant map calculation
  type(c_damap) N_c
  type(c_universal_taylor) H_res
  logical use_vector_field
  integer mf,mfo, i, j, fact
  complex(rp) del, alpha_xx_half, gamma
  
  real(rp) jx(2)
  type(c_ray) f1,f2
  type(c_ray) fix1,fix2
  real(rp) x1(6), x2(6)
  real(rp) g, phi0, phi

  real(rp) jx_fp(2), param(3)
  real(rp) o_sfp(10,6), o_ufp(10,6)

  ! for cubic solutions
  character(len=*), parameter :: fmt = '(3(1x,f8.5,sp,f8.5,"i"))'
  type(t_cubic_solution) :: cubic_sol
  real(rp) a1, a2, a3, a4
  real(rp) amp(2)
  
  !
  
  jx_fp(:) = 0.0_rp
  
  use_info = .true.
  use_quaternion=.true.
  call kanalnummer(mf,"fp_est.txt")
  call kanalnummer(mfo,"fp.txt")
  
  n_cai=-i_                             !  J= x_1 x_2 (h+ h-)    h+=(x - i px)/sqrt(2)
!  n_cai=-2*i_                          !  J= x_1 x_2 (h+ h-)/2  h+=(x - i px)

!  write(6,*) " (map_order must be >= 3)  map order = ", map_order
 
  call ptc_ini_no_append
  
  switch_to_drift_kick=.false.
  check_excessive_cutting=.false.
  
  call lat_to_ptc_layout (lat) 
  ring => lat%branch(0)%ptc%m_t_layout
  p=>ring%start

  closed=0
  my_estate=>state
  state=only_4d0 ! care only (x,px,y,py), ignore (z,pz). In PTC (5,6) defined as (e,t)
  
  call init(state,map_order,0)
  call find_orbit_x(closed,state,1.e-7_rp,fibre1=p)
!  print *, "Closed orbit (4): "
!  write(6,format4) closed(1:4) ! closed orbit
  
  !!  In PTC all objects that contain a TPSA variable 
  !!  or might contain one (..._8)  must be allocated and killed on exit
  !!
  call alloc(id,one_turn_map)
  call alloc(xs)
  call alloc(normal_form)  
  call alloc(N_c)

  
  xs0=closed ! probe =X(1:6)
  id=1       ! identity map
  xs=xs0+id  ! FPP -> PTC, closed orbit added

  call propagate(xs,state,fibre1=p)
  one_turn_map=xs ! PTC -> FPP
  
 ! Resonant Map for the coefficients
  
 ! normal_form%positive=.false.! for tunes conversion to positive or not
  
  normal_form%nres=(c_%no+1)/5
!  write(6,*) " # of resonance terms ",normal_form%nres
  do i=1, (c_%no+1)/5
     normal_form%m(1,i)=5*i ! only leave 1/5 order 
  enddo

  call c_normal(one_turn_map,normal_form,canonize=.true.)
!  write(6,*) normal_form%tune(1:3)
!  write(6,*) normal_form%spin_tune

  call clean(normal_form%h_nl,normal_form%h_nl,prec=1.d-10)
  
  fact=nint(normal_form%tune(1)*5)  ! this determines whether is near 3/5 or 4/5
!  print *, "factor: ", fact
  
  use_vector_field =.false.
  if(use_vector_field) then
     normal_form%h_l%v(1)=( i_*twopi*fact/5.0)*dz_c(1)+normal_form%h_l%v(1)
     normal_form%h_l%v(2)=(-i_*twopi*fact/5.0)*dz_c(2)+normal_form%h_l%v(2)
     normal_form%h_l%v(3)=0.0_rp
     normal_form%h_l%v(4)=0.0_rp
     N_c=exp(normal_form%h_l)
     N_c=exp(normal_form%h_nl,N_c)
  else
     normal_form%h_l=0
     normal_form%h_l%v(1)=( i_*twopi*fact/5.0)*dz_c(1) 
     normal_form%h_l%v(2)=(-i_*twopi*fact/5.0)*dz_c(2)
     normal_form%h_l%v(3)=( i_*twopi*normal_form%tune(2))*dz_c(3) 
     normal_form%h_l%v(4)=(-i_*twopi*normal_form%tune(2))*dz_c(4) 

     normal_form%h_l=ci_phasor()*normal_form%h_l
     N_c=normal_form%atot**(-1)*one_turn_map*normal_form%atot
     id=exp(normal_form%h_l)
     N_c=id*N_c
     N_c=ci_phasor()*N_c*c_phasor()
  endif

  normal_form%h=c_logf_spin(N_c)  ! this shows message "no convergence in c_logf_spin"

  call d_field_for_demin(normal_form%h, H_res)
  call clean(H_res,H_res,prec=1.d-5)
!  call print(H_res)
!  write(mf,*)
!  write(mf,*) " -2*pi*H_r "
!  write(mf,*)
!  call print(H_res,mf)


! Now calculate the fixed points:

  gamma= c_get_coeff(h_res,[5,0,0,0])  
  alpha_xx_half= c_get_coeff(h_res,[2,2,0,0])
  del=c_get_coeff(h_res,[1,1,0,0])
  
  phi0=atan2(aimag(gamma),real(gamma) )
  g=sqrt(aimag(gamma)**2+real(gamma)**2)
!  write(6,*) "del: ", del
!  write(6,*) "axx/2: ", alpha_xx_half
!  write(6,*) "g: ", g


  !computing the fixed points stable and unstable analytical using 5th order results
  a1=5/2*g
  a2=2*real(alpha_xx_half)
  a3=0.0_rp
  a4=del;

  amp=0

  cubic_sol = solve_cubic(a1, a2, a3, a4)

!  write(*, fmt) cubic_sol%x1, cubic_sol%x2, cubic_sol%x3
  
  if (aimag(cubic_sol%x1) .eq. 0 .and. real(cubic_sol%x1) > 0 .and. real(cubic_sol%x1) < 0.01 ) then
     amp(1) = real(cubic_sol%x1)
  elseif (aimag(cubic_sol%x2) .eq. 0 .and. real(cubic_sol%x2) > 0 .and. real(cubic_sol%x2) < 0.01 ) then
     amp(1) = real(cubic_sol%x2)
  elseif (aimag(cubic_sol%x3) .eq. 0 .and. real(cubic_sol%x3) > 0 .and. real(cubic_sol%x3) < 0.01 ) then
     amp(1) = real(cubic_sol%x3)
  end if

  cubic_sol = solve_cubic(-a1, a2, a3, a4)
!  write(*, fmt) cubic_sol%x1, cubic_sol%x2, cubic_sol%x3
  
  if (aimag(cubic_sol%x1) .eq. 0 .and. real(cubic_sol%x1) > 0 .and. real(cubic_sol%x1) < 0.01 ) then
     amp(2) = real(cubic_sol%x1)
  elseif (aimag(cubic_sol%x2) .eq. 0 .and. real(cubic_sol%x2) > 0 .and. real(cubic_sol%x2) < 0.01 ) then
     amp(2) = real(cubic_sol%x2)
  elseif (aimag(cubic_sol%x3) .eq. 0 .and. real(cubic_sol%x3) > 0 .and. real(cubic_sol%x3) < 0.01 ) then
     amp(2) = real(cubic_sol%x3)
  end if

  if (amp(1) .eq. 0 .or. amp(2) .eq. 0) then
     print *, 'No islands in the phase space'
     goto 100
  end if

!  print *, ' Amplitude of the fixed points (A=sqrt(J)):'
!  Write(6,*) amp

  jx = amp

  jx_fp = jx**2
  write(6,'(a,2es15.5, 2(a,es15.5))') " jx: ", jx_fp, " phi0: ", phi0, " g: ", g
     
!  optimize directly on alphaxx and g
  param(1) =real(alpha_xx_half)
  param(2) = g
  param(3) = phi0
  
  do j=1,5
     phi = twopi*(j-2)
     
     fix1=0
     fix2=0
     fix1%x(1)=jx(1)*exp(-i_*(phi0+phi)/5.d0)
     fix1%x(2)=jx(1)*exp(i_*(phi0+phi)/5.d0)
     fix2%x(1)=jx(2)*exp(-i_*(pi+phi0+phi)/5.d0)
     fix2%x(2)=jx(2)*exp(i_*(pi+phi0+phi)/5.d0)

     f1=0
     f1= fix1 
     f2=0
     f2= fix2

     ! transfer fixed points from normalized space to real space  for estimated points
     fix1=f1
     fix2=f2
     fix1=c_phasor().o.fix1
     fix2=c_phasor().o.fix2
     fix1=normal_form%atot.o.fix1
     fix2=normal_form%atot.o.fix2
     x1=0
     x2=0
     x1=fix1%x(1:6) + closed
     x2=fix2%x(1:6) + closed
     
     write(mf,format4) x1(1:4)
     write(mf,format4) x2(1:4)

     o_sfp(j+5,:)=x1 ! estimate
     o_ufp(j+5,:)=x2
     
!     write(6,*) 'Estimate SFPs in phasors ', j, ': '
!     write(6,format4) f1%x(1:2) 
     call Newton_search(normal_form%h,f1)
!     write(6,*) 'After Newton search FPs ', j, ': '
!     write(6,format4) f1%x(1:2)
     
!     write(6,*) 'Estimate UFPs in phasors ', j, ': '
!     write(6,format4) f2%x(1:2) 
     call Newton_search(normal_form%h,f2)
!     write(6,*) 'After Newton search FPs ', j, ': '
!     write(6,format4)f2%x(1:2)

     ! transfer fixed points from normalized space to real space  
     fix1=f1
     fix2=f2
     fix1=c_phasor().o.fix1
     fix2=c_phasor().o.fix2
     fix1=normal_form%atot.o.fix1
     fix2=normal_form%atot.o.fix2

     x1=0
     x2=0
     x1=fix1%x(1:6) + closed
     x2=fix2%x(1:6) + closed
     
!     write(6,*) " Fixed points of resonance "
!     write(6,format4) x1(1:4)
!     write(6,format4) x2(1:4)    
!     write(6,*) "  "
     write(mfo,format4) x1(1:4)
     write(mfo,format4) x2(1:4)
     
     o_sfp(j,:)=x1 ! exact
     o_ufp(j,:)=x2
     
  end do
  
100 continue
  
  call kill(id,one_turn_map)
  call kill(xs)
  call kill(normal_form)  
  call kill(N_c)


  goto 1001
  
1000 call ptc_end(graphics_maybe=1,flat_file=.false.)

1001 continue
  
end subroutine ptc_calc_fps_5nux
