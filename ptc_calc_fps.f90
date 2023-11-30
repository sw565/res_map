subroutine ptc_calc_fps (lat, map_order, jx_fp, o_sfp, o_ufp, closed)

  use bmad
  use ptc_interface_mod
  use ptc_layout_mod
  use pointer_lattice, latptc => lat, dummy => pi, dummy1 =>twopi

  implicit none
  
  type(lat_struct), target :: lat
  integer map_order
  
  type(probe) xs0 
  type(probe_8) xs
  type(layout), pointer :: ring
  type(fibre),pointer:: p, ptc_fibre
  type(internal_state), target :: state
  
 !resonant map calculation
  type(c_damap) one_turn_map, id, a0, N_c
  type(c_normal_form) normal_form
  type(c_universal_taylor) H_res

  integer mf, mfo, i, j, fact
  complex(rp) del, alpha_xx_half, gamma
  
  real(rp) jx(2), discriminant
  type(c_ray) f1,f2
  type(c_ray) fix1,fix2
  real(rp) x1(6), x2(6)
  real(rp) g, phi0, phi

  real(rp) jx_fp(2), n_res/3.0/
  real(rp) o_sfp(10,6), o_ufp(10,6)
  real(rp) closed(6)
  
  !
 
  use_info = .true.
  use_quaternion=.true.
  call kanalnummer(mf,"fp_est.txt")
  call kanalnummer(mfo,"fp.txt")
  
  n_cai=-i_                             !  J= x_1 x_2 (h+ h-)    h+=(x - i px)/sqrt(2)
!  n_cai=-2*i_                          !  J= x_1 x_2 (h+ h-)/2  h+=(x - i px)

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

  !!  In PTC all objects that contain a TPSA variable 
  !!  or might contain one (..._8)  must be allocated and killed on exit
  !!
  call alloc(id,one_turn_map)
  call alloc(xs)
  call alloc(normal_form)  
  call alloc(a0)
  call alloc(N_c)

  
  xs0=closed ! probe =X(1:6)
  id=1       ! identity map
  xs=xs0+id  ! FPP -> PTC, closed orbit added

  call propagate(xs,state,fibre1=p)
  one_turn_map=xs ! PTC -> FPP
  
  call c_normal(one_turn_map, normal_form)
  call c_fast_canonise(normal_form%atot,a0)

! check whether it is 3nux or 4nux
  if ( abs(normal_form%tune(1)*3 - nint(normal_form%tune(1)*3)) .lt. 0.08) then
     n_res = 3.0_rp
  elseif( abs(normal_form%tune(1)*4 - nint(normal_form%tune(1)*4)) .lt. 0.1) then
     n_res = 4.0_rp
  end if
!  write(6,*) 'Resonance order: ', nint(n_res)
  
!  Resonant Map for the coefficients
  normal_form%nres=(c_%no+1)/n_res
  do i=1, (c_%no+1)/n_res
     normal_form%m(1,i)=n_res*i ! only leave resonant terms
  enddo

  call c_normal(one_turn_map,normal_form,canonize=.true.)
  call clean(normal_form%h_nl,normal_form%h_nl,prec=1.d-10)
  
  fact=nint(normal_form%tune(1)*n_res)  ! determines whether is near 1/3 (1/4) or 2/3 (3/4)
  
  normal_form%h_l=0
  normal_form%h_l%v(1)=( i_*twopi*fact/n_res)*dz_c(1) 
  normal_form%h_l%v(2)=(-i_*twopi*fact/n_res)*dz_c(2)
  normal_form%h_l%v(3)=( i_*twopi*normal_form%tune(2))*dz_c(3) 
  normal_form%h_l%v(4)=(-i_*twopi*normal_form%tune(2))*dz_c(4) 

  normal_form%h_l=ci_phasor()*normal_form%h_l
  N_c=normal_form%atot**(-1)*one_turn_map*normal_form%atot
  id=exp(normal_form%h_l)
  N_c=id*N_c
  N_c=ci_phasor()*N_c*c_phasor()

  normal_form%h=c_logf_spin(N_c)  
  call d_field_for_demin(normal_form%h, H_res)
  call clean(H_res,H_res,prec=1.d-5)
!  write(6,*)
!  write(6,*) " -2*pi*H_r "
!  write(6,*)
!  call print(H_res)

! Now calculate the fixed points:
  del=c_get_coeff(h_res,[1,1,0,0])
  alpha_xx_half= c_get_coeff(h_res,[2,2,0,0])  

  if (n_res .eq. 3.0_rp) then
     gamma= c_get_coeff(h_res,[3,0,0,0])
     phi0=atan2(aimag(gamma),real(gamma) )
     g=abs(gamma)

     discriminant=1.d0-8.d0*alpha_xx_half*del/9.d0/g**2
  
     jx(1)=-(3.d0*g/4.d0/alpha_xx_half)*(1.d0-sqrt(discriminant)) 
     jx(2)=(3.d0*g/4.d0/alpha_xx_half)*(1.d0+sqrt(discriminant))

     jx_fp = jx**2
     
  elseif (n_res .eq. 4.0_rp) then
     
     gamma= c_get_coeff(h_res,[4,0,0,0])
     phi0=atan2(aimag(gamma),real(gamma) )
     g=abs(gamma)
     
     jx(1)=-del/(2*alpha_xx_half + 4*g)
     jx(2)=-del/(2*alpha_xx_half - 4*g)

     jx_fp = jx

     if (jx(1) < 0 .or. jx(2) < 0) then
        print *, 'No islands in the phase space'
        goto 100
     end if
     
     jx(1) = sqrt(jx(1))
     jx(2) = sqrt(jx(2))
  
  end if

  write(6,'(a,2es15.5)') 'Estimated jx: ', jx_fp
  write(6,'(a,1es15.5,a,2es15.5,a)') ' axx: ', real(alpha_xx_half), ' g: (', real(gamma), aimag(gamma),')'
  write(6,'(a,1es15.5,a,1es15.5)') '|g|: ', abs(gamma), ' phi0: ', phi0
  write(6,*) 'Note axx and g need to be divided by a factor of -pi in order to compare with formula parameters!'
  write(6,'(a,1es15.5,a,1es15.5)') '|G|: ', abs(gamma)/pi, ' alpha_xx: ', -real(alpha_xx_half)/pi
  
  do j=1, nint(n_res)
     phi = twopi*(j-2)
     
     fix1=0
     fix2=0
     fix1%x(1)=jx(1)*exp(-i_*(phi0+phi)/n_res)
     fix1%x(2)=jx(1)*exp(i_*(phi0+phi)/n_res)
     fix2%x(1)=jx(2)*exp(-i_*(pi+phi0+phi)/n_res)
     fix2%x(2)=jx(2)*exp(i_*(pi+phi0+phi)/n_res)

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

     o_sfp(j+nint(n_res),:)=x1 ! estimate
     o_ufp(j+nint(n_res),:)=x2
     
     call Newton_search(normal_form%h,f1)    
     call Newton_search(normal_form%h,f2)

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

     write(mfo,format4) x1(1:4)
     write(mfo,format4) x2(1:4)

     o_sfp(j,:)=x1 ! exact
     o_ufp(j,:)=x2
     
  end do
  
100 continue
  
  call kill(id,one_turn_map)
  call kill(xs)
  call kill(normal_form)  
  call kill(a0)
  call kill(N_c)
  close(mf)
  close(mfo)
  
  call ptc_end(graphics_maybe=1,flat_file=.false.)
 
end subroutine ptc_calc_fps
