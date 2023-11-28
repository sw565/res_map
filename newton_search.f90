subroutine Newton_search(h_co,f)

  use bmad
  use ptc_interface_mod
!  use madx_ptc_module
  use pointer_lattice, dpe => dp
  
  implicit none
  type(c_vector_field)  f_rot,h_co
  type(c_ray) f,f1h
  real(dp) epsn0,epsn,epsnb
  integer nmax0,ncut,i,k
  logical potential_exit
  type(c_damap) id,rot
  complex(dp) dx

  ncut=c_%no+1 
  epsn0=abs(f%x(1))+abs(f%x(2))
  epsn0=epsn0*1.e-6_dp

  call alloc(f_rot)
  call alloc(id,rot)

  nmax0=100

  f1h=f     !  (6a) 
  potential_exit=.false.
  epsnb=1.e38_dp
  do k=1,nmax0

     id=1   ! (6b)
     do i=1,c_%nd2
        id%v(i)=id%v(i)+f1h%x(i)   ! (6c)
     enddo

     !  evaluating the vector field around the fixed point    
     do i=1,c_%nd2
        f_rot%v(i)=h_co%v(i).o.id   ! (6d)
     enddo
     ! putting in a c_damap because FPP can only inverse maps
     id=1
     do i=1,c_%nd2
        id%v(i)=f_rot%v(i)   
     enddo


     !    F(x)=0, the x is the fixed point
     !  .oo. does a TPSA inversion of the map, 
     !  ie, taking into account the constant part
     ! This inversion is exact if and only 
     ! if the calculation is done to infinite order
     ! Therefore this is part of a Newton search

     rot=1
     !do i=3,c_%nd2
     ! rot%v(i)=0.0_dp
     !enddo
     ! removing all "vertical dependence "  (6f)
     id=id.o.rot   
     ! adding identity in the y-py planes to allow inversion
     !do i=3,c_%nd2
     ! id%v(i)=1.0_dp.cmono.i ! (6g)  
     !enddo

     !doing a linear Newton search if ncut = 2
     id=id.cut.ncut   !  (6h) 
     id=id.oo.(-1)    !  (6i)    Inversion

     !   The fixed point is updated by taking the constant part of the map id.
     epsn=0
     do i=1, c_%nd2
        dx=(id%v(i).sub.'0')  !
        epsn=abs(dx)+epsn
        f1h%x(i)=f1h%x(i)+dx   ! (6j) updating the fixed point
     enddo

     if(potential_exit) then  
        if(abs(dx)==0.0_dp.or.epsn>=epsnb) exit
     endif
     if(epsn<epsn0) then
        potential_exit=.true.
     endif
     epsnb=epsn

     !     write(6,format8) (f1h%x(1:4))  ! check convergence

  enddo

  if(k>nmax0) then
     write(6,*) " Search was not succesful in", nmax0 , " steps"
  else 
!     write(6,*) " Search was succesful in", k, " steps"
  endif


  f=f1h
  call kill(f_rot)
  call kill(id,rot)


end subroutine Newton_search
