module m_cardano
  use bmad
  implicit none
  !Cardano's solutions:
  ! x1 = s + t - b/(3a) ,  real root
  ! x2 = -(s+t)/2 - b/(3a) + i(sqrt(3)/2)(s-t) , first complex root
  ! x3 = -(s+t)/2 - b/(3a) - i(sqrt(3)/2)(s-t) , second complex root
  !where:
  ! s = ( r + sqrt(q**3 + r**2) )**(1/3)
  ! t = ( r - sqrt(q**3 + r**2) )**(1/3)
  ! q = (3ac - b**2)/(9a**2)
  ! r = (9abc - 27da**2 - 2b**3)/(54a**3)
  private 
  public t_cubic_solution, solve_cubic


  type t_cubic_solution
     complex :: x1
     complex :: x2
     complex :: x3
  end type t_cubic_solution


contains
  pure type(t_cubic_solution) function solve_cubic(a, b, c, d) result(res)
    real(rp), intent(in)    :: a, b, c, d
    real(rp) q, r,  r_temp
    complex   :: s, t, c_temp, c_const
    real(rp) disc
    
    q       = (3.0*a*c - b**2)/(9.0*a**2) 
    r       = (9.0*a*b*c - 27.0*d*a**2 - 2.0*b**3)/(54.0*a**3)

    disc = q**3 + r**2
    if (disc .gt. 0) then
       r_temp    = r + sqrt(disc)
       s       = sign(1.0, r_temp) * abs(r_temp)**(1.0/3.0)
       r_temp    =  r - sqrt(disc)
       t       = sign(1.0, r_temp) * abs(r_temp)**(1.0/3.0)
    else
       c_temp = r + sqrt(cmplx(disc,0))
       s       = (c_temp)**(1.0/3.0)
       c_temp =  r - sqrt(cmplx(disc,0))
       t       = (c_temp)**(1.0/3.0)
    end if
    
    res%x1  = s + t - b/(3.0*a)
    c_const = cmplx(0, sqrt(3.0)/2.0)
    res%x2  = -(s+t)/2.0 - b/(3.0*a) + c_const*(s-t)
    res%x3  = -(s+t)/2.0 - b/(3.0*a) - c_const*(s-t)
    
  end function solve_cubic
end module m_cardano
