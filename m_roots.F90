module roots
    implicit none
    
    contains
    ! *******************************************************************
    subroutine poleq4(n,c,t,m)
    ! solve quatric equations x^4+c(:,1)*x^3+c(:,2)*x^2+c(:,3)*x+c(:,4)=0
    ! Ferrari's solution
    ! Created by Cezar G. Diaconu - Tohoku University
    ! arguments:
    !   n - input  -  number of equations to be solved - integer
    !   c - input  -  coeficient matrix - size (n,4)   - real (8)
    !   c - output -  roots matrix      - size (n,4)   - real (8)
    !   t - output -  type of roots     - size (n,2)   - integer
    !   for example
    !       t(:,1)=1 - one real c(:,1)
    !       t(:,1)=2 - two real c(:,1) and c(:,2)
    !       t(:,1)=3 - complex  c(:,1) - real part
    !                           c(:,2) - imaginary part
    !   similar for t(:,2) and c(:,3), c(:,4)
    !   m - input  - mask the equations - size(n)      - logical
    integer, intent(in)                            :: n
    real(8), intent(inout), dimension(n,4), target :: c
    real(8), dimension(n,3)                        :: s3
    real(8), dimension(n,2)                        :: p,q
    real(8), dimension(:,:),pointer                :: cf,cff
    real(8), dimension(n)                          :: d1,d2,d
    integer, intent(out), dimension(n,2),target    :: t
    integer, dimension(:),pointer                  :: t1
    logical, dimension(n)                          :: lo
    logical, intent(in), dimension(n), optional    :: m
    real(8)                                        :: e,er
    e=0.0d0     ! error 1
    er=1.0d-10  ! error 2
    t1=>t(:,1)
    if(present(m)) then
     where(m)
      s3(:,1)=-c(:,2)
      s3(:,2)=c(:,1)*c(:,3)-4*c(:,4)
      s3(:,3)=4*c(:,2)*c(:,4)-c(:,3)**2-c(:,4)*(c(:,1)**2)
     end where
    else
     s3(:,1)=-c(:,2)
     s3(:,2)=c(:,1)*c(:,3)-4*c(:,4)
     s3(:,3)=4*c(:,2)*c(:,4)-c(:,3)**2-c(:,4)*(c(:,1)**2)
    endif
    if(present(m)) then
     call poleq3 (n,s3,t1,m)
    else
     call poleq3 (n,s3,t1)
    endif
    if (present(m)) then
     where (m)
      d=.25d0*(c(:,1)**2)-c(:,2)
      d1=d+s3(:,1)
      d2=(s3(:,1)*.5d0)**2-c(:,4)
      where (t1 == 2)
        where (d1 < e .or. d2 < e)
            d1=d+s3(:,2)
            d2=(s3(:,2)*.5d0)**2-c(:,4)
            d2=sqrt(abs(d2))
            q(:,1)=0.5d0*s3(:,2)-d2
            q(:,2)=0.5d0*s3(:,2)+d2
        elsewhere
            d2=sqrt(abs(d2))
            q(:,1)=0.5d0*s3(:,1)-d2
            q(:,2)=0.5d0*s3(:,1)+d2
        end where
      elsewhere (t1 == 3)
        where (d1 < e .or. d2 < e)
            d1=d+s3(:,2)
            d2=(s3(:,2)*.5d0)**2-c(:,4)
            where (d1 < e .or. d2 < e)
                d1=d+s3(:,3)
                d2=(s3(:,3)*.5d0)**2-c(:,4)
                d2=sqrt(abs(d2))
                q(:,1)=0.5d0*s3(:,3)-d2
                q(:,2)=0.5d0*s3(:,3)+d2
            elsewhere
                d2=sqrt(abs(d2))
                q(:,1)=0.5d0*s3(:,2)-d2
                q(:,2)=0.5d0*s3(:,2)+d2
            end where
        elsewhere
            d2=sqrt(abs(d2))
            q(:,1)=0.5d0*s3(:,1)-d2
            q(:,2)=0.5d0*s3(:,1)+d2
        end where
      elsewhere
        d2=sqrt(abs(d2))
        q(:,1)=0.5d0*s3(:,1)-d2
        q(:,2)=0.5d0*s3(:,1)+d2
      end where
      d1=sqrt(abs(d1))
      p(:,1)=0.5d0*c(:,1)-d1
      p(:,2)=0.5d0*c(:,1)+d1
     end where
    else
     d=.25d0*(c(:,1)**2)-c(:,2)
     d1=d+s3(:,1)
     d2=(s3(:,1)*.5d0)**2-c(:,4)
     where (t1 == 2)
        where (d1 < e .or. d2 < e)
            d1=d+s3(:,2)
            d2=(s3(:,2)*.5d0)**2-c(:,4)
            d2=sqrt(abs(d2))
            q(:,1)=0.5d0*s3(:,2)-d2
            q(:,2)=0.5d0*s3(:,2)+d2
        elsewhere
            d2=sqrt(abs(d2))
            q(:,1)=0.5d0*s3(:,1)-d2
            q(:,2)=0.5d0*s3(:,1)+d2
        end where
     elsewhere (t1 == 3)
        where (d1 < e .or. d2 < e)
            d1=d+s3(:,2)
            d2=(s3(:,2)*.5d0)**2-c(:,4)
            where (d1 < e .or. d2 < e)
                d1=d+s3(:,3)
                d2=(s3(:,3)*.5d0)**2-c(:,4)
                d2=sqrt(abs(d2))
                q(:,1)=0.5d0*s3(:,3)-d2
                q(:,2)=0.5d0*s3(:,3)+d2
            elsewhere
                d2=sqrt(abs(d2))
                q(:,1)=0.5d0*s3(:,2)-d2
                q(:,2)=0.5d0*s3(:,2)+d2
            end where
        elsewhere
            d2=sqrt(abs(d2))
            q(:,1)=0.5d0*s3(:,1)-d2
            q(:,2)=0.5d0*s3(:,1)+d2
        end where
     elsewhere
        d2=sqrt(abs(d2))
        q(:,1)=0.5d0*s3(:,1)-d2
        q(:,2)=0.5d0*s3(:,1)+d2
     end where
     d1=sqrt(abs(d1))
     p(:,1)=0.5d0*c(:,1)-d1
     p(:,2)=0.5d0*c(:,1)+d1
    endif
    cf=>c(:,1:2)
    cff=>c(:,3:4)
    if (present(m)) then
     where (m)
      lo=abs(p(:,1)*p(:,2)+q(:,1)+q(:,2) - c(:,2)) < er .and. &
         abs(p(:,1)*q(:,2)+p(:,2)*q(:,1) - c(:,3)) < er
      where (lo)
        cf(:,1)=p(:,1)
        cf(:,2)=q(:,1)
        cff(:,1)=p(:,2)
        cff(:,2)=q(:,2)
      elsewhere
        cf(:,1)=p(:,1)
        cf(:,2)=q(:,2)
        cff(:,1)=p(:,2)
        cff(:,2)=q(:,1)
      end where
     end where
    else
     lo=abs(p(:,1)*p(:,2)+q(:,1)+q(:,2) - c(:,2)) < er .and. &
        abs(p(:,1)*q(:,2)+p(:,2)*q(:,1) - c(:,3)) < er
     where (lo)
        cf(:,1)=p(:,1)
        cf(:,2)=q(:,1)
        cff(:,1)=p(:,2)
        cff(:,2)=q(:,2)
     elsewhere
        cf(:,1)=p(:,1)
        cf(:,2)=q(:,2)
        cff(:,1)=p(:,2)
        cff(:,2)=q(:,1)
     end where
    endif
    if (present(m)) then
      t1=>t(:,1)
      call poleq2(n,cf,t1,m)
      t1=>t(:,2)
      call poleq2(n,cff,t1,m)
    else
      t1=>t(:,1)
      call poleq2(n,cf,t1)
      t1=>t(:,2)
      call poleq2(n,cff,t1)
    endif
    end subroutine poleq4
    ! ****************************************************************
    subroutine poleq3(n,c,t,m)
    ! solve cubic equations x^3+c(:,1)*x^2+c(:,2)*x+c(:,3)=0
    ! Cardan's solution
    ! Created by Cezar G. Diaconu - Tohoku University
    ! arguments:
    !   n - input  -  number of equations to be solved - integer
    !   c - input  -  coeficient matrix - size (n,3)   - real (8)
    !   c - output -  roots matrix      - size (n,3)   - real (8)
    !   t - output -  type of roots     - size (n)     - integer
    !       t(:)=1 - one real c(:,1)
    !                two complex c(:,2) - real part
    !                            c(:,3) - imaginary part
    !       t(:)=2 - two real c(:,1) and c(:,2)
    !       t(:)=3 - three real c(:,1), c(:,2), c(:,3)
    !   m - input  - mask the equations - size(n)      - logical
    integer, intent(in)                      :: n
    real(8), intent(inout), dimension(n,3)   :: c
    real(8),                dimension(n)     :: q,p,qq,a,b
    integer, intent(out),   dimension(n)     :: t
    real(8)                                  :: pi3,e
    logical, intent(in), dimension(n), optional  :: m
    e=1.0d-16 ! error
    pi3=4.0d0*atan(1.0d0)/3
    if (present(m)) then
     where (m)
      p=-(c(:,1)**2)/3+c(:,2)
      c(:,1)=c(:,1)/3
      q=2*(c(:,1)**3)-c(:,1)*c(:,2)+c(:,3)
      qq=(p/3)**3+(q/2)**2
      where (qq > 0.0d0)
        qq=sqrt(qq)
        q=-.5d0*q
        a=q+qq
        where (a >= 0)
            a=a**(1.d0/3.0d0)
        elsewhere
            a=-(abs(a)**(1.d0/3.0d0))
        end where
        b=q-qq
        where (b >= 0)
            b=b**(1.d0/3.0d0)
        elsewhere
            b=-(abs(b)**(1.d0/3.0d0))
        end where
        c(:,3)=0.5d0*(a-b)*sqrt(3.0d0)
        c(:,2)=-0.5d0*(a+b)-c(:,1)
        c(:,1)=a+b-c(:,1)
        t = 1
      elsewhere (qq == 0.0d0)
        a=-.5*q
        where (a >= 0)
            a=a**(1.0d0/3.0d0)
        elsewhere
            a=-(abs(a)**(1d0/3.0d0))
        end where
        c(:,2)=-a-c(:,1)
        c(:,1)=2*a-c(:,1)
        t = 2
      elsewhere ! trigonometric solution (irreducible case)
        p=-p/3
        a=acos(-0.5d0*q/sqrt(p**3))/3.0d0
        p=2*sqrt(p)
        c(:,3)=-p*cos(a-pi3)-c(:,1)
        c(:,2)=-p*cos(a+pi3)-c(:,1)
        c(:,1)= p*cos(a)-c(:,1)
        t = 3
      end where
     end where
    else
     p=-(c(:,1)**2)/3+c(:,2)
     c(:,1)=c(:,1)/3
     q=2*(c(:,1)**3)-c(:,1)*c(:,2)+c(:,3)
     qq=(p/3)**3+(q/2)**2
     where (qq > 0.0d0)
        qq=sqrt(qq)
        q=-.5d0*q
        a=q+qq
        where (a >= 0)
            a=a**(1.d0/3.0d0)
        elsewhere
            a=-(abs(a)**(1.d0/3.0d0))
        end where
        b=q-qq
        where (b >= 0)
            b=b**(1.d0/3.0d0)
        elsewhere
            b=-(abs(b)**(1.d0/3.0d0))
        end where
        c(:,3)=0.5d0*(a-b)*sqrt(3.0d0)
        c(:,2)=-0.5d0*(a+b)-c(:,1)
        c(:,1)=a+b-c(:,1)
        t = 1
     elsewhere (qq == 0.0d0)
        a=-.5*q
        where (a >= 0)
            a=a**(1.0d0/3.0d0)
        elsewhere
            a=-(abs(a)**(1d0/3.0d0))
        end where
        c(:,2)=-a-c(:,1)
        c(:,1)=2*a-c(:,1)
        t = 2
     elsewhere ! trigonometric solution (irreducible case)
        p=-p/3
        a=acos(-0.5d0*q/sqrt(p**3))/3.0d0
        p=2*sqrt(p)
        c(:,3)=-p*cos(a-pi3)-c(:,1)
        c(:,2)=-p*cos(a+pi3)-c(:,1)
        c(:,1)= p*cos(a)-c(:,1)
        t = 3
     end where
    endif
    end subroutine poleq3
    ! ****************************************************************
    subroutine poleq2(n,c,t,m)
    ! solve quadratic equations x^2+c(:,1)*x+c(:,2)=0
    ! Created by Cezar G. Diaconu - Tohoku University
    ! arguments:
    !   n - input  -  number of equations to be solved - integer
    !   c - input  -  coeficient matrix - size(n,2)    - real(8)
    !   c - output -  roots matrix      - size(n,2)    - real(8)
    !   t - output -  type of roots     - size(n)      - integer
    !       t(:)=1 - one real c(:,1)
    !       t(:)=2 - two real c(:,1) and c(:,2)
    !       t(:)=3 - complex  c(:,1) - real part
    !                         c(:,2) - imaginary part
    !   m - input  - mask the equations - size(n)      - logical
    integer, intent(in)                    :: n
    real(8), intent(inout), dimension(n,2) :: c
    real(8),                dimension(n)   :: d
    real(8)                                :: e
    integer, intent(out),   dimension(n)   :: t
    logical, intent(in), dimension(n), optional  :: m
    e=1.0d-16 ! error
    if (present(m)) then
     where (m)
      d=c(:,1)**2-4*c(:,2)
      where (abs(d) < e)
        c(:,1) = - 0.5d0*c(:,1)
        c(:,2) = c(:,1)
        t = 1
      elsewhere (d > e)
        d=sqrt(abs(d))
        c(:,2) = 0.5d0*(- c(:,1) + d)
        c(:,1) = 0.5d0*(- c(:,1) - d)
        t = 2
      elsewhere
        c(:,1) = - 0.5d0*c(:,1)
        c(:,2) =   0.5d0*sqrt(abs(d))
        t = 3
      end where
     end where
    else
      d=c(:,1)**2-4*c(:,2)
      where (abs(d) < e)
        c(:,1) = - 0.5d0*c(:,1)
        c(:,2) = c(:,1)
        t = 1
      elsewhere (d > e)
        d=sqrt(abs(d))
        c(:,2) = 0.5d0*(- c(:,1) + d)
        c(:,1) = 0.5d0*(- c(:,1) - d)
        t = 2
      elsewhere
        c(:,1) = - 0.5d0*c(:,1)
        c(:,2) =   0.5d0*sqrt(abs(d))
        t = 3
      end where
    endif
    end subroutine poleq2
end module roots    
