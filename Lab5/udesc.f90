function ui(ii,xx)
implicit none
  integer, intent(in)                      :: ii
  integer                                  :: jj
  real(pr), dimension(1:N,1:3), intent(in) :: xx
  real(pr)                                 :: ui, sum_ui
  real(pr)                                 :: dx, dy, dz, r2
  sum_ui = 0._pr
  do jj = 1, N
    if (jj==ii) then
      cycle
    end if
    dx = xx(ii,1) - xx(jj,1)
    dy = xx(ii,2) - xx(jj,2)
    dz = xx(ii,3) - xx(jj,3)
    dx = dx - L*anint(dx/L)
    dy = dy - L*anint(dy/L)
    dz = dz - L*anint(dz/L)
    r2 = dx*dx + dy*dy + dz*dz
    if (r2==0._pr) then
      write(*,*) 'Hubo distancia cero entre las particulas', ii, jj
      cycle
    else if (r2>=r2cut) then
      cycle
    end if
    sum_ui = sum_ui + 1._pr/r2*r2*r2*r2*r2*r2 - 1._pr/r2*r2*r2
  end do
  ui = 4._pr*sum_ui
end function ui

function u(xx)
implicit none
  real(pr), dimension(1:N,1:3), intent(in) :: xx
  real(pr)                                 :: u, sum_u
  integer                                  :: ii
  sum_u = 0._pr
  do ii = 1, N
    sum_u = sum_u + ui(ii,xx)
  end do
  u = sum_u
end function u

function delta_u(ii,xx,xnew)
implicit none
  integer, intent(in)                      :: ii
  real(pr), dimension(1:N,1:3), intent(in) :: xx
  real(pr), dimension(1:N,1:3)             :: xxx
  real(pr), dimension    (1:3), intent(in) :: xnew
  real(pr)                                 :: delta_u
  real(pr)                                 :: u_ini, u_fin
  u_ini = ui(ii,xx)
  xxx(:,:) = xx(:,:)
  xxx(ii,:) = xx(ii,:) + xnew(:)
  u_fin = ui(ii,xxx)
  delta_u = u_fin - u_ini
end function delta_u
