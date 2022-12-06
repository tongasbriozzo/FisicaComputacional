!*******************************************************************************
!
!  Módulo de Dinámica Browniana
!
!  Este módulo contiene las funciones y subrutinas necesarias para el
!  funcionamiento de los programas de Dinámica Browniana
!
!  Gaston Briozzo, Física Computacional, FaMAFyC, UNC, 07/07/2021
!
!*******************************************************************************

module dbmod
use precision, only: pr
use mtmod,     only: grnd
implicit none

CONTAINS

!*******************************************************************************
! La funcion dcm_fi calcula el desplazamiento cuadratico medio de un sistema de
! N particulas entre las posiciones final xxf e inicial xxi
!*******************************************************************************
  function dcm_fi(NN,xxf,xxi)
  implicit none
    integer,                       intent(in) :: NN
    real(pr), dimension(1:NN,1:3), intent(in) :: xxf, xxi
    real(pr)                                  :: dcm_fi
    integer                                   :: ii
    real(pr)                                  :: dx, dy, dz, r2
    r2 = 0._pr
    do ii = 1, NN
      dx = xxf(ii,1) - xxi(ii,1)
      dy = xxf(ii,2) - xxi(ii,2)
      dz = xxf(ii,3) - xxi(ii,3)
      r2 = r2 + dx*dx + dy*dy + dz*dz
    end do
    dcm_fi = r2/real(NN,pr)
  end function dcm_fi

!*******************************************************************************
! La funcion u(xx) nos da la energia potencial de un sistema de N particulas
! distribuidas segun las posiciones xx bajo un potencial de Lennard-Jones
!*******************************************************************************
  function U(NN,LL,r2c,xx)
  implicit none
    integer,                       intent(in) :: NN
    real(pr),                      intent(in) :: LL, r2c
    real(pr), dimension(1:NN,1:3), intent(in) :: xx
    real(pr)                                  :: U
    integer                                   :: ii, jj
    real(pr)                                  :: dx, dy, dz, r2, sum_u
    sum_u = 0._pr
    do ii = 2, NN
      do jj = 1, ii-1
        dx = xx(jj,1) - xx(ii,1)
        dy = xx(jj,2) - xx(ii,2)
        dz = xx(jj,3) - xx(ii,3)
        dx = dx - LL*anint(dx/LL)
        dy = dy - LL*anint(dy/LL)
        dz = dz - LL*anint(dz/LL)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2==0._pr) then
!          write(*,*) 'Hubo distancia cero entre las particulas', ii, jj
          cycle
        else if (r2>=r2c) then
          cycle
        end if
        sum_u = sum_u + 1._pr/(r2*r2*r2*r2*r2*r2) - 1._pr/(r2*r2*r2)
      end do
    end do
    U = 4._pr*sum_u
  end function U

!*******************************************************************************
! La función vel calcula las velocidades en un sistema de N particulas dadas las
! posiciones final xxf, inicial xxi y el intervalo de tiempo dt
!*******************************************************************************
  function vel(NN,dtt,xxf,xxi)
  implicit none
    integer,                       intent(in) :: NN
    real(pr),                      intent(in) :: dtt
    real(pr), dimension(1:NN,1:3), intent(in) :: xxf, xxi
    real(pr), dimension(1:NN,1:3)             :: vel
    integer                                   :: ii
    do ii = 1, NN
      vel(ii,1) = (xxf(ii,1)-xxi(ii,1))/dtt
      vel(ii,2) = (xxf(ii,2)-xxi(ii,2))/dtt
      vel(ii,3) = (xxf(ii,3)-xxi(ii,3))/dtt
    end do
  end function vel

!*******************************************************************************
! La funcion K(vv) nos da la energia cinetica de un sistema de N particulas con
! velocidades dadas por el vector vv
!*******************************************************************************
  function K(NN,vv)
  implicit none
    integer,                       intent(in) :: NN
    real(pr), dimension(1:NN,1:3), intent(in) :: vv
    real(pr)                                  :: K
    integer                                   :: ii
    real(pr)                                  :: sum_k
    sum_k = 0._pr
    do ii = 1, NN
      sum_k = sum_k + vv(ii,1)*vv(ii,1) + vv(ii,2)*vv(ii,2) + vv(ii,3)*vv(ii,3)
    end do
    K = 0.5_pr*sum_k
  end function K

!*******************************************************************************
! La funcion fij nos da el vecor fuerza ejercida sobre la particula ii por la
! particula jj, bajo un potencial de Lennard-Jones
!*******************************************************************************
  function fij(NN,LL,r2c,ii,jj,xx)
  implicit none
    integer,                       intent(in) :: NN, ii, jj
    real(pr),                      intent(in) :: LL, r2c
    real(pr), dimension(1:NN,1:3), intent(in) :: xx
    real(pr), dimension     (1:3)             :: fij
    real(pr)                                  :: dx, dy, dz, r2
    if (jj==ii) then
      fij = 0._pr
      return
    end if
    dx = xx(ii,1) - xx(jj,1)
    dy = xx(ii,2) - xx(jj,2)
    dz = xx(ii,3) - xx(jj,3)
    dx = dx - LL*anint(dx/LL)
    dy = dy - LL*anint(dy/LL)
    dz = dz - LL*anint(dz/LL)
    r2 = dx*dx + dy*dy + dz*dz
    if (r2==0._pr) then
!      write(*,*) 'Hubo distancia cero entre las particulas', ii, jj
      fij = 0._pr
      return
    else if (r2>=r2c) then
      fij = 0._pr
      return
    end if
    fij(1) = 48._pr*dx*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    fij(2) = 48._pr*dy*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    fij(3) = 48._pr*dz*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    return
  end function fij

!*******************************************************************************
! La funcion fij_np es una modificacion de la funcion fij que no asume
! condiciones periodicas de contorno
!*******************************************************************************
  function fij_np(NN,r2c,ii,jj,xx)
  implicit none
    integer,                       intent(in) :: NN, ii, jj
    real(pr),                      intent(in) :: r2c
    real(pr), dimension(1:NN,1:3), intent(in) :: xx
    real(pr), dimension     (1:3)             :: fij_np
    real(pr)                                  :: dx, dy, dz, r2
    if (jj==ii) then
      fij_np = 0._pr
      return
    end if
    dx = xx(ii,1) - xx(jj,1)
    dy = xx(ii,2) - xx(jj,2)
    dz = xx(ii,3) - xx(jj,3)
    r2 = dx*dx + dy*dy + dz*dz
    if (r2==0._pr) then
      fij_np = 0._pr
      return
    else if (r2>=r2c) then
      fij_np = 0._pr
      return
    end if
    fij_np(1) = 48._pr*dx*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    fij_np(2) = 48._pr*dy*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    fij_np(3) = 48._pr*dz*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    return
  end function fij_np

!*******************************************************************************
! La funcion fi nos da la fuerza total ejercida sobre la particula ii en un
! sistema de NN particulas bajo un potencial de Lennard-Jones
!*******************************************************************************
  function fi(NN,LL,r2c,ii,xx)
  implicit none
    integer,                       intent(in) :: NN, ii
    real(pr),                      intent(in) :: LL, r2c
    real(pr), dimension(1:NN,1:3), intent(in) :: xx
    real(pr), dimension     (1:3)             :: fi
    integer                                   :: jj
    real(pr), dimension     (1:3)             :: sum_fi
    sum_fi(:) = 0._pr
    do jj = 1, NN
      sum_fi = sum_fi + fij(NN,LL,r2c,ii,jj,xx)
    end do
    fi = sum_fi
  end function fi

!*******************************************************************************
! La funcion fi_np es una modificacion de la funcion fi que no asume condiciones
! de contorno periodicas
!*******************************************************************************
  function fi_np(NN,r2c,ii,xx)
  implicit none
  integer,                       intent(in) :: NN, ii
  real(pr),                      intent(in) :: r2c
  real(pr), dimension(1:NN,1:3), intent(in) :: xx
  real(pr), dimension     (1:3)             :: fi_np
  integer                                   :: jj
  real(pr), dimension     (1:3)             :: sum_finp
    sum_finp(:) = 0._pr
    do jj = 1, NN
      sum_finp = sum_finp + fij_np(NN,r2c,ii,jj,xx)
    end do
    fi_np = sum_finp
  end function fi_np

!*******************************************************************************
! La funcion pf da la suma de los productos escalares entre las fuerzas que se
! ejercen las particulas ii y jj y la distancia entre estas
!*******************************************************************************
  function pf(NN,LL,r2c,xx)
  implicit none
    integer,                       intent(in) :: NN
    real(pr),                      intent(in) :: LL, r2c
    real(pr), dimension(1:NN,1:3), intent(in) :: xx
    real(pr)                                  :: pf
    integer                                   :: ii, jj
    real(pr), dimension     (1:3)             :: ff
    real(pr)                                  :: dx, dy, dz, sum_pf
    sum_pf = 0._pr
    do ii = 2, NN
      do jj = 1, ii-1
        dx = xx(ii,1) - xx(jj,1)
        dy = xx(ii,2) - xx(jj,2)
        dz = xx(ii,3) - xx(jj,3)
        dx = dx - LL*anint(dx/LL)
        dy = dy - LL*anint(dy/LL)
        dz = dz - LL*anint(dz/LL)
        ff = fij(NN,LL,r2c,ii,jj,xx)
        sum_pf = sum_pf + dx*ff(1) + dy*ff(2) + dz*ff(3)
      end do
    end do
    pf = sum_pf
  end function pf

!*******************************************************************************
! La funcion Temperatura nos da la temperatura cinetica de un sistema de NN
! particulas con velocidades dadas por el vector vv
!*******************************************************************************
  function Temperatura(NN,vv)
  implicit none
    integer,                       intent(in) :: NN
    real(pr), dimension(1:NN,1:3), intent(in) :: vv
    real(pr)                                  :: Temperatura
    integer                                   :: ii
    real(pr)                                  :: sum_t
    sum_t = 0._pr
    do ii = 1, NN
      sum_t = sum_t + vv(ii,1)*vv(ii,1) + vv(ii,2)*vv(ii,2) + vv(ii,3)*vv(ii,3)
    end do
    Temperatura = sum_t/real(3*NN,pr)
  end function Temperatura

!*******************************************************************************
! La funcion poc calcula el parametro de orden cristalino de un sistema de N
! particulas distribuidas segun el vetor xx para el vector kk
!*******************************************************************************
  function poc(NN,kk,xx)
  implicit none
    integer,                       intent(in) :: NN
    real(pr), dimension     (1:3), intent(in) :: kk
    real(pr), dimension(1:NN,1:3), intent(in) :: xx
    real(pr)                                  :: poc
    integer                                   :: ii
    real(pr)                                  :: sum_r, sum_i, kr
    sum_r = 0._pr
    sum_i = 0._pr
    do ii = 1, NN
      kr = kk(1)*xx(ii,1) + kk(2)*xx(ii,2) + kk(3)*xx(ii,3)
      sum_r = sum_r + cos(kr)
      sum_i = sum_i + sin(kr)
    end do
    poc = (sum_r*sum_r + sum_i*sum_i)/(NN*NN)
  end function poc

!*******************************************************************************
! La función n_ale devuelve un versor aleatorio tomado de una distribución
! gaussiana con <n>=0 y <n*n>=I, en una dirección aleatoria
!*******************************************************************************
  function n_ale()
  implicit none
    real(pr), dimension(1:3) :: n_ale
    n_ale(1) = gasdev()
    n_ale(2) = gasdev()
    n_ale(3) = gasdev()
  end function n_ale

!*******************************************************************************
! La subrutina maxwell toma un conjunto de N velocidades y las transforma para
! resultar una distribución de Maxwell con temperatura TT
!*******************************************************************************
  subroutine maxwell(NN,TT,VV)
  implicit none
    integer,                       intent(in)    :: NN
    real(pr),                      intent(in)    :: TT
    real(pr), dimension(1:NN,1:3), intent(inout) :: vv
    integer                                      :: ii
    do ii = 1, NN
      vv(ii,:) = n_ale()*sqrt(TT)
    end do
    return
  end subroutine maxwell

!*******************************************************************************
! La subrutina xcm toma la el conjunto de posicionesde un sistema de N
! particulas y modifica una al azar de modo que la velocidad resultante
! del centro de masa sea nula
!*******************************************************************************
  subroutine xcm(NN,xx)
  implicit none
    integer,                       intent(in)    :: NN
    real(pr), dimension(1:NN,1:3), intent(inout) :: xx
    integer                                      :: ii, nnxx
    real(pr), dimension     (1:3)                :: sum_xx
    sum_xx(:) = 0._pr
    do ii = 1, NN
      sum_xx(:) = sum_xx(:) + xx(ii,:)
    end do
    nnxx = 1 + int(real(NN-1,pr)*grnd())
    xx(nnxx,:) = xx(nnxx,:) - sum_xx(:)
    return
  end subroutine xcm

!*******************************************************************************
! la subrutina disrad da la funcion distribucion radial gg para un sistema de NN
! particulas distribuidas segun el vector xx
!*******************************************************************************
  subroutine disrad(NN,nnr,LL,pii,rhoo,xx,gg)
  implicit none
    integer,                       intent(in)  :: NN, nnr
    real(pr),                      intent(in)  :: LL, pii, rhoo
    real(pr), dimension(1:NN,1:3), intent(in)  :: xx
    real(pr), dimension(1:nnr),    intent(out) :: gg
    integer                                    :: ii, jj, kk
    real                                       :: dl, dx, dy, dz, dd, rmax, rmin
    gg(:) = 0._pr
    dl = 0.5_pr*LL/real(nnr,pr)
    do ii = 2, NN
      do jj = 1, ii-1
        dx = xx(ii,1) - xx(jj,1)
        dy = xx(ii,2) - xx(jj,2)
        dz = xx(ii,3) - xx(jj,3)
        dx = dx - LL*anint(dx/LL)
        dy = dy - LL*anint(dy/LL)
        dz = dz - LL*anint(dz/LL)
        dd = sqrt(dx*dx + dy*dy + dz*dz)
        do kk = 1, nnr
          if (real(kk-1,pr)*dl<dd.and.dd<=real(kk,pr)*dl) then
            gg(kk) = gg(kk) + 2._pr
          end if
        end do
      end do
    end do
    do kk = 1, nnr
      rmin = real(kk-1,pr)*dl
      rmax = real(kk,pr)*dl
      gg(kk) = gg(kk)/NN/(4._pr*pii*rhoo*(rmax*rmax*rmax-rmin*rmin*rmin)/3._pr)
    end do
    return
  end subroutine disrad

!*******************************************************************************
! La subrutina pbc aplica condiciones periodicas de contorno a un sistema de n
! particulas en una caja cubica de lado L
!*******************************************************************************
  subroutine pbc(NN,LL,xx)
  implicit none
    integer,                       intent(in)    :: NN
    real(pr),                      intent(in)    :: LL
    real(pr), dimension(1:NN,1:3), intent(inout) :: xx
    integer                                      :: ii, jj
    do ii = 1, NN
      do jj = 1, 3
        if (xx(ii,jj)<0._pr.or.LL<=xx(ii,jj)) then
          xx(ii,jj) = xx(ii,jj) - LL*floor(xx(ii,jj)/LL)
        end if
      end do
    end do
    return
  end subroutine pbc

!*******************************************************************************
! La subrutina langevin aplica el método de Langevin para obtener las posiciones
! de un sistema dinámico de N partículas en el instante siguiente.
!*******************************************************************************
  subroutine langevin(NN,LL,r2c,ffx,fnx,ffv,fnv,xx,xxnew,vvnew)
  implicit none
    integer,                       intent(in)  :: NN
    real(pr),                      intent(in)  :: LL, r2c, ffx, fnx, ffv, fnv
    real(pr), dimension(1:NN,1:3), intent(in)  :: xx
    real(pr), dimension(1:NN,1:3), intent(out) :: xxnew, vvnew
    integer                                    :: ii
    real(pr), dimension     (1:3)              :: ff, nv
    do ii = 1, NN
      ff = fi(NN,LL,r2c,ii,xx)
      nv = n_ale()
      xxnew(ii,:) = xx(ii,:) + ffx*ff(:) + fnx*nv(:)
      vvnew(ii,:) =            ffv*ff(:) + fnv*nv(:)
    end do
  end subroutine langevin

!*******************************************************************************
! La subrutina langevin_np es una variante de la subrituna langevin que no asume
! condiciones periódicas de contorno.
!*******************************************************************************
  subroutine langevin_np(NN,r2c,ffx,fnx,ffv,fnv,xx,xxnew,vvnew)
  implicit none
    integer,                       intent(in)  :: NN
    real(pr),                      intent(in)  :: r2c, ffx, fnx, ffv, fnv
    real(pr), dimension(1:NN,1:3), intent(in)  :: xx
    real(pr), dimension(1:NN,1:3), intent(out) :: xxnew, vvnew
    integer                                    :: ii
    real(pr), dimension     (1:3)              :: ff, nv
    do ii = 1, NN
      ff = fi_np(NN,r2c,ii,xx)
      nv = n_ale()
      xxnew(ii,:) = xx(ii,:) + ffx*ff(:) + fnx*nv(:)
      vvnew(ii,:) =            ffv*ff(:) + fnv*nv(:)
    end do
  end subroutine langevin_np

!*******************************************************************************
! La subrutina verlet implementa en metodo de Verlet para actualizarlas
! posiciones y velocidades de un sistema de N particulas.
!*******************************************************************************
  subroutine verlet(NN,LL,r2c,dt,xx,vv,xxnew,vvnew)
  implicit none
    integer,                       intent(in)  :: NN
    real(pr),                      intent(in)  :: LL, r2c, dt
    real(pr), dimension(1:NN,1:3), intent(in)  :: xx, vv
    real(pr), dimension(1:NN,1:3), intent(out) :: xxnew, vvnew
    integer                                    :: ii
    do ii = 1, NN
      xxnew(ii,:) = xx(ii,:) + vv(ii,:)*dt + 0.5_pr*fi(NN,LL,r2c,ii,xx)*dt*dt
    end do
    do ii = 1, NN
      vvnew(ii,:) = vv(ii,:) + 0.5_pr*(fi(NN,LL,r2c,ii,xxnew)+fi(NN,LL,r2c,ii,xx))*dt
    end do
    return
  end subroutine verlet

!*******************************************************************************
! La subrutina verlet_np es una variante de la subrituna verlet que no asume
! condiciones periódicas de contorno.
!*******************************************************************************
  subroutine verlet_np(NN,r2c,dt,xx,vv,xxnew,vvnew)
  implicit none
    integer,                       intent(in)  :: NN
    real(pr),                      intent(in)  :: r2c, dt
    real(pr), dimension(1:NN,1:3), intent(in)  :: xx, vv
    real(pr), dimension(1:NN,1:3), intent(out) :: xxnew, vvnew
    integer                                    :: ii
    do ii = 1, NN
      xxnew(ii,:) = xx(ii,:) + vv(ii,:)*dt + 0.5_pr*fi_np(NN,r2c,ii,xx)*dt*dt
    end do
    do ii = 1, NN
      vvnew(ii,:) = vv(ii,:) + 0.5_pr*(fi_np(NN,r2c,ii,xxnew)+fi_np(NN,r2c,ii,xx))*dt
    end do
    return
  end subroutine verlet_np

!*******************************************************************************
! La subrutina SC acomoda las posiciones de un sistema de N particulas para
! obtener una red cúbica simple (SC)
!*******************************************************************************
  subroutine SC(NN,LL,nll,xx)
  implicit none
    integer,                       intent(in)    :: NN, nll
    real(pr),                      intent(in)    :: LL
    real(pr), dimension(1:NN,1:3), intent(inout) :: xx
    integer                                      :: ii, jj, kk
    real(pr)                                     :: lado
    lado = LL/real(nll,pr)
    do ii = 1, nll
      do jj = 1, nll
        do kk = 1, nll
          xx(kk + nll*(jj-1) + nll*nll*(ii-1), 1) = lado*real(ii-1,pr)
          xx(kk + nll*(jj-1) + nll*nll*(ii-1), 2) = lado*real(jj-1,pr)
          xx(kk + nll*(jj-1) + nll*nll*(ii-1), 3) = lado*real(kk-1,pr)
        end do
      end do
    end do
    return
  end subroutine SC

!*******************************************************************************
! La subrutina FCC acomoda las posiciones de un sistema de N particulas para
! obtener una red cubica centrada en las caras (FCC)
!*******************************************************************************
  subroutine FCC(NN,LL,nll,xx)
  implicit none
    integer,                       intent(in)    :: NN, nll
    real(pr),                      intent(in)    :: LL
    real(pr), dimension(1:NN,1:3), intent(inout) :: xx
    integer                                      :: ii, jj, kk
    real(pr)                                     :: lado
    lado = LL/real(2*nll,pr)
    do ii = 1, 2*nll, 2
      do jj = 1, 2*nll, 2
        do kk = 1, 2*nll, 2
          xx(1 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 1) = lado*real(ii  ,pr)
          xx(1 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 2) = lado*real(jj  ,pr)
          xx(1 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 3) = lado*real(kk  ,pr)
          xx(2 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 1) = lado*real(ii  ,pr)
          xx(2 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 2) = lado*real(jj+1,pr)
          xx(2 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 3) = lado*real(kk+1,pr)
          xx(3 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 1) = lado*real(ii+1,pr)
          xx(3 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 2) = lado*real(jj  ,pr)
          xx(3 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 3) = lado*real(kk+1,pr)
          xx(4 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 1) = lado*real(ii+1,pr)
          xx(4 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 2) = lado*real(jj+1,pr)
          xx(4 + 4*(kk-1)/2 + 4*nll*(jj-1)/2 + 4*nll*nll*(ii-1)/2, 3) = lado*real(kk  ,pr)
        end do
      end do
    end do
    return
  end subroutine FCC

!*******************************************************************************
! La funcion gasdev, extraida del Numerical Recipes, da una distribucion de
! velocidades de Maxwell-Boltzmann
!*******************************************************************************
  function gasdev()
  implicit none
    real(pr) :: gasdev
    integer  :: iset
    real(pr) :: fac, gset, rsq, v1, v2, ran1
    SAVE        iset, gset
    DATA        iset/0/
    if (iset/=1) iset = 0
    if (iset/=0) iset = 1
    if (iset.eq.0) then
1     v1 = 2._pr*grnd()-1._pr
      v2 = 2._pr*grnd()-1._pr
      rsq = v1**2 + v2**2
      if(rsq.ge.1..or.rsq.eq.0.)goto 1
      fac = sqrt(-2._pr*log(rsq)/rsq)
      gset = v1*fac
      gasdev = v2*fac
      iset = 1
    else
      gasdev = gset
      iset = 0
    endif
    return
  end function gasdev

end module dbmod
