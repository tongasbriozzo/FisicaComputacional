!*******************************************************************************
!
!  Programa VERSOR
!
!  Este programa esta diseñado para evaluar la calidad de nuestro generador de
!  versores aleatorios, obteniendo histogramas para su módulo, sus componentes
!  y su valor cuadrático medio, asi como su desviación estandar
!
!  Gaston Briozzo, Física Computacional, FaMAFyC, UNC, 07/07/2021
!
!*******************************************************************************

program versor
use precision, only: pr
use mtmod,     only: grnd
implicit none

  integer                         :: i, j, nh
  integer, parameter              :: nsim=1000000
  real(pr)                        :: n_mod, theta, phi
  real(pr)                        :: n1, n2, n4
  real(pr), dimension       (1:3) :: n
  real(pr), dimension(-50:49,1:3) :: his_c
  real(pr), dimension(-50:49)     :: his_m, his_2
  real(pr), parameter             :: pi=acos(-1._pr)

! Definimos los formatos
  101 format(A3,22x,A5,20x,A5,20x,A10,15x,A12)
  102 format(5(E23.16,2x))
  103 format(A5,20x,A12,13x,A12,13X,A12,13x,A6,19x,A18)
  104 format(6(E23.16,2x))

! Inicializamos las variables
  n1 = 0._pr          ! Valor medio del módulo del versor
  n2 = 0._pr          ! Valor medio del cuadrado módulo del versor
  n4 = 0._pr          ! Valor medio del (módulo del versor)^4
  his_c(:,:) = 0._pr  ! Histograma de las componentes del versor
  his_m(:)   = 0._pr  ! Histograma del módulo del versor
  his_2(:)   = 0._pr  ! Histograma del módulo al cuadrado del versor

! Generamos nsim versores independientes
  do i = 1, nsim

    n(1) = gasdev()
    n(2) = gasdev()
    n(3) = gasdev()

    n_mod = sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))

    n1 = n1 + n_mod
    n2 = n2 + n_mod*n_mod
    n4 = n4 + n_mod*n_mod*n_mod*n_mod

!   Realizamos el histograma sobre las componentes del versor
    do j = 1, 3
      do nh = -50, 49
        if (real(nh)<=10._pr*n(j).and.10._pr*n(j)<real(nh+1,pr)) then
          his_c(nh,j) = his_c(nh,j) + 1._pr
          cycle
        end if
      end do
    end do

!   Realizamos el histograma sobre el módulo del versor
    do nh = -50, 49
      if (real(nh)<=10._pr*n_mod.and.10._pr*n_mod<real(nh+1,pr)) then
        his_m(nh) = his_m(nh) + 1._pr
        cycle
      end if
    end do

!   Realizamos el histograma sobre el módulo al cuadrado del versor
    do nh = 0, 49
      if (real(nh)<=10._pr*n_mod*n_mod.and.10._pr*n_mod*n_mod<real(nh+1,pr)) then
        his_2(nh) = his_2(nh) + 1._pr
        cycle
      end if
    end do

  end do

! Normalizamos los valores obtenidos
  n1 = n1/real(nsim,pr)
  n2 = n2/real(nsim,pr)
  n4 = n4/real(nsim,pr)
  his_c = his_c/real(nsim,pr)
  his_m = his_m/real(nsim,pr)
  his_2 = his_2/real(nsim,pr)

! Anotamos los valores medios y sus varianzas
  open(201,file='med.d')
  write(201,101) '<n>', '<n^2>', '<n^4>', 'Varianza n', 'Varianza n^2'
  write(201,102) n1, n2, n4, n2-n1*n1, n4-n2*n2
  close(201)

! Anotamos los histogramas
  open(202,file='his.d')
  write(202,103) 'Valor', 'Componente x', 'Componente y', 'Componente z', 'Módulo', 'Módulo al cuadrado'
  do nh = -50, 49
    write(202,104) real(nh,pr)*0.1_pr+0.05_pr, his_c(nh,1), his_c(nh,2), his_c(nh,3), his_m(nh), his_2(nh)
  end do
  close(202)

CONTAINS

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

end program versor
