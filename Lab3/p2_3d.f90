! ______________________________________________________________________________
! CAMINATA ALEATORIA 3D
!
! Este programa simula una caminata aleatoria tridimensional, con probabilidad
! 1/6 de dar un paso unitario hacia adelante, atras, izda, drcha, arriba o abajo.
! Para obtener números al azar se emplean las subrutinas ran0, ran2, Marsaglia
! y Zaman y Mersenne Twister. Se simulan d2 caminatas de d1 pasos, cada una
! partiendo del origen, y se calculan trayectorias, desplazamiento cuadrático
! medio y estadías por cuadrante.
!
! Gaston Briozzo, FaMAFyC, UNC, 03/05/2021
! ______________________________________________________________________________

program p2_3d
use precision, only: pr
use randommod, only: ran0, ran2
use mzranmod,  only: rmzran
use mtmod,     only: grnd
implicit none

  integer, parameter                :: d1 = 10000, d2 = 1000000
  integer                           :: i, j, k, sem1, sem2
  integer                           :: xp, yp, zp
  integer, dimension(1:4)           :: x, y, z
  real(pr), dimension(1:4,1:d1)     :: dcm, dcmf
  real(pr), dimension(1:4)          :: u

! Definimos los formatos que utilizaremos en los documentos de texto
  1 format (13(A4, 5x))
  2 format (13(I7, 2x))
  3 format (A4, 5x, 4(A6, 15x))
  4 format (I7, 2x, 4(E23.16, 2x))

! Inicializamos los parámetros
  sem1 = 2073741823  ! Es primo
  sem2 = 1073741827  ! Es primo

! Desplazamientos cuadráticos medios
  dcmf(:,:) = 0._pr

! Abrimos los documentos de texto y describimos los datos a guardar
  open(321,file='cam.d')
  write(321,1) 'Paso' , 'x_r0', 'y_r0', 'z_r0', 'x_r2', 'y_r2', 'z_r2', 'x_mz', 'y_mz', 'z_mz', 'x_mt', 'y_mt', 'z_mt'
  open(322,file='dcm1.d')
  write(322,3) 'Paso', 'dcm_r0', 'dcm_r2', 'dcm_mz', 'dcm_mt'
  open(331,file='dcm.d')
  write(331,3) 'Paso', 'dcm_r0', 'dcm_r2', 'dcm_mz', 'dcm_mt'

! Repetimos la caminata d2 veces
  do i = 1, d2

!   Posición de partida
    x(:) = 0
    y(:) = 0
    z(:) = 0
!   Desplazamiento cuadrático en esta caminata
    dcm(:,:) = 0._pr

!   Calculamos el paso j de cada caminata
    do j = 1, d1

!     obtenemos los números u al azar (pseudo)
      u(1) = ran0(sem1)
      u(2) = ran2(sem2)
      u(3) = rmzran()
      u(4) = grnd()

! De acuerdo al número u obtenido, decidimos hacia donde irá nuestro caminante
      do k = 1, 4
      call paso(u(k),xp,yp,zp)
        x(k) = x(k) + xp
        y(k) = y(k) + yp
        z(k) = z(k) + zp
      end do

!     Calculamos los desplazamientos cuadráticos medios (dcm) y promediamos en el paso j, sobre las i caminatas
      dcm(:,j)  = real(x(:)*x(:) + y(:)*y(:) + z(:)*z(:),pr)
      dcmf(:,j) = dcmf(:,j) + dcm(:,j)

!     Si estamos en la primer caminata, la guardamos para tener una muestra no promediada
      if (i==1) then
        write(321,2) j, x(1), y(1), z(1), x(2), y(2), z(2), x(3), y(3), z(3), x(4), y(4), z(4)
        write(322,4) j, dcm(1,j), dcm(2,j), dcm(3,j), dcm(4,j)
      end if

    end do

!   Imprimimos i en la terminal para supervizar el avance del programa
    write(*,*) 'i=', i

  end do

! Normalizamos los desplazamientos cuadráticos medios
  dcmf(:,:) = dcmf(:,:)/real(d2,pr)

! Anotamos los datos obtenidos
  do j = 1, d1
    write(331,4) j,    dcm(1,j),    dcm(2,j),    dcm(3,j),    dcm(4,j)
  end do

CONTAINS

! ______________________________________________________________________________
! Subrutina paso
! Esta subrrutina seleccina, segun la entrada uu, hacia donde se dara el paso
! ______________________________________________________________________________

  subroutine paso(uu,xx,yy,zz)
  implicit none
    real(pr), intent(in) :: uu
    integer              :: xx, yy, zz
    xx = 0
    yy = 0
    zz = 0
    if (uu<=(1._pr/6._pr)) then
      xx = -1
    else if (uu<=(2._pr/6._pr)) then
      xx = +1
    else if (uu<=(3._pr/6._pr)) then
      yy = -1
    else if (uu<=(4._pr/6._pr)) then
      yy = +1
    else if (uu<=(5._pr/6._pr)) then
      zz = -1
    else
      zz = +1
    end if
  end subroutine paso

end program p2_3d
