! ______________________________________________________________________________
!
! Programa time
!
! Este programa hace uso de los documentos de texto generados por el programa
! err.f90 para encontrar para cada método el conjunto de parámetros nx y nt que,
! dando un error global relativo menor a 10^-4, optimizan el tiempo de cpu.
!
! Gaston Briozzo, FaMAFyC, UNC, 27/04/2021
! ______________________________________________________________________________

program time
use precision, only: pr
implicit none

  integer  :: i, nx, nt
  integer  :: nx_fe, nt_fe, nx_imp, nt_imp, nx_cn, nt_cn
  real(pr) :: E_fe, E_imp, E_cn, E_fe_f, E_imp_f, E_cn_f
  real(pr) :: t_fe_i, t_imp_i, t_cn_i
  real(pr) :: t_fe_f, t_imp_f, t_cn_f
  real(pr) :: t_fe,   t_imp,   t_cn

! inicializamos los tiempos iniciales que nos serviran para comparar parámetros
  t_fe_i  = 1._pr
  t_imp_i = 1._pr
  t_cn_i  = 1._pr

! Abrimos los documentos de texto
  open(131,file='eg1800.d',status='old')
  open(132,file='et1800.d',status='old')

! Descartamo la primer línea de los documentos, ya que presenta los datos
  read(131,*)
  read(132,*)

! Leemos todas las líneas siguientes de los documentos
  do i = 2, 81405, 1

    read(131,*) nx, nt, E_fe, E_imp, E_cn
    read(132,*) nx, nt, t_fe_f, t_imp_f, t_cn_f

!   foward Euler
    if (E_fe<0.0001_pr) then  ! Solo consideraremos los parametrso que resulten en un error global relativo menor a 10^-4
      if (t_fe_f<t_fe_i) then ! Estamos buscando el tiempo óptimo, por lo que solo nos interesan los tiempos menores a los que ya hallamos obtenido
      ! si satisfacemos ambas condiciones, guardamos los parámetros  
        nx_fe = nx
        nt_fe = nt
        E_fe_f = E_fe
        t_fe   = t_fe_f
        t_fe_i = t_fe_f
      end if
    end if

!   implícito
    if (E_imp<0.0001_pr) then
      if (t_imp_f<t_imp_i) then
        nx_imp = nx
        nt_imp = nt
        E_imp_f = E_imp
        t_imp   = t_imp_f
        t_imp_i = t_imp_f
      end if
    end if

!   Crank-Nicolson
    if (E_cn<0.0001_pr) then
      if (t_cn_f<t_cn_i) then
        nx_cn = nx
        nt_cn = nt
        E_cn_f = E_cn
        t_cn   = t_cn_f
        t_cn_i = t_cn_f
      end if
    end if

  end do

! Escribimos en terminal los resultados
  write(*,*) nx_fe,  nt_fe,  E_fe_f,  t_fe
  write(*,*) nx_imp, nt_imp, E_imp_f, t_imp
  write(*,*) nx_cn,  nt_cn,  E_cn_f,  t_cn

end program
