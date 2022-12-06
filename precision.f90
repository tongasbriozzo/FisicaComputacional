module precision
implicit none

  integer, parameter :: sp = selected_real_kind(p=6)
  integer, parameter :: dp = selected_real_kind(p=13)
  integer, parameter :: tp = selected_real_kind(p=22)		  
  integer, parameter :: pr = dp

end module precision