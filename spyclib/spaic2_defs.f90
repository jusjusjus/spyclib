! Written by Ulf Saalmann
! CONSTANTS
double precision, parameter :: pi=3.141592653589793238462d0
double precision, parameter :: dr=0.05, rmax=100, rasymp=50
integer, parameter :: nr=int(rmax/dr), nasymp=int(rasymp/dr)
double precision, parameter :: ecmin = 0.002d0

! GLOBAL PARAMS
double precision :: vclus=0.2d0, fluc=0.1d0, rpot=15d0
integer :: num_para=-1
