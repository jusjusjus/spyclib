! Written by Ulf Saalmann
! CONSTANTS
double precision, parameter :: pi=3.141592653589793238462d0
double precision, parameter :: dr=0.05, rmax=100, rasymp=50
integer, parameter :: nr=int(rmax/dr), nasymp=int(rasymp/dr)
double precision, parameter :: ecmin = 0.05d0

! GLOBAL PARAMS
double precision :: rclus=13d0, drclus=1d0, vclus=0.222d0, fluc=0.25d0, rpot=20d0
integer :: num_para=-1
