make clean
make main
module load mit/matlab/2018a
module add gcc/6.3.0
mex -setup
mex -v -g in_environment.cpp visilibity.o
mex -v -g visibility_polygon.cpp visilibity.o
mex inpoly.c