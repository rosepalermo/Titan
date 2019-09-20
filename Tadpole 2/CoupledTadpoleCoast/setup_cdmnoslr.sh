make clean
make main
mex -setup
mex -v -g in_environment.cpp visilibity.o
mex -v -g visibility_polygon.cpp visilibity.o