% mex all things that are needed to run tadpole & coastal erosion models

mex -setup

%visilibity
mex CoupledTadpoleCoast/in_environment.cpp CoupledTadpoleCoast/visilibity.o


mex mexD8.c
mex mexDinf.c
mex mexDms.c
mex mexGetMats.c
mex mexLandslide.c
mex mexSlope.c
mex inpoly.c
