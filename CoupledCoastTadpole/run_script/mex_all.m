% mex all things that are needed to run tadpole & coastal erosion models

mex -setup

%visilibity
mex Visilibity/in_environment.cpp Visilibity/visilibity.o
mex Visilibity/visibility_polygon.cpp Visilibity/visilibity.o


mex Tadpole/mexD8.c
mex Tadpole/mexDinf.c
mex Tadpole/mexDms.c
mex Tadpole/mexGetMats.c
mex Tadpole/mexLandslide.c
mex Tadpole/mexSlope.c
mex Coast/inpoly.c
