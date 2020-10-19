mex -I. -DWINDOWS=1 -c visilibity.cpp
mex -DWINDOWS=1 visilibity.obj in_environment.cpp
mex -DWINDOWS=1 visilibity.obj visibility_polygon.cpp