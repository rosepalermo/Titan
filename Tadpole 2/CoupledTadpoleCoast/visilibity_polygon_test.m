% To start debugging go the the terminal and type:
% /Applications/MATLAB_R2019a.app/bin/matlab -Dgdb
%
% Then follow the instructions here: https://www.mathworks.com/help/matlab/matlab_external/debugging-on-linux-platforms.html
% 
% (gdb) handle SIGSEGV SIGBUS nostop noprint
% (gdb) handle SIGUSR1 stop print
% (gdb) run -nojvm
%  >>   dbmex on
%  >>   visilibity_polygon_test
% (gdb) break visibility_polygon.cpp:152
% (gdb) break visilibity.cpp:2906
% (gdb) break visilibity.cpp:3121
% (gdb) c

load('visilibity_polygon_test.mat')
V = visibility_polygon(Pobs,env,epsilon,0.05);

% mex -setup C++
% mex -v in_environment.cpp visilibity.o
% mex -v shortest_path.cpp visilibity.o
% mex -v visibility_graph.cpp visilibity.o
% mex -v visibility_polygon.cpp visilibity.o