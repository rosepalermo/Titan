function [lake,dx] = test_square_repelem(ncells,sidelength)

lake1cell = zeros(3,3);
lake1cell(2,2) = 1;
lake = repelem(lake1cell,ncells,ncells);
dx = sidelength/ncells;


end