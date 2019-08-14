t_max = 1000;
t = 1:t_max;%each time-step represents a Titan year
vol_L = zeros(1,t_max);%volume of methane in each reservoir
vol_K = zeros(1,t_max);
vol_P = zeros(1,t_max);
vol_S = zeros(1,t_max);
vol_A = zeros(1,t_max);
dep_L = zeros(1,t_max);%depth of methane in North and South reservoirs
dep_K = zeros(1,t_max);
dep_P = zeros(1,t_max);
dep_S = zeros(1,t_max);
SA_L = zeros(1,t_max);%surface area of methane in North and South reservoirs 
SA_K = zeros(1,t_max);
SA_P = zeros(1,t_max);
SA_S = zeros(1,t_max);
E_L = zeros(1,t_max);%evaporation rates, based on longitudinal radiation dependence
E_K = zeros(1,t_max);
E_P = zeros(1,t_max);
E_S = zeros(1,t_max);
P_L = zeros(1,t_max);%precipitation rates, also based on longitudinal dependence
P_K = zeros(1,t_max);
P_P = zeros(1,t_max);
P_S = zeros(1,t_max);
total_methane = 10000;
dep_L(1) = 10;%assumed to be constant (cylindrical)
SA_L(1) = 126000;%actual number can come from arcmap shapefiles
vol_L(1) = SA_L(1)*dep_L(1);%assumes cylinder
dep_K(1) = 20;
SA_K(1) = 100;
vol_K(1) = SA_K(1)*dep_K(1);
dep_P(1) = 20;
SA_P(1) = 100;
vol_P(1) = SA_P(1)*dep_P(1);
dep_S(1) = 20;
SA_S(1) = 100;
vol_S(1) = SA_S(1)*dep_S(1);
vol_A(1) = total_methane - vol_K(1) - vol_L(1) - vol_P(1) - vol_S(1);
E_L(1) = LG_evap; %comes from latitudinal averages calculation
E_L = ones(1,t_max)*LG_evap; %assumes that the surface area of LG, and therefore the corresponding evap and precip rates as well, stays constant
E_K(1) = 150;
E_P(1) = 150;
E_S(1) = 150;
P_L(1) = LG_precip; %also comes from latitudinal averages calculation
P_L = ones(1,t_max)*LG_precip; %like the earlier command for E_L = ones*LG_evap, this command will not apply once bathymetry data is incorporated
P_K(1) = 200;
P_P(1) = 200;
P_S(1) = 200;
for i = 2:t_max
    vol_L(i) = vol_L(i-1) - E_L(i-1) + P_L(i-1);
    vol_K(i) = vol_K(i-1) - E_K(i-1) + P_K(i-1);
    vol_P(i) = vol_P(i-1) - E_K(i-1) + P_P(i-1);
    vol_S(i) = vol_S(i-1) - E_S(i-1) + P_S(i-1);
    vol_A(i) = total_methane - vol_K(i) - vol_L(i) - vol_P(i) - vol_S(i);
    SA_L(i) = SA_L(1);%assuming cylindrical geometry (constant SA)_
    SA_K(i) = SA_K(1);
    SA_P(i) = SA_P(1);
    SA_S(i) = SA_S(1);
    dep_L(i) = vol_L(i)/SA_L(i);%assuming cylinder
    dep_K(i) = vol_K(i)/SA_K(i);
    dep_P(i) = vol_P(i)/SA_P(i);
    dep_S(i) = vol_S(i)/SA_S(i);
    % If not assuming cylinder and instead using lookup with a bathymetry
    % table, the depth variable becomes irrelevant and the surface area
    % variable will just be looked up on the table. Perhaps something like 
    % SA_L(i) = v(2,find(v(1,:)==vol_L(i)), where v is a lookup table with 
    % first row being volumes and second row being corresponding surface areas. 
    % If assuming perfect cylinder, precip and evap rates would remain
    % constant. If taking bathymetry into account, those rates would have
    % to be recalculated using longitudinal weighting for each iteration.
end


plot(t,vol_L)
xlabel('Time(Titan years)');
ylabel('Volume(m^3)');
title('Change in Volume of Ligeia Mare over time');
% hold on
% plot(t,P_N)
% legend('E','P')