%%First wrote on 7th july 2017
% This is used for getting slope of 1/Nc-1/Nc(0) vs time plot. It
% varies with time. it is equal to the Z exponent.
% 1/Nc-1/Nc(0) = constant*(time)^z is the kinetics equation
clear;
load('unified_data.mat');% Run other matalb fiels to get this file
num_data=418;%EDIT
run=3;%EDIT k for run number
N=11;% degree of polynomial fit

%% Getting data from the uploaded unified data
time_var=sprintf('time%d',run); 
mass_var=sprintf('mass%d',run); 
time = eval(time_var);
mass = eval(mass_var);
for snapshot= 1:418
    [Nc(snapshot,1),foo] = size(mass{snapshot,1});
end
Nc_inv=1./Nc;
Nc_inv_for_plot=Nc_inv-Nc_inv(1,1);% Here zero term 1/Nc(0) is subtracted
%% Getting slope and plotting
syms x y;
%Get trend line for the  1/Nc-1/Nc(0) vs time log log plot for this fv
% % [p,s] = polyfit(time,Nc_inv_for_plot,N);
% % for n = 1:N+1
% %     y = y + p(n) * x^n;
% % end
get trend line for y from excel curve fit of 1/Nc-1/Nc(0) vs time plot;

slope_var=diff(y,x);
x=log10(time);
slope=eval(slope_var);
plot(x,slope);
%or the next one
plot(time,slope);