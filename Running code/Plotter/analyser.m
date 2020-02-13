%% This analyser program contains two parts 
%%% Don't Forget to look at all "%EDIT" in this file form first line to last line 
% It uses get_the_full_data.m file 
clear;
numdata=418;%EDIT
fv='0.01';%EDIT
run_array=[1:3 20:26 36:40];%EDIT
%%% **The next function call and the **FOLLOWING save statement is required only once at the beginning of analysis

%a_get_the_full_data( numdata,fv,run_array);% function call to get the data in one place
%save('unified_data.mat');%COMMENT OUT this too if a_get_full_data is commented.
% IF any errors occur, Run the code with above two lines uncommented.

load('unified_data.mat');
snapshot_for_mass_spec=311;%EDIT in a_Plot_mass_spectrum file and this file. number for different snapshots (all runs combined)
run_for_slope=38;%EDIT number for different runs slope vs time plot
save('unified_data.mat');

%% append_unified_data into one variable
load('unified_data.mat');
for k = run_array 
        mass_var=sprintf('mass%d',k);
        time_var=sprintf('time%d',k); 
        mass_combined_matrix=eval(mass_var); %This equation is a product of a mistake.But it works, Haha yikes!
    for snapshot= 1:numdata  
        [Nc_combined(snapshot,k),~] = size(mass_combined_matrix{snapshot,1});% Nc matrix
        time_combined(:,k) = eval(time_var); %Time matrix for different folder numbers in rows
        % correspondence is present between these two variables.
    end
end
%% PART 1 %%
% It does slope vs time analysis
% % % This is used for getting slope of 1/Nc-1/Nc(0) vs time plot (and this
% % % plot itself)
% % % Slope varies with time. it is equal to the Z exponent.
% % % 1/Nc-1/Nc(0) = constant*(time)^z is the presumed kinetics equation
% % % Nc_to_plot=1/Nc-1/Nc(0);

%%%% 1byNc-1byNc(0) vs time plots
% for k = run_array 
%     loglog(time_combined(:,k),1./Nc_combined(:,k)-1./Nc_combined(1,k),'-s')%1/Nc-1/Nc(0) vs time plots
%     if(k==run_array(1))
%         hold on;
%         grid on;
%     end
% end
%%% printing the plot data into excel
% bar =1;
% for k=run_array
%      xlswrite('nc_inv_and_time.xlsx',1./Nc_combined(:,k)-1./Nc_combined(1,k),bar);
%      bar=bar+1;
% end
% xlswrite('nc_inv_and_time.xlsx',time_combined./eval(fv)^(-2.5),bar) % Printing time into last sheet i.e. bar(last)+1
%%%eval(fv)^2.5 is the scaling factor


%%%slope vs time plots
Nc_this_run= Nc_combined(:,run_for_slope);
time_this_run = time_combined(:,run_for_slope);
filename='appended_unified_data.xlsx';
Nc_inv=1./Nc_this_run;
sheet = 1;
xlRange = 'A1';
xlswrite(filename,time_this_run,sheet,xlRange);
xlRange = 'B1';
xlswrite(filename,Nc_inv,sheet,xlRange);


time_var=sprintf('time%d',run_for_slope);  
time = eval(time_var);
syms x y;
%%% Actually y=log10(1./Nc_this_run-1./Nc_this_run(1));
%get trend line for y from above excel curve fit of 1/Nc-1/Nc(0) vs time plot;
%/**THIS LINE IS TO BE UPDATED WITH POLYOMIAL FROM EXCEL**/
y = -0.057*x^3 + 1.741*x^2 - 7.440*x + 4.185; %EDIT 
slope_var=diff(y,x);
x=log10(time(2:end));
slope=eval(slope_var);
%plot(x,slope,'s');% Time is in log10 scale
%or the next one
plot(time(2:end),slope,'s'); %Here time is in normal scale
t_by_to=time/4032.5;%EDIT
%%NOTE: There is an asymptotic behaviour for the slope in the above plot

%% PART 2 %%
%Mass Spectrum plots
%**********Do it in Excel using binning***********%
% This part is used to get the data required for excel
%a_Plot_mass_spectrum(run_array,snapshot_for_mass_spec);% Carefully check here and change 
%snap_shot_for_mass_spec value in command window
