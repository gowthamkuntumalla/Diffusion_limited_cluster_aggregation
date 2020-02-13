%%% Don't Forget to look at all "%EDIT" in this file form first line to last line 
function a_Plot_mass_spectrum( run_array,snapshot_for_mass_spec )
%Mass Spectrum plots
%**********Do it in Excel using binning***********%
% This part is used to get the data required for excel
load('unified_data.mat');
%%snapshot_for_mass_spec=350;%EDIT HERE
for k=run_array
    mass_var=sprintf('mass%d',k);
    %Rg_var=sprintf('Rg%d',k);
    time_var=sprintf('time%d',k); 
    mass_combineda=eval(mass_var); %This equation is a product of a mistake.But it works, Haha yikes!
    mass_combined{k,1} = mass_combineda{snapshot_for_mass_spec,1};
    time_combined{k,1} = eval(time_var);
end
conc_mass=zeros(0);
for k = run_array
   % [size_temp, foo]= size(mass_combined{k,1});
   conc_mass=vertcat(conc_mass,mass_combined{k,1});
end
filename=sprintf('unified_data_%d.xlsx',snapshot_for_mass_spec);
xlswrite(filename,conc_mass);%,time_combined);%% Use the unified_mass data at a given snapshot in file for binning 

end

