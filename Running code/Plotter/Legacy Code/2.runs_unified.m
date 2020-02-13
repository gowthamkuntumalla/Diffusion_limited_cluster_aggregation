% data concatenator from different runs
% run this after runs_into_1folder
clear;
for k = [1:3 20:26 36:40] 
    data_file_name=sprintf('%s_%d.mat','data_of_run',k);
    load(data_file_name);
    mass_var=sprintf('mass%d',k);
    Rg_var=sprintf('Rg%d',k);
    time_var=sprintf('time%d',k);   
    assignin('base',mass_var,mass);% !!! this function helped a lot
    assignin('base',Rg_var,Rg);
    assignin('base',time_var,time);
    clear mass Rg time;
end
clear mass_var Rg_var time_var;
clear k data_file_name;
save('unified_data.mat');
clear
