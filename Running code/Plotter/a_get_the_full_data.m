%%% Don't Forget to look at all "%EDIT" in this file form first line to last line 
function a_get_the_full_data( numdata,fv,run_array )
% a_get_the_full_data 
%  This function collects data from all the runs folders and 
save('func_inputs.mat');
% %% Runs data into one folder
cd('..');
clearvars -except numdata fv run_array;
for k = run_array 
    cur_dir=cd;
    foldername=sprintf('%s\\%d',cur_dir,k);
    cd(foldername);
    save('temp.mat');% save the variables
    %get the data of mass and time
    a_import_run_data(numdata); % note clear command is invoked in plotter
    load('temp.mat'); %loads k and cur_dir
    load('workspace_data.mat');
    data_file_name=sprintf('%s_%d.mat','data_of_run',k);
    save(data_file_name,'time','mass','Rg');
    fv_fold=sprintf('fv=%s',fv);
    dest_dir=sprintf('%s\\%s',cur_dir,fv_fold);
    src_dir=sprintf('%s\\%s',foldername,data_file_name);
    copyfile (src_dir, dest_dir);
    delete temp.mat;
    clearvars -except numdata fv run_array;
    cd('..');
end

cur_dir=cd;
fv_fold=sprintf('fv=%s',fv);
cd(sprintf('%s\\%s',cur_dir,fv_fold));%Get back to the fv folder

%% runs_unified to a single .mat file
% data concatenator from different runs
% run this after runs_into_1folder
 clearvars -except numdata fv run_array;
 load('func_inputs.mat');
for k = run_array
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
end

