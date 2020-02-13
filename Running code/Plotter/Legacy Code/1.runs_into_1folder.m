cd('..');
clear;
for k = [1:3 20:26 36:40] 
cur_dir=cd;
foldername=sprintf('%s\\%d',cur_dir,k);
cd(foldername);
save('temp.mat');% save the variables
%get the data of mass and time
plotter; % note clear command is invoked in plotter
load('temp.mat'); %loads k and cur_dir
data_file_name=sprintf('%s_%d.mat','data_of_run',k);
save(data_file_name,'time','mass','Rg');
dest_dir=sprintf('%s\\%s',cur_dir,'fv=0.01');%EDIT
src_dir=sprintf('%s\\%s',foldername,data_file_name);
copyfile (src_dir, dest_dir);
delete temp.mat;
clear;
cd('..');
end
cur_dir=cd;
cd(sprintf('%s\\%s',cur_dir,'fv=0.01'));%EDIT