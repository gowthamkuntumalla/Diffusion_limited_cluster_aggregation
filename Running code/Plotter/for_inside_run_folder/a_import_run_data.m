function import_run_data( numfiles )
clearvars -except numfiles;
%%% Import time data
%numfiles=418;
time=zeros(numfiles,1);
for i=1:numfiles
    filename=sprintf('data_time%d',i);    
    time(i)=importdata(filename); 
end
%%% Import Proc_data
%Place this file in the folder where proc_data files are present
%plotter file comes after this and leverages this file's output
%Detailed explanation goes here
%numfiles=418;
mass=cell(numfiles,1); % Cell data struct
Rg=cell(numfiles,1);% radius of gyration
proc_data=cell(1,numfiles);

for i =1:numfiles
    filename=sprintf('proc_data%d',i);  
    FID=fopen(filename,'rt');
    proc_data{i}=textscan(FID,'%d %f %f %f %f %f %f %f','Headerlines',1); %#ok<SAGROW>
    mass{i,1}= proc_data{1,i}{1,2};
    Rg{i,1}= proc_data{1,i}{1,6};
    fclose(FID);
end
save('workspace_data'); %save workspace variables
end

