function import_proc_data(numfiles)
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
save('workspace_proc_data'); %save workspace variables
end





