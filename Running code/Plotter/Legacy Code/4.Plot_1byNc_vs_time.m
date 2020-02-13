clear;
% %part1 analysis
% % make sure following are run before the actual run
% %impNOTE
runs_into_1folder;% get all the runs into one folder
% %impNOTE
runs_unified;% get the unified data
% %impNOTE
%%%%%%%some times clear is needed
load('unified_data.mat');
for snapshot= 1:418
    for k = [1:3 20:26 36:40]  %EDIT folder numbers
        mass_var=sprintf('mass%d',k);
        time_var=sprintf('time%d',k); 
        mass_combineda=eval(mass_var); %This equation is a product of a mistake.But it works, Haha yikes!
        [mass_combined(k,snapshot),foo] = size(mass_combineda{snapshot,1});% Nc matrix
        time_combined(k,:) = eval(time_var); %Time matrix for different folder numbers in rows
        % Row correspondence is present between these two variables.
    end
end
clear time_var k foo mass_var snapshot;
color='rmbcywogB';
markers = '+o*.x-#!&';
time_combined=time_combined.';
mass_combined=mass_combined.';
onebyNc=1./mass_combined;
for k = [1:3 20:26 36:40]%EDIT
    loglog(time_combined(1:end,k),1./mass_combined(1:end,k),'-s')%,[markers(i) color(i)])
    if(k==1)%EDIT the first file name
        hold on;
        grid on;
    end
end

