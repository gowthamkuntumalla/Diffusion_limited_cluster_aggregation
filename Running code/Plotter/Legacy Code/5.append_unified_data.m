clear;
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
mass_combined=mass_combined.';
time_combined=time_combined.';
mass_append = mass_combined(:,3);%EDIT column number for different runs
time_append = time_combined(:,3);%EDIT column number for different runs
% for k = [2:3 20:26 36:40]  %EDIT
%     mass_append=vertcat(mass_append,mass_combined(:,k)); % 15 runs now
%     time_append=vertcat(time_append,time_combined(:,k));
% end
filename='appended_unified_data.xlsx';
mass_inv_append=1./mass_append;
sheet = 1;
xlRange = 'A1';
xlswrite(filename,time_append,sheet,xlRange);
xlRange = 'B1';
xlswrite(filename,mass_inv_append,sheet,xlRange);
