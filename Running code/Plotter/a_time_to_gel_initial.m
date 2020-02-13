%load('unified_data.mat');
i=310;%EDIT
run_array=[1:3 20:26 36:40];%EDIT
avg=0;
for k=run_array
    avg=avg+time_combined(i,k);
end
[~,length]=size(run_array);
avg=avg/length;