clear;
%****Note: EDIT the following function before running them**^
%impNOTE
%runs_into_1folder;% get all the runs into one folder
%impNOTE
%runs_unified;% get the unified data
%impNOTE
load('unified_data.mat');
snapshot = 414;% the snapshot for mass spectrum etc
for k = [1:3 20:26 36:40]  %EDIT
    mass_var=sprintf('mass%d',k);
   %Rg_var=sprintf('Rg%d',k);
    time_var=sprintf('time%d',k); 
    mass_combineda=eval(mass_var); %This equation is a product of a mistake.But it works, Haha yikes!
    mass_combined{k,1} = mass_combineda{snapshot,1};
    time_combined{k,1} = eval(time_var);
end
conc_mass=mass_combined{1,1};%EDIT
for k = [2:3 20:26 36:40]
   % [size_temp, foo]= size(mass_combined{k,1});
   conc_mass=vertcat(conc_mass,mass_combined{k,1});
end
%Exportingthe concatenated mass for this snapshot files to Excel 
filename=sprintf('unified_data_%d.xlsx',snapshot);
xlswrite(filename,conc_mass);%,time_combined);


%%%% UPDATE : FOR NOW USE EXCEL FOR BINNING %%%%
%BINNING
A=conc_mass;
% botEdge=min(A); % X axis lower boundary
% topEdge=max(A); % X axis upper boundary
% numEdgesa=.15*size(A);% Make it less than size of N
% numEdges=floor(numEdgesa(1,1))+1;
% numBins=numEdges-1;
% binEdges=linspace(botEdge,topEdge,numEdges);% produce a equally spaced vector with size=numBins with the boundaries as shown
% binEdges=transpose(binEdges);
% binEdges=0:0.1:2001;
% binEdges=transpose(binEdges);
% [foo, numBins]=size(binEdges);
% [a_h,whichBin,bin] = histcounts(log10(A),binEdges);
numBins=50;
[bin_id,C]= kmeans(conc_mass,numBins);
for i = 1:numBins %+1 is to make sizes equal for plotting
    flagBinMembers = (bin_id == i); % collecting items of bin i (logic bit)
    binMembers     = A(flagBinMembers);% selects logic 1
    [bin_size(i),bar]=size(binMembers);
    binMeanWeight(i)  = mean(binMembers);%#ok<SAGROW> % take mean weight of each bin
end
%a_h=a_h.';
%a_h(end+1) = NaN;

% % MASS plots
loglog(C,(bin_size.').*(binMeanWeight.')/10000000);
   