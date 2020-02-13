% Note Newer codes use these codes for analysis
% Written by Gowtham Kuntumalla, IIT Bombay in June 2017
% Project: Aggregation kinetics of Gelation, WUSTL
% This code is used for plotting results as defined in each PART
% First get the individual snapshot masses and Rg values to run this code.
% Use code import_proc_data.m for the above purpose
% look at file workspace_var.mat. Workspace variables should be loaded
% % Important terms:
% % Snapshot = An instant of time during simulation
% % 1.	Number of monomers initially placed, Nm
% % 2.	Volume fraction( or initial monomer packing density in the simulation box), fv
% % 3.  Number of clusters at a given time instant, Nc
% % 4.	Number of monomers in a given cluster, N
% % 5.	Average number of monomers (or mass) per cluster at a give time instant in the system, <N>
% % 6.	Radius of gyration of a given cluster, Rg
% % 7.	Average radius of gyration in a snapshot, <Rg> - This gives an idea about how the system is evolving

% In paper_runs proc_data'i' goes from 1 (Nc~50,000) and 2 to 18 in 500 steps(i.e.Nc= ~10,000 to ~2000)
% and 19 to 418 in 5 steps (i.e. till Nc=1)
%** ctrl-r and ctrl-t to comment , uncomment lines for specific plots
clear;
numfiles=418;
num_monomers=1000000;
Nm=num_monomers;
import_proc_data(numfiles);% function call
import_time(numfiles);
load('workspace_proc_data.mat'); % load the workspace variables
load('data_time.mat');

%############################################################%
%% PART 1 (It gives Aggregate monomer count N vs aggregate radius
% of gyration Rg normalized by monomer radius a on log-log plots.)
%EDITS%
%look at the mass, rad_gyr variable to get an idea of what to edit
%it is the row term in {THIS,..}
% rad_gyr=Rg{1,1};% edit this for different regimes
% mass_append=mass{1,1}; % edit this for different regimes
% for i=2:418% edit this too
%     rad_gyr=vertcat(rad_gyr,Rg{i,1}); %#ok<*AGROW>
%     mass_append=vertcat(mass_append,mass{i,1});
% end
%scatter(log10(rad_gyr*2),log10(mass_append)); % The plot

%############################################################%
%% PART 2 (For Frequency and weight in single snapshot of system)
%%EDIT%%
%%it is the row term in {THIS,}

% A=(mass{350,1}); % EDIT snapshot number, It has all clusters' mass info
% 
% % BINNING
% botEdge=min(A); % X axis lower boundary
% topEdge=max(A); % X axis upper boundary
% numEdgesa=.30*size(A);% Make it less than size of N
% numEdges=floor(numEdgesa(1,1))+1;
% numBins=numEdges-1;
% binEdges=linspace(botEdge,topEdge,numEdges);% produce a equally spaced vector with size=numBins with the boundaries as shown
% binEdges=transpose(binEdges);
% [h,whichBin,bin] = histcounts(A,binEdges);
% for i = 1:numBins+1 %+1 is to make sizes equal for plotting
%     flagBinMembers = (bin == i); % collecting items of bin i (logic bit)
%     binMembers     = A(flagBinMembers);% selects logic 1
%     binMeanWeight(i,1)  = mean(binMembers);%#ok<SAGROW> % take mean weight of each bin
% end
% h=h.';
% h(end+1) = NaN;
% 
% % % MASS plots
% scatter(binEdges,(h.*binMeanWeight)/Nm)
% % loglog(N,weight/Nm) % The plot
%  title('Mass spectrum at Nc=') % edit
%  xlabel('log10(N)')
%  ylabel('log10(weight/Nm)')

%%% FREQUENCY plots
% scatter(log10(N),log10(freq)) % The plot
% title('frequency spectrum at Nc=') % edit
% xlabel('log10(N)')
% ylabel('log10(freq)')

%############################################################%
%% PART 3 (Time dependency plots)
% num_clusta=zeros(numfiles,2);
% avg_mass=zeros(numfiles,1);
% for i=1:numfiles
%     [num_clusta(i,1),num_clusta(i,2)]=size(mass{i});
%     avg_mass(i)=mean(mass{i});
% end
% num_clust=num_clusta(:,1);
% x=time;
% y=1./num_clust;
% loglog(time,1./num_clust,'-s')
% title('Time dependence of 1/Nc at fv =')
% xlabel('Time')
% ylabel('Normalized number of clusters in the system(1/Nc)')
% grid on
% % scatter(time,avg_mass)
% % title('Time dependence of avg Mass') % edit
% % xlabel('Time')
% % ylabel('Avg mass')

%****************USE 'cftool' for curve fitting************%