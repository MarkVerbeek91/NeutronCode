%clear all

% plotting data from the NeutronCode
%data   = csvread('C:/Users/Mark/Dropbox/Projecten/NeutronCode/DATA.csv');

% r,phi[r],ParticleEnergy[r],f[r],g[0][r],A[r]
%data_c = csvread('C:/Users/Mark/Dropbox/Projecten/NeutronCode/DATA_cross_sections.csv');

% data from laptop
data = csvread('C:/Users/Mark/Documents/Development/NeutronCalculation/DATA.csv');
data_c = csvread('C:/Users/Mark/Documents/Development/NeutronCalculation/DATA_cross_sections.csv');


% r,data.phi[r],data.ParticleEnergy[r],data.f[r],data.g[0][r],data.A[r]

%% %potential
figure 
potential = plot(data(:,1),data(:,2))
legend(potential, 'Phi')

%% 
figure 
energies = plot(data(:,1),data(:,3))
legend(energies, 'E')
% figure
% surfival = plot(data_c(:,1),data_c(:,2))
% legend(surfival, 'SIIEE')
figure
surfival = plot(data(:,1),data(:,4),data(:,1),data(:,5))
legend(surfival, 'f(r)', 'g(r,r1)')
ylim([0, 1]);
%%
% figure
% cross = loglog(data_c(:,1),data_c(:,2),data_c(:,1),data_c(:,3),data_c(:,1),data_c(:,4))
% legend(cross, 'CX','Ion','Tot')

%%
% figure
% cross = semilogy(data_c(:,1),data_c(:,5))
% axis([50000 500000 4.708856E-031 6.272065E-030])
% legend(cross, 'Fusion')



% 
% , ...
%      data(:,3),data(:,4), ...
%      data(:,5),data(:,6), ...
%      data(:,7),data(:,8), ...
%      data(:,9),data(:,10), ...
%      data(:,11),data(:,12), ...
%      data(:,13),data(:,14), ...
%      data(:,15),data(:,16), ...
%      data(:,17),data(:,18), ...
%      data(:,19),data(:,20), ...
%      data(:,21),data(:,22), ...
%      data(:,23),data(:,24), ...
%      data(:,25),data(:,26), ...
%      data(:,27),data(:,28), ...
%      data(:,29),data(:,30))
