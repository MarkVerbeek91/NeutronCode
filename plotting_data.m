% plotting data from the NeutronCode
str = 'DATA.csv';
data = csvread(str);

% r,data.phi[r],data.ParticleEnergy[r],data.f[r],data.g[0][r],data.A[r]
figure %potential
plot(data(:,1),data(:,2))
figure %particle energies
plot(data(:,1),data(:,3))
figure
plot(data(:,3),data(:,4),data(:,3),data(:,5))
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
