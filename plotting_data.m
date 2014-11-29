%clear all

loc = 'C:/Users/Mark/Documents/Development/NeutronCalculation/';

%% plotting Potential
% filename = strcat(loc, 'Potential.csv');
% data_Potential = csvread(filename);
% figure 
% Phi = plot(data_Potential(:,1),data_Potential(:,2));
% legend(Phi, 'Phi')

%% plotting SIIEE
% filename = strcat(loc, 'SIIEE.csv');
% data_SIIEE = csvread(filename);
% figure 
% SIIEE = plot(data_SIIEE(:,1),data_SIIEE(:,2));
% legend(SIIEE, 'SIIEE')

%% plotting Cross section of Charge Exchange
% filename = strcat(loc, 'CrosssecCX.csv');
% data_CrosssecCX = csvread(filename);
% figure 
% CrosssecCX = loglog(data_CrosssecCX(:,1),data_CrosssecCX(:,2));
% legend(CrosssecCX, 'CrosssecCX')

% filename1 = strcat(loc, 'CrosssecCX1.csv');
% data_CrosssecCX1 = csvread(filename1);
% CrosssecCX = loglog(data_CrosssecCX1(:,1),data_CrosssecCX(:,2)-data_CrosssecCX1(:,2));
% legend(CrosssecCX, 'CrosssecCX')

% 
% %% plotting Cross section of Iononisation
% filename = strcat(loc, 'CrosssecIon.csv');
% data_CrosssecIon = csvread(filename);
% % figure 
% % CrosssecIon = loglog(data_CrosssecIon(:,1),data_CrosssecIon(:,2));
% % legend(CrosssecIon, 'CrosssecIon')
% 
% %% plotting Cross section of Iononisation
% filename = strcat(loc, 'CrosssecTot.csv');
% data_CrosssecTot = csvread(filename);
% % figure 
% % CrosssecTot = loglog(data_CrosssecTot(:,1),data_CrosssecTot(:,2));
% % legend(CrosssecTot, 'CrosssecTot')
% 
% %% plotting all Cross section in one. 
% Crosssecs = loglog(data_CrosssecCX(:,1),data_CrosssecCX(:,2), data_CrosssecIon(:,1),data_CrosssecIon(:,2),data_CrosssecTot(:,1),data_CrosssecTot(:,2));
% legend(Crosssecs, 'Cross section CX','Cross section Ion','Cross section Tot')

%% plotting the survival funtion f
% filename = strcat(loc, 'f.csv');
% data_f = csvread(filename);
% figure 
% f = plot(data_f(:,1),data_f(:,2));
% legend(f, 'f')
% xlim([0 0.25]);
% ylim([0 1]);


%% plotting the survival function g(0,r)
filename = strcat(loc, 'g.csv');
data_g = csvread(filename);
filename = strcat(loc, 'Mathematica\g1.csv');
data_g1 = csvread(filename);

data_g3 = data_g(:,2) - data_g1(:,2);

figure 
g = plot(data_g(:,1),data_g(:,2),data_g1(:,1),data_g1(:,2),data_g(:,1),data_g3);
%plot(data_g(:,1),data_g(:,2))
legend(g, 'G code - G mathematica')
xlim([0 0.25]);
ylim([0 1]);

%% plotting A
% filename = strcat(loc, 'A.csv');
% data_A = csvread(filename);
% figure 
% A = plot(data_A(:,1),data_A(:,2));
% legend(A, 'A')
% xlim([0.05 0.25]);

%% verivacation of the Kernel
% filename = strcat(loc, 'k.csv');
% data_K = csvread(filename);
% filename = strcat(loc, 'k1.csv');
% data_K1 = csvread(filename);
% 
% data_K(:,2) = data_K(:,2)-data_K1(:,2);
% 
% figure 
% K = plot(data_K(:,1),data_K(:,2));
% legend(K, 'Kernel')
% 

%% plotting the Source function
% filename = strcat(loc,'S4.csv');
% data = csvread(filename);
% plot(data(:,1) , data(:,2),'.')

%% plotting spectrums of the ions

filename = strcat(loc,'spectrum.csv');
data = csvread(filename);
semilogy(data(:,1),abs(data(:,2)))
%plot(data(:,1),data(:,2))


