%clear all

%LoadData

showPotential   = false;
showSIIEE       = false;
showCrossSections   = false;
showSurvial_f   = false;
showSurvial_g   = false;
showA           = true;
showK           = false;
showSourceFunc  = false;
showSpectrum    = false;


%loc = 'C:/Users/Mark/Documents/Development/NeutronCalculation/';
loc = 'C:\Users\Mark\SkyDrive\Documenten\MATLAB\NeutronCode\new\';

filename = strcat(loc, 'Potential.csv');
data_Potential = csvread(filename);

filename = strcat(loc, 'SIIEE.csv');
data_SIIEE = csvread(filename);

filename = strcat(loc, 'CrosssecCX.csv');
data_CrosssecCX = csvread(filename);

filename = strcat(loc, 'CrosssecIon.csv');
data_CrosssecIon = csvread(filename);

filename = strcat(loc, 'CrosssecTot.csv');
data_CrosssecTot = csvread(filename);

filename = strcat(loc, 'f.csv');
data_f = csvread(filename);

filename = strcat(loc, 'g.csv');
data_g = csvread(filename);

filename = strcat(loc, 'A.csv');
data_A = csvread(filename);

filename = strcat(loc, 'k.csv');
data_K = csvread(filename);

%% plotting Potential
if showPotential
    filename = strcat(loc, 'Potential.csv');
    data_Potential = csvread(filename);
    figure 
    Phi = plot(data_Potential(:,1),data_Potential(:,2));
    legend(Phi, 'Phi')
end
    
%% plotting SIIEE
if showSIIEE
    filename = strcat(loc, 'SIIEE.csv');
    data_SIIEE = csvread(filename);
    figure 
    SIIEE = plot(data_SIIEE(:,1),data_SIIEE(:,2));
    legend(SIIEE, 'SIIEE')
end

%% plotting Cross section of Charge Exchange
if showCrossSections
    filename = strcat(loc, 'CrosssecCX.csv');
    data_CrosssecCX = csvread(filename);
    figure 
    CrosssecCX = loglog(data_CrosssecCX(:,1),data_CrosssecCX(:,2));
    legend(CrosssecCX, 'CrosssecCX')

    filename1 = strcat(loc, 'CrosssecCX1.csv');
    data_CrosssecCX1 = csvread(filename1);
    CrosssecCX = loglog(data_CrosssecCX1(:,1),data_CrosssecCX(:,2)-data_CrosssecCX1(:,2));
    legend(CrosssecCX, 'CrosssecCX')


    %% plotting Cross section of Iononisation
    filename = strcat(loc, 'CrosssecIon.csv');
    data_CrosssecIon = csvread(filename);
    figure 
    CrosssecIon = loglog(data_CrosssecIon(:,1),data_CrosssecIon(:,2));
    legend(CrosssecIon, 'CrosssecIon')

    %% plotting Cross section of Iononisation
    filename = strcat(loc, 'CrosssecTot.csv');
    data_CrosssecTot = csvread(filename);
    figure 
    CrosssecTot = loglog(data_CrosssecTot(:,1),data_CrosssecTot(:,2));
    legend(CrosssecTot, 'CrosssecTot')

    %% plotting all Cross section in one. 
    Crosssecs = loglog(data_CrosssecCX(:,1),data_CrosssecCX(:,2), data_CrosssecIon(:,1),data_CrosssecIon(:,2),data_CrosssecTot(:,1),data_CrosssecTot(:,2));
    legend(Crosssecs, 'Cross section CX','Cross section Ion','Cross section Tot')

end

%% plotting the survival funtion f
if showSurvial_f
    filename = strcat(loc, 'f.csv');
    data_f = csvread(filename);
    figure 
    f = plot(data_f(:,1),data_f(:,2));
    legend(f, 'f')
    xlim([0 0.25]);
    ylim([0 1]);
end

%% plotting the survival function g(0,r)
if showSurvial_g
    %filename = strcat(loc, 'Mathematica\g1.csv');
    %data_g1 = csvread(filename);

    %data_g3 = data_g(:,2) - data_g1(:,2);

    figure 
    % g = plot(data_g(:,1),data_g(:,2),data_g1(:,1),data_g1(:,2),data_g(:,1),data_g3);
    % legend(g, 'G code - G mathematica')
    plot(data_g(:,1),data_g(:,2))
    
    xlim([0 0.25]);
    ylim([0 1]);
end

%% plotting A
if showA
    filename = strcat(loc, 'A.csv');
    data_A = csvread(filename);
    figure 
    A = plot(data_A(:,1),data_A(:,2));
    legend(A, 'A')
    xlim([0.05 0.25]);
end
%% verivacation of the Kernel
if showK
    filename = strcat(loc, 'k.csv');
    data_K = csvread(filename);
    filename = strcat(loc, 'k1.csv');
    data_K1 = csvread(filename);

    data_K(:,2) = data_K(:,2)-data_K1(:,2);

    figure 
    K = plot(data_K(:,1),data_K(:,2));
    legend(K, 'Kernel')
end

%% plotting the Source function
if showSourceFunc
    filename = strcat(loc,'S4.csv');
    data = csvread(filename);
    plot(data(:,1) , data(:,2),'.')
end
    
%% plotting spectrums of the ions
if showSpectrum
    filename = strcat(loc,'spectrum1.csv');
    data = csvread(filename);
    filename = strcat(loc,'spectrum2.csv');
    data1 = csvread(filename);
    semilogy(data(:,1),abs(data(:,2)),data(:,1),data1(:,2))
    plot(data(:,1),data(:,2))
end

%% plotting spectrum surface
if false

    data1 = zeros(3999,24);

    for i=1:19
        filename = strcat(loc,'f_min_',int2str(i),'.csv')
        data = csvread(filename);
        data1(:,i) = data(:,2);    
    end

    %%
    semilogy(data1)

    %surf(data1)
    %set(gca, 'ZScale', 'log')
end