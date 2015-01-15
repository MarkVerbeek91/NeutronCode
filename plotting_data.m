%clear all

%LoadData

showPotential       = false;
showSIIEE           = false;
showCrossSections   = false;
showSurvial_f       = false;
showSurvial_g       = false;
showA               = false;
showK               = false;
showSourceFunc      = false;
showSpectrum        = true;

loc = pwd;

%% plotting Potential
if showPotential
    FileName = strcat(loc, '\Potential.csv');
    if exist(FileName, 'file')
        data_Potential = csvread(FileName);
    end
    
    figure 
    Phi = plot(data_Potential(:,1),data_Potential(:,2));
    legend(Phi, 'Phi')
end
    
%% plotting SIIEE
if showSIIEE
    FileName = strcat(loc, '\SIIEE.csv');
    if exist(FileName, 'file')
        data_SIIEE = csvread(filename);
    end

    figure 
    SIIEE = plot(data_SIIEE(:,1),data_SIIEE(:,2));
    legend(SIIEE, 'SIIEE')
end

%% plotting Cross section of Charge Exchange
if showCrossSections
    FileName = strcat(loc, '\CrosssecCX.csv');
    if exist(FileName, 'file')
        data_Potential = csvread(FileName);
    end
    
    FileName = strcat(loc, '\CrosssecIon.csv');
    if exist(FileName, 'file')
        data_SIIEE = csvread(filename);
    end

    FileName = strcat(loc, '\CrosssecTot.csv');
    if exist(FileName, 'file')
        data_Potential = csvread(FileName);
    end
        
    figure 
    Crosssecs = loglog(data_CrosssecCX(:,1),data_CrosssecCX(:,2), data_CrosssecIon(:,1),data_CrosssecIon(:,2),data_CrosssecTot(:,1),data_CrosssecTot(:,2));
    legend(Crosssecs, 'Cross section CX','Cross section Ion','Cross section Tot')
end

%% plotting the survival funtion f
if showSurvial_f
    FileName = strcat(loc, '\f.csv.csv');
    if exist(FileName, 'file')
        data_SIIEE = csvread(filename);
    end

    figure 
    f = plot(data_f(:,1),data_f(:,2));
    legend(f, 'f')
    xlim([0 0.25]);
    ylim([0 1]);
end

%% plotting the survival function g(0,r)
if showSurvial_g
    FileName = strcat(loc, '\g.csv');
    if exist(FileName, 'file')
        data_Potential = csvread(FileName);
    end

    figure 
    plot(data_g(:,1),data_g(:,2))
    
    xlim([0 0.25]);
    ylim([0 1]);
end

%% plotting A
if showA
    FileName = strcat(loc, '\A.csv');
    if exist(FileName, 'file')
        data_SIIEE = csvread(filename);
    end
    
    figure
    A = plot(data_A(:,1),data_A(:,2));
    legend(A, 'A')
    xlim([0.05 0.25]);
end
%% verivacation of the Kernel
if showK
    FileName = strcat(loc, '\k.csv');
    if exist(FileName, 'file')
        data_SIIEE = csvread(filename);
    end

    figure 
    K = plot(data_K(:,1),data_K(:,2));
    legend(K, 'Kernel')
end

%% plotting the Source function
if showSourceFunc
    FileName = strcat(loc,'S4.csv');
    data = csvread(FileName);
    plot(data(:,1) , data(:,2),'.')
end

%% plotting spectrum surface
if showSpectrum

    data_total1 = zeros(3999,24);
    
    radii = 1:1:24;
    
    figure
    hold on
    for i=1:24
        FileName = strcat(loc, '\IonSpectrumInwards', int2str(i),'.csv');
        if exist(FileName, 'file')
            data1 = csvread(FileName);
        end
        
        FileName = strcat(loc, '\IonSpectrumOutwards', int2str(i),'.csv');
        if exist(FileName, 'file')
            data2 = csvread(FileName);
        end
        
        semilogy(data1(:,1),data1(:,2))
        data_total1(:,i) = data1(:,2);    
    end

    hold off
    set(gca, 'YScale', 'log')
    %%
    %semilogy(data(:,1),data1)

    %surf(data1)
    %set(gca, 'ZScale', 'log')
end

%%
surf(radii,data1(:,2)*10000,data_total1)
set(gca, 'ZScale', 'log')
