%Empty workspace and close figures
close all;
clear;


%Define the SNR
SNR = 1;

%Define strInterFerence (strength of inter-cell interference)
strInterFerence = 1e-1;

%Define range of number of BS antennas
M = 1:100;

%Extract the maximal number of BS antennas
Mmax = max(M);

%Select number of Monte Carlo realizations of the Rayleigh fading 
numberOfRealizations = 100000;

%Generate random UE angles from 0 to 2*pi
varphiDesired = 2*pi*rand(1,numberOfRealizations);
varphiInterfering = 2*pi*rand(1,numberOfRealizations);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance


%Preallocate matrices for storing the simulation results
SE_LoS = zeros(Mmax,numberOfRealizations);
SE_NLoS_montecarlo = zeros(Mmax,numberOfRealizations);


%Generate uncorrelated Rayleigh fading realizations
fadingRealizationsDesired = (randn(Mmax,numberOfRealizations)+1i*randn(Mmax,numberOfRealizations))/sqrt(2);
fadingRealizationsInterfering = (randn(Mmax,numberOfRealizations)+1i*randn(Mmax,numberOfRealizations))/sqrt(2);


%Go through different channel realizations
for k = 1:numberOfRealizations
    
    %Compute the argument, that appears in (Equation 4.7), for different UE angles
    argument = (2*pi*antennaSpacing*( sin(varphiDesired(k))  - sin(varphiInterfering(k))) );
    
    %Compute the SE under line-of-sight (LoS) propagation, for a particular
    %set of UE angles, using (Equation 4.6). Note that the g-function is
    %implemented slightly differently from (Equation 4.7)
    SE_LoS(:,k) = log2(1 + SNR*M  ./ (  (SNR*strInterFerence/(1-cos(argument)))*(1-cos(argument*M)) ./ M + 1 ) );
    
    %Compute the SE under NLoS propagation for one fading realization
    SE_NLoS_montecarlo(:,k) = log2(1 + SNR*cumsum(abs(fadingRealizationsDesired(:,k)).^2) ./ (  strInterFerence*SNR*(abs(cumsum(conj(fadingRealizationsDesired(:,k)).*fadingRealizationsInterfering(:,k))).^2)./cumsum(abs(fadingRealizationsDesired(:,k)).^2) +1 ) );
    
end

%Compute the lower bound on SE under NLoS propagation in (Equation 4.12)
SE_NLoS_lower = log2(1 + SNR*(M-1)  ./ ( SNR*strInterFerence + 1 ) );



%Compute the SE under NLoS propagation using (Equation 4.8)
SE_NLoS =  (exp(1/(SNR*strInterFerence))*expint(1/(SNR*strInterFerence))/log(2))*(1./(1-1/strInterFerence).^M -1 );

%First summation in (Equation 4.8) by considering all the M values at once
for m = 1:Mmax
    
    disp([num2str(m) ' antenna numbers out of ' num2str(Mmax)]);
    
    %First term that depend only on m
    term1 = 1/strInterFerence/(1-1/strInterFerence)^m;
    
    %Second summation
    for l = 0:Mmax-m
        
        %Second term that depend only on m and l
        term2 = term1 * (-1).^(M-m-l+1) ./ (SNR.^(M-m-l)  .* factorial(abs(M-m-l)) * log(2));
        
        %Third and fourth summation
        term3 = exp(1/SNR)*expint(1/SNR);
        
        for n = 1:l
            
            for j = 0:n-1
                
                term3 = term3 + 1/n/factorial(j)/SNR^j;
                
            end
            
        end
        
        %Determine which of the considered M values that contain the current term
        feasible = (m<=M) & (l<=(M-m));
        
        %Add the terms to the SE computation
        SE_NLoS(feasible) = SE_NLoS(feasible) + term2(feasible)*term3;
        
        
    end
    
end


%% Plot the simulation results
figure;
hold on; box on;

plot(M,mean(SE_LoS,2),'k-','LineWidth',1);
plot(M(1),mean(SE_NLoS_montecarlo(1,:),2),'bd-.','LineWidth',1);
plot(M,SE_NLoS_lower,'r--','LineWidth',1);
plot(M(10:10:end),mean(SE_NLoS_montecarlo(10:10:end,:),2),'bd','LineWidth',1);
plot(M,SE_NLoS,'b-.','LineWidth',1);

xlabel('Number of antennas (M)');
ylabel('Average SE [bit/s/Hz]');

legend('LoS','NLoS: Rayleigh fading','NLoS: Rayleigh fading (lower bound)','Location','SouthEast');
