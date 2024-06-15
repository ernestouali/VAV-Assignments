clear all
close all
clc

%% Load data and plot time histories
load('Data')

t=Data(:,1);
F=Data(:,2);
x=Data(:,3:end);
np=size(x,1);
nj=size(x,2);

dt=t(2)-t(1);
fsamp=1/dt;

disp(' ')
disp(['Number of sensors: ' num2str(nj)])
disp(['dt [s]: ' num2str(dt)])
disp(['np [-]: ' num2str(np)])

figure(1)
sf1(1)=subplot(3,1,1);
plot(t,F)
grid on
ylabel('F [N]')
title('Force')
sf1(2)=subplot(3,1,2);
plot(t,x(:,1),'b',t,x(:,2),'r')
grid on
ylabel('x_1,x_2 [m]')
title('Displacement')
sf1(3)=subplot(3,1,3);
plot(t,x(:,3),'b',t,x(:,4),'r')
grid on
xlabel('Time [s]')
ylabel('x_3,x_4 [m]')
title('Displacement')
linkaxes(sf1,'x')
sgtitle('Measurements')


%% Definition of Experimental TF

% Force and Measurement spectra
[xfft,frq]=ffg(x,np,dt); % We get the complex fft of the displacements
[Ffft,frq]=ffg(F,np,dt); % We get the complex fft of the applied force

% transfer functions
Hjkexp = xfft./Ffft; % Experimental transfer function is the ratio between input and output

figure(2) % This is the figure which plots the fft. Note the absence of nodes
sf2(1)=subplot(2,1,1);
for jj=1:nj
    plot(frq,abs(Hjkexp(:,jj)))
    hold on
    legenda{jj}=['H_' num2str(jj) '_k'];
end
grid on
title('Magnitude')
sgtitle('Experimental trasfer functions H_j_k')
ylabel(['|H_j_k| [m/N]'])
legend(legenda)
sf2(2)=subplot(2,1,2);
for jj=1:nj
    plot(frq,angle(Hjkexp(:,jj))*180/pi)
    hold on
end
grid on
ylabel('\angleH_j_k [deg]')
xlabel('Freq. [Hz]')
yticks([-180 -90 0 90 180])

title('Phase')
linkaxes(sf2,'x')
xlim([0 5])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parte da qui


%% Identification parameters single range
% disp(' ')
% nrange=input('Number of ranges to identify? ');

% for irange=1:nrange

    figure(2)
    disp(' ')
    disp('Define the range of frequencies')
    fini=input('Freq. min [Hz]: ');
    ffin=input('Freq. max [Hz]: ');    
    iini=min(find(round(frq*1000)/1000>=fini));
    ifin=max(find(round(frq*1000)/1000<=ffin));
    disp(' ')
    disp('Frequency range')
    disp(['- fmin [Hz]: ' num2str(frq(iini))])
    disp(['- fmax [Hz]: ' num2str(frq(ifin))])    
    npid=ifin-iini+1;
    disp(['Number of points for the identification: ' num2str(npid)])
    


    %%
    % Function reduced to the range of interest
    rfHjki=frq(iini:ifin);          % We extract the reduce frequency estimation range
    Hjkiexp=Hjkexp(iini:ifin,:);    % We extract the reduced frequency transfer function
    %%
    % First guess of parameters: SIMPLIFIED METHODS
    disp(' ')
    jj=input('Which diagram for identification? '); % Which measurement do you want to use to calculate the TF?
    
    %% Point 3  Modal parameters guessing
    % Natural frequency
    [vmax,iwmax]=max(abs(Hjkiexp(:,jj)));  % We identify the maximum value in our frequency range
    f0i=rfHjki(iwmax); %[Hz]                 We then compute the corresponding frequency value
    w0i0=2*pi*f0i; %[rad/s]                  Resonance frequency
    % Adimensional damping ration: we use the phase derivative
    derFIjki=(angle(Hjkiexp(iwmax+1,jj))-angle(Hjkiexp(iwmax-1,jj)))/(2*pi*(rfHjki(iwmax+1)-rfHjki(iwmax-1)));
    csii0=-1/(w0i0*derFIjki); % Guessed adimensional damping ratio
    r0i=2*w0i0*csii0;
    % Mode shapes
    Aj0=-imag(Hjkiexp(iwmax,jj))*w0i0*r0i; 
    % all the other constants are set equal to 0 (m_q = 1)

    disp(' ')
    disp(['Init guess f0 [Hz]: ' num2str(f0i)])
    disp(['Init guess csi [-]: ' num2str(csii0)])
    disp(['Init guess A0 [Hz]: ' num2str(Aj0)])
            
    % Filling of vector xpar
    xpar0=[csii0; w0i0; Aj0; zeros(5,1)]; % all the unknowns (2 mod parameters + 6 constants)
    % for the initial guess of the optimization
    
    % Identification: single channel    
    options=optimset('fminsearch');
    options=optimset(options,'TolFun',1e-8,'TolX',1e-8);
    xpar=fminsearch(@(xpar) errHjki_cw(xpar,rfHjki,Hjkiexp(:,jj)) ,xpar0,options);

    % @(xpar) parameter set we're trying to find
    % Error function (difference between experimental and optimized FRF)


    % Plot results of idetification
    vpar=[1; 2*xpar(1)*xpar(2); xpar(2)^2; xpar(3:8)];
    %    [m;   c = 2 m w0 csi; k = w0^2 m; A;B;C;D;E;F]

    Hjkiid=funHjki(vpar,rfHjki);
%%
    figure(3)
    sf3(1)=subplot(2,1,1);
    plot(rfHjki,abs(Hjkiexp(:,jj)),'*b-',rfHjki,abs(Hjkiid),'r--','linewidth',1.2)
    grid on
    xlabel('Frequency [Hz]')
    ylabel(['|H_' num2str(jj) '_k| [m/N]'])
    title('Magnitude')
    legend('Experimental','Identified')
    sf3(2)=subplot(2,1,2);
    plot(rfHjki,angle(Hjkiexp(:,jj)),'*b-',rfHjki,angle(Hjkiid),'r--','linewidth',1.2)
    grid on
    ylabel(['\angleH_' num2str(jj) '_k [rad]'])
    xlabel('Frequency [Hz]')
    title('Phase')
    legend('Experimental','Identified')
    linkaxes(sf3,'x')
    sgtitle(['Trasfer functions H_' num2str(jj) '_k comparison: Experimental vs Identified'])

%%
% end
