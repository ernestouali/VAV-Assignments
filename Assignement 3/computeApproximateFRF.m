function [h, res_freq, mode_shape, best_param] = computeApproximateFRF(fini,ffin, jj, frq_axis, Hjk_exp)
% INPUT :   fini Lower f in band
%           ffin Higher f in band
%           jj FRF identifier
%           frq_axis : frequency axis
%           Hjk_exp : experimental FRF computed from the data file
% OUTPUT :  h estimated adimensional damping ratio  [-]
%           res_freq estimated resonance frequency [Hz]
%           estimated mode shape [m]
%           best_param dim(2+6*nm,1)  estimated best choice of parameters for the
%           approximated FRF

figure()
sgtitle("Residual minimization: FRF_" + num2str(jj));
h = zeros(4,1);
res_freq = zeros(4,1);
mode_shape = zeros(4,1);
best_param = zeros(9,4);
r0i = zeros(4,1);
w0i = zeros(4,1);


for i = 1:4        % for each mode
    %==========================
    %           Q2
    %==========================
    
    % Indexes indentification
    iini = min(find(round(frq_axis*1000)/1000 >= fini(i)));
    ifin = max(find(round(frq_axis*1000)/1000 <= ffin(i)));
    
    % Count of points for identification
    npid = ifin - iini + 1;
    
    % FRF evaluated in [fini, ffin]
    sub_freq_axis = frq_axis(iini:ifin);
    Hjkiexp = Hjk_exp(iini:ifin,:);
    
    % Natural frequency
    [vmax, iwmax] = max(abs(Hjkiexp(:, jj)));        % Maximum of |FRF| inside [fini, ffin]
    res_freq(i) = sub_freq_axis(iwmax);         % Frequency [Hz] of the maximum
    w0i(i) = 2*pi*res_freq(i);                  % Resonance frequency [rad/s]
    
    % Adimensional damping ration (phase derivative)
    derFIjki = (angle(Hjkiexp(iwmax+1, jj)) - angle(Hjkiexp(iwmax-1, jj)))/(2*pi*(sub_freq_axis(iwmax+1) - sub_freq_axis(iwmax-1)));
    h(i) = -1/(w0i(i)*derFIjki);        % Guessing h
    r0i(i) = 2 * w0i(i) * h(i);        % c/m
    
    % Mode shapes
    mode_shape(i) = -imag(Hjkiexp(iwmax, jj)) * w0i(i) * r0i(i); 
    
    
    %==========================
    %           Q3
    %==========================

    % Filling of vector xpar
    xpar0 = [h(i); w0i(i); mode_shape(i); zeros(5,1)];    % Unknowns: 2 mod parameters + 6 constants
    
    % Identification: single channel    
    options = optimset('fminsearch');
    options = optimset(options,'TolFun',1e-8,'TolX',1e-8);
    xpar = fminsearch(@(xpar) errHjki_cw(xpar, sub_freq_axis, Hjkiexp(:, jj)), xpar0, options);
    
    % @(xpar) parameter set we're trying to find
    % Error function (difference between experimental and optimized FRF)
    
    % Plot results of identification
    best_param(:,i) = [1;    2*xpar(1)*xpar(2);    xpar(2)^2;    xpar(3:8)];
    %               [m;    c = 2 m w0 csi;       k = w0^2 m;  A;B;C;D;E;F]
    
    %Reconstructing the experimental FRF after solving the minimization
    Hjkiid = funHjki(best_param(:,i), sub_freq_axis);
    
    % Plots
    subplot(2, 4, i);
    plot(sub_freq_axis, abs(Hjkiexp(:, jj)), '*b-', sub_freq_axis, abs(Hjkiid), 'r--', 'linewidth', 1.2)
    grid on;
    title("Magnitude of mode " + num2str(i));
    ylabel("[m/N]");
    xlabel("[Hz]");
    legend("Experimental", "Identified");
    
    subplot(2, 4, i+4);
    plot(sub_freq_axis, angle(Hjkiexp(:, jj)), '*b-', sub_freq_axis, angle(Hjkiid), 'r--', 'linewidth', 1.2)
    grid on;
    title("Phase of mode " + num2str(i));
    ylabel("[°]");
    xlabel("[Hz]");
    yticks([-pi, -pi/2, 0, pi/2, pi]);
    yticklabels(["-\pi", "-\pi/2", "0", "\pi/2", "\pi"]);
    legend("Experimental", "Identified");
end

end