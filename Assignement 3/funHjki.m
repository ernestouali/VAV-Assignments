function Hjki=funHjki(vpar,rfHjkiexp)
% Reconstructing the experimental FRF from parameters vpar with the
% frequency range rfHjkiexp
% vpar -> vector of simplified identification formula parameters
if max(size(vpar))~=9
    disp('Parameter vector vpar: length error!')
    disp([' Length=' num2str(length(vpar)) ' instead of 9'])
    return    
end

mii=vpar(1);
cii=vpar(2);
kii=vpar(3);
Aj=vpar(4);
Bj=vpar(5);
Cj=vpar(6);
Dj=vpar(7);
Ej=vpar(8);
Fj=vpar(9);

% Computation of the simplified formula
Om=2*pi*rfHjkiexp;
Hjki=(Aj+1i*Bj)./(-mii*Om.^2+1i*Om*cii+kii)+(Cj+1i*Dj)+(Fj+1i*Ej)./(Om.^2);






