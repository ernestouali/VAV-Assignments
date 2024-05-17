 function [num_cmplx,freq]=ffg(Dati_NMV,N,dt)
 grf=fft(Dati_NMV,N);
 mod=[grf(1,:)/N; 2/N*abs(grf(2:floor(N/2),:))];
 freq=[(1/dt)/N*(0:(N/2-1))]';
 fase=atan2(imag(grf(1:floor(N/2),:)),real(grf(1:floor(N/2),:)));
 num_cmplx=mod.*exp(1i*fase);

