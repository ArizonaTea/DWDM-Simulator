function out=ESA(sig,dt,resBandwidth,name)

% out=ESA(sig,dt,resolution,name)

slen=length(sig);

% generate a frequency vector 
f_vec=0:slen/2-1;
f_vec=(f_vec/slen)/dt;

fft_sig = fft(sig);
spctrm=fft_sig;
fft_sig = abs(fft_sig(1:slen/2));

% resolution bandwidth
flen=max(1,round(resBandwidth*dt*slen));
out.flen=flen;
% averaging rectangular filter
fft_sig=filter(ones(flen,1)/flen,1,circshift(fft_sig,-round(flen/2)));

fft_sig=fft_sig/max(fft_sig); % normalize 
out.spctrm=spctrm;
out.freq=f_vec;
out.abs_fft_sig=fft_sig;

fig_handle=plot(f_vec/1e9,20*log10(fft_sig));
title(['Spectrum electrical domain: ' name '; resBW= ' num2str(resBandwidth/1e6) ' MHz' ])
xlabel('Frequency (GHz)')
ylabel('|E(f)| (dB)')
% the axis setting is somewhat arbitrary; adapt for your needs
axis([0 1/(dt*4)/1e9 -80 10])
grid on
set(fig_handle,'linewidth',2);
out.fig_handle=fig_handle;
end