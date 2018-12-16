function out=OSA(sig,dt,resBandwidth,f_ref,name)
% sig=s_mod;
% dt=1/(param.f_ocs*param.oversampling);
% resBandwidth=1e8;
% f_ref=200e9;
% name='SB';

% out=OSA(sig,dt,resBandwidth,f_ref,name)

slen=length(sig);

% generate a frequency vector 
f_vec=0:slen/2-1;
% f_vec=0:slen-1;
f_vec=(f_vec/slen)/dt;

fft_sig = fft(sig ,slen);
spctrm=fft_sig;
fft_sig = abs(fft_sig(1:slen/2));
% fft_sig = abs(fft_sig(1:slen));
fft_sig=fft_sig.^2; % optical power

% resolution bandwidth
flen=max(1,round(resBandwidth*dt*slen));
% flen=max(1,round(500e9*dt*slen));
out.flen=flen;
% averaging rectangular filter
fft_sig=filter(ones(flen,1)/flen,1,circshift(fft_sig,-round(flen/2)));
% fft_sig=filter(ones(flen,1)/flen,1,fft_sig);

fft_sig=fft_sig/max(fft_sig); % normalize 
out.spctrm=spctrm;
out.freq=f_vec;
out.opt_power_density=fft_sig;

fig_handle=plot((f_vec-f_ref)/1e9,10*log10(fft_sig));
title(['Spectrum optical domain: ' name '; resBW= ' num2str(resBandwidth/1e6) ' MHz' ])
xlabel('Frequency (GHz)')
ylabel('|E(f)|^2  (dB)')
%axis([-f_ref/2e9 f_ref/2e9 -100 10])
axis([-400 400 -80 10])
grid on
set(fig_handle,'linewidth',2);
out.fig_handle=fig_handle;
end