function out=Constellation(sig,dt,BaudRate,f_opt,show_transition,name)

% out=Constellation(sig,dt,BaudRate,f_opt,show_transisiton,name)

slen=length(sig);
t=(1:slen)'*dt;

% generate a frequency vector
f_vec=-slen/2:slen/2-1;
f_vec=(f_vec/slen)/dt;
f_vec=fftshift(f_vec);
% down convert
dwn=sig.*exp(1j*2*pi*f_opt*t);

%filter in order to remove high frequency part
H_Gauss=exp(-log(2)/2*(f_vec/BaudRate).^2)';
down_flt=ifft(fft(dwn).*H_Gauss);

I_vec=real(down_flt);
Q_vec=imag(down_flt);
if show_transition==1
    I_vec_sample=I_vec;
    Q_vec_sample=Q_vec;
else
    hmax=0;
    samples_symbol=1/(BaudRate*dt);
    for i=1:round(samples_symbol)
        samplevec=round(i:samples_symbol:slen);
        hist_IQ=hist(atan2(I_vec(samplevec),Q_vec(samplevec)),100);
        hmax_new=max(hist_IQ);
        if(hmax_new>hmax)
            hmax=hmax_new;
            I_vec_sample=I_vec(samplevec);
            Q_vec_sample=Q_vec(samplevec);
        end
    end
end



fig_handle=plot(I_vec_sample,Q_vec_sample,'.');
title(['Constellation diagram: ' name])
xlabel('I')
ylabel('Q')
axis([min(I_vec_sample)*1.1 max(I_vec_sample)*1.1 min(Q_vec_sample)*1.1 max(Q_vec_sample)*1.1])
grid on
out.fig_handle=fig_handle;
end