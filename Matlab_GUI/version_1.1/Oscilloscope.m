function out = Oscilloscope(sig,dt,name,ya,yb)

% out=Oscilloscope(signal,dt,name)
% Note that this oscilloscope will assume the perfect sampling rate. That
% is, the scope can catch every sample point

slen=length(sig);
xaxis=0:dt:(slen-1)*dt;
xaxis=xaxis';

fig_handle=plot(xaxis,sig);
title(['Oscilloscope: ' name ]);
xlabel('Time (s)')
ylabel('E(t) mW')
axis([0 slen*dt ya yb])
grid on
set(fig_handle,'linewidth',2);
out.fig_handle=fig_handle;

end

