function E_out=SSMF(E_in,L,D,dt,f_opt,Pin,alpha,gamma)
% SSMF(s_mod,param.fiber_length,param.dispersion,dt,param.f_opt,param.Pin,param.alfa,param.gamma);
% L: fiber length in m
% D: dispersion in ps/nm/km; e.g. SMF @ 1550: 17 ps/nm/km 
% wavelength lambda is set to 1550 nm, important for calculation of beta2
% out of D
% dt: time between succeeding sampling points

E_in=E_in(:);
% Fiber transfer function
% H_fiber(f,z) = exp(-i*beta(f,f0)*z), 
% beta(f,f0)=n(f,f0)*2*pi*f0/c_light-beta(0,f0); (move with v_group at f0)
% calculation of refraction index
% based on Sellmeier equation

c0=3e8;
lambda=c0/f_opt; % wavelength in m
z=L;

% rough approximation of best step size
dz=L;
if Pin>1e-4
    dz=D*1e7/17/Pin;
end
% dz=1000;

% Calculation of dispersion effect in frequency domain via fft

f_max=1/(dt*2);
slen=length(E_in);
% generate a frequency vector 
f_vec=-slen/2:slen/2-1;
f_vec=(f_vec/slen)/dt;
f_vec=fftshift(f_vec)';
beta2=-D*lambda^2/2/pi/c0; 
omega=2*pi*(f_vec);

H=exp(-1i*1/2*beta2*(omega).^2*dz);
H(slen/2:end)= 0;  % single sided!!!  

% dispersion and nonlinearity: split step fourier = fragment fiber in 
% short sections << effective fiber length and applying operator 
% for dispersion and operator for nonlinearity successively for each
% section
% nonlinearity: time domain, dispersion: frequency domain
% dispersion
for k=1:round(z/dz) % split step fourier
    disp(strcat(num2str(k),' / ', num2str(z/dz)));
    % nonlinearity (total field approach "SPM" includes: SPM, XPM, FWM)
    E_out=SPM(E_in,dz,alpha,gamma,Pin);
    U_in=fft(E_out);
    U_out=U_in.*H;
    E_out=ifft(U_out); % ifft to time domain
    E_out=E_out(:);
    E_in=E_out;
    Pout=Pin*exp(-alpha*dz);
    Pin=Pout;
end

end



function E_out=SPM(E_in,z,alpha_SPM,gamma,P)
% normalizing E_in to power = 1
E_in=E_in/sqrt(mean(abs(E_in).^2));
z_eff=1/alpha_SPM*(1-exp(-alpha_SPM*z/2));
PHI_nl=P*gamma*z_eff*abs(E_in).^2;
E_out=E_in.*exp(1i*PHI_nl);
E_out=E_out(:);
end