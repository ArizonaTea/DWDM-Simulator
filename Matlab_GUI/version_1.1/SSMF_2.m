function E_out=SSMF_2(E_in,L,D,dt,f_opt,Pin,alpha,gamma)
% SSMF(s_mod,param.fiber_length,param.dispersion,dt,param.f_opt,param.Pin,param.alfa,param.gamma);
% L: fiber length in m
% D: dispersion in ps/nm/km; e.g. SMF @ 1550: 17 ps/nm/km 
% wavelength lambda is set to 1550 nm, important for calculation of beta2
% out of D
% dt: time between succeeding sampling points

sig=E_in;

C_light=299792458; % speed of light in meters %wavelength of light in meters
Wavelength=C_light/f_opt; %wavelength of light in meters
% Sampling_points=param.oversampling;
dispersion_slope=0.09e3;
launch_power=Pin;


length_fiber_m=L;
length_section_m=1000;
NN=(length_fiber_m/length_section_m);

slen=length(E_in);
% generate a frequency vector 
f_vec=-slen/2:slen/2-1;
f_vec=(f_vec/slen)/dt;
f_vec=fftshift(f_vec)';

% alpha=param.alpha; %Ref page#55 eqn 2.5.3 Fiber optic Comm by GP Agrawal
% gamma=param.gamma;
Beta2=-1*D*(Wavelength^2)/(2*pi*C_light);
Beta3=(dispersion_slope-(4*pi*C_light*Beta2)/(Wavelength^3))/(((2*pi*C_light)/(Wavelength^2))^2);
Beta3=0;

omega=2*pi*f_vec;

sig=sig/sqrt(mean(abs(sig).^2)); %normalize sig to power=1
sig=sig*sqrt(launch_power);


for jj=1:NN
    % linera part 1
    spectrum=fft(sig);
    Linear_operator_fre_dom=exp(1i*((Beta2)/2*omega.^2+(Beta3/6)*omega.^3)*(length_section_m/2))*exp(-(alpha/2*(length_section_m/2)));
    spectrum=spectrum.*Linear_operator_fre_dom;
    sig=ifft(spectrum);
    % non-linear part
    Nonlinear_operator=exp(1i*gamma*(abs(sig).^2)*(length_section_m));
    sig=sig.*Nonlinear_operator;
%     a=launch_power*gamma*(mean(abs(sig).^2))*(length_section_m)
    % linear part 2
    spectrum=(fft(sig));
    spectrum=spectrum.*Linear_operator_fre_dom;
    sig=ifft(spectrum);
end

E_out=sig;

end


% E_in=E_in(:);
% % Fiber transfer function
% % H_fiber(f,z) = exp(-i*beta(f,f0)*z), 
% % beta(f,f0)=n(f,f0)*2*pi*f0/c_light-beta(0,f0); (move with v_group at f0)
% % calculation of refraction index
% % based on Sellmeier equation
% 
% c0=3e8;
% lambda=c0/f_opt; % wavelength in m
% z=L;
% 
% % rough approximation of best step size
% dz=L;
% if Pin>1e-4
%     dz=D*1e7/17/Pin;
% end
% % dz=1000;
% 
% % Calculation of dispersion effect in frequency domain via fft
% 
% f_max=1/(dt*2);
% slen=length(E_in);
% % generate a frequency vector 
% f_vec=-slen/2:slen/2-1;
% f_vec=(f_vec/slen)/dt;
% f_vec=fftshift(f_vec)';
% beta2=-D*lambda^2/2/pi/c0; 
% omega=2*pi*(f_vec);
% 
% H=exp(-1i*1/2*beta2*(omega).^2*dz);
% H(slen/2:end)= 0;  % single sided!!!  
% 
% % dispersion and nonlinearity: split step fourier = fragment fiber in 
% % short sections << effective fiber length and applying operator 
% % for dispersion and operator for nonlinearity successively for each
% % section
% % nonlinearity: time domain, dispersion: frequency domain
% % dispersion
% for k=1:round(z/dz) % split step fourier
%     disp(strcat(num2str(k),' / ', num2str(z/dz)));
%     % nonlinearity (total field approach "SPM" includes: SPM, XPM, FWM)
%     E_out=SPM(E_in,dz,alpha,gamma,Pin);
%     U_in=fft(E_out);
%     U_out=U_in.*H;
%     E_out=ifft(U_out); % ifft to time domain
%     E_out=E_out(:);
%     E_in=E_out;
%     Pout=Pin*exp(-alpha*dz);
%     Pin=Pout;
% end
% 
% end
% 
% 
% 
% function E_out=SPM(E_in,z,alpha_SPM,gamma,P)
% % normalizing E_in to power = 1
% E_in=E_in/sqrt(mean(abs(E_in).^2));
% z_eff=1/alpha_SPM*(1-exp(-alpha_SPM*z/2));
% PHI_nl=P*gamma*z_eff*abs(E_in).^2;
% E_out=E_in.*exp(1i*PHI_nl);
% E_out=E_out(:);
% end