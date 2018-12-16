% more_carriers

%this code just adds two more carriers in order to demonstrate the
%nonlinear fiber properties.

spacing=param.spacing;
decorrelate=samples_symbol*20;

% Carriers generation
s_optl_1=sin(2*pi*(param.f_ocs-spacing)*t);
s_optl_2=sin(2*pi*(param.f_ocs-spacing*2)*t);
s_optl_3=sin(2*pi*(param.f_ocs-spacing*3)*t);
s_optl_4=sin(2*pi*(param.f_ocs-spacing*4)*t);

s_optu_1=sin(2*pi*(param.f_ocs+spacing)*t);
s_optu_2=sin(2*pi*(param.f_ocs+spacing*2)*t);
s_optu_3=sin(2*pi*(param.f_ocs+spacing*3)*t);
s_optu_4=sin(2*pi*(param.f_ocs+spacing*4)*t);

% Generate NRZ signal

% some useful values
NRZdatarate=param.NRZdatarate;
dt=1/(param.f_ocs*param.oversampling);


sam_sym=round(1/(NRZdatarate*dt));
numsymbol=round(param.noss*round(1/(param.BaudRate*dt))/sam_sym)+1;
slen=floor((param.noss*samples_symbol));
t=(0:slen-1)'*dt; % time vector

% Random data generation
data_NRZ=round(rand(numsymbol,1));
data_opt_NRZ=zeros(slen,1);


% Spread data to sample_symbol sampling points
for i=1:slen
    data_opt_NRZ(i)=data_NRZ(1+floor((i-1)/sam_sym));
end


% Pulse shaping
% Generate a frequency vector for filtering
f_vec=-slen/2:slen/2-1;
f_vec=(f_vec/slen)/dt;
f_vec=fftshift(f_vec);


% Gaussian filter
H_Gauss=exp(-log(2)/2*(f_vec/(param.BaudRate)).^2);
Fs_NRZ=fft(data_opt_NRZ);
Fs_NRZ=Fs_NRZ.*H_Gauss';


% Apply filter on signal
data_opt_NRZ=real(ifft(Fs_NRZ));
data_opt_NRZ=data_opt_NRZ/max(data_opt_NRZ); %normalize;



% Intensity modulation&mix
s_opt1=data_opt_NRZ.*s_optl_1;
s_mod=s_mod+circshift(s_opt1,decorrelate);
s_opt1=data_opt_NRZ.*s_optl_2;
s_mod=s_mod+circshift(s_opt1,decorrelate*2);
s_opt1=data_opt_NRZ.*s_optl_3;
s_mod=s_mod+circshift(s_opt1,decorrelate*3);
s_opt1=data_opt_NRZ.*s_optl_4;
s_mod=s_mod+circshift(s_opt1,decorrelate*4);

s_opt1=data_opt_NRZ.*s_optu_1;
s_mod=s_mod+circshift(s_opt1,-1*decorrelate);
s_opt1=data_opt_NRZ.*s_optu_2;
s_mod=s_mod+circshift(s_opt1,-1*decorrelate*2);
s_opt1=data_opt_NRZ.*s_optu_3;
s_mod=s_mod+circshift(s_opt1,-1*decorrelate*3);
s_opt1=data_opt_NRZ.*s_optu_4;
s_mod=s_mod+circshift(s_opt1,-1*decorrelate*4);

% figure(1)
% plot(data_opt_NRZ);
% figure(2)
% plot(s_optl_1);








