% more_carriers

%this code just adds two more carriers in order to demonstrate the
%nonlinear fiber properties.

spacing=param.spacing;
decorrelate=samples_symbol*20;

s_optl_1=sin(2*pi*(param.f_ocs-spacing)*t);
s_optl_2=sin(2*pi*(param.f_ocs-spacing*2)*t);
s_optl_3=sin(2*pi*(param.f_ocs-spacing*3)*t);
s_optl_4=sin(2*pi*(param.f_ocs-spacing*4)*t);

s_optu_1=sin(2*pi*(param.f_ocs+spacing)*t);
s_optu_2=sin(2*pi*(param.f_ocs+spacing*2)*t);
s_optu_3=sin(2*pi*(param.f_ocs+spacing*3)*t);
s_optu_4=sin(2*pi*(param.f_ocs+spacing*4)*t);

% create 8 more QAM signals on optical carrier
s_opt1=s_optl_1.*exp(1i*2*pi*param.Vpi*data_opt_I/2)-s_optl_1.*exp(-1i*2*pi*param.Vpi*data_opt_I/2);
s_opt2=s_optl_1.*exp(1i*2*pi*param.Vpi*data_opt_Q/2)-s_optl_1.*exp(-1i*2*pi*param.Vpi*data_opt_Q/2);
s_mod=s_mod+circshift(s_opt1+1i*s_opt2,decorrelate); %decorrelation of data

s_opt1=s_optl_2.*exp(1i*2*pi*param.Vpi*data_opt_I/2)-s_optl_2.*exp(-1i*2*pi*param.Vpi*data_opt_I/2);
s_opt2=s_optl_2.*exp(1i*2*pi*param.Vpi*data_opt_Q/2)-s_optl_2.*exp(-1i*2*pi*param.Vpi*data_opt_Q/2);
s_mod=s_mod+circshift(s_opt1+1i*s_opt2,decorrelate*2);

s_opt1=s_optl_3.*exp(1i*2*pi*param.Vpi*data_opt_I/2)-s_optl_3.*exp(-1i*2*pi*param.Vpi*data_opt_I/2);
s_opt2=s_optl_3.*exp(1i*2*pi*param.Vpi*data_opt_Q/2)-s_optl_3.*exp(-1i*2*pi*param.Vpi*data_opt_Q/2);
s_mod=s_mod+circshift(s_opt1+1i*s_opt2,decorrelate*3);

s_opt1=s_optl_4.*exp(1i*2*pi*param.Vpi*data_opt_I/2)-s_optl_4.*exp(-1i*2*pi*param.Vpi*data_opt_I/2);
s_opt2=s_optl_4.*exp(1i*2*pi*param.Vpi*data_opt_Q/2)-s_optl_4.*exp(-1i*2*pi*param.Vpi*data_opt_Q/2);
s_mod=s_mod+circshift(s_opt1+1i*s_opt2,decorrelate*4);

s_opt1=s_optu_1.*exp(1i*2*pi*param.Vpi*data_opt_I/2)-s_optu_1.*exp(-1i*2*pi*param.Vpi*data_opt_I/2);
s_opt2=s_optu_1.*exp(1i*2*pi*param.Vpi*data_opt_Q/2)-s_optu_1.*exp(-1i*2*pi*param.Vpi*data_opt_Q/2);
s_mod=s_mod+circshift(s_opt1+1i*s_opt2,-1*decorrelate); %decorrelation of data

s_opt1=s_optu_2.*exp(1i*2*pi*param.Vpi*data_opt_I/2)-s_optu_2.*exp(-1i*2*pi*param.Vpi*data_opt_I/2);
s_opt2=s_optu_2.*exp(1i*2*pi*param.Vpi*data_opt_Q/2)-s_optu_2.*exp(-1i*2*pi*param.Vpi*data_opt_Q/2);
s_mod=s_mod+circshift(s_opt1+1i*s_opt2,-1*decorrelate*2);

s_opt1=s_optu_3.*exp(1i*2*pi*param.Vpi*data_opt_I/2)-s_optu_3.*exp(-1i*2*pi*param.Vpi*data_opt_I/2);
s_opt2=s_optu_3.*exp(1i*2*pi*param.Vpi*data_opt_Q/2)-s_optu_3.*exp(-1i*2*pi*param.Vpi*data_opt_Q/2);
s_mod=s_mod+circshift(s_opt1+1i*s_opt2,-1*decorrelate*3);

s_opt1=s_optu_4.*exp(1i*2*pi*param.Vpi*data_opt_I/2)-s_optu_4.*exp(-1i*2*pi*param.Vpi*data_opt_I/2);
s_opt2=s_optu_4.*exp(1i*2*pi*param.Vpi*data_opt_Q/2)-s_optu_4.*exp(-1i*2*pi*param.Vpi*data_opt_Q/2);
s_mod=s_mod+circshift(s_opt1+1i*s_opt2,-1*decorrelate*4);







