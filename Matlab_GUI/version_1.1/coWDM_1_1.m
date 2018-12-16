function varargout = coWDM_1_1(varargin)
% COWDM_1_1 MATLAB code for coWDM_1_1.fig
%      COWDM_1_1, by itself, creates a new COWDM_1_1 or raises the existing
%      singleton*.
%
%      H = COWDM_1_1 returns the handle to a new COWDM_1_1 or the handle to
%      the existing singleton*.
%
%      COWDM_1_1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COWDM_1_1.M with the given input arguments.
%
%      COWDM_1_1('Property','Value',...) creates a new COWDM_1_1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before coWDM_1_1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to coWDM_1_1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help coWDM_1_1

% Last Modified by GUIDE v2.5 21-Apr-2016 04:48:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @coWDM_1_1_OpeningFcn, ...
                   'gui_OutputFcn',  @coWDM_1_1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before coWDM_1_1 is made visible.
function coWDM_1_1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to coWDM_1_1 (see VARARGIN)

% Choose default command line output for coWDM_1_1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

I = imread('image_3.jpg');
axes(handles.axes1);
imshow(I);

CC=clock;

TT_text=[num2str(CC(1)) '/' num2str(CC(2)) '/' num2str(CC(3)) '  ' num2str(CC(4)) ':' num2str(CC(5)) 10 'Welcome ^_^' 10 'Set parameters and press Run!'];

set(handles.text17,'string',TT_text);
pause(0.5);

global MMMM;
MMMM=1;

% UIWAIT makes coWDM_1_1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = coWDM_1_1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MMMM;

% Initialize the figure
axes(handles.axes2);
plot(0,0,'w');
axes(handles.axes3);
plot(0,0,'w');
axes(handles.axes4);
plot(0,0,'w');

% Log creation
TT_text=['Simulation Begins'];

set(handles.text17,'string',TT_text);
pause(1);

% Control Parameter
ST_noise=str2double(get(handles.edit6,'String'));
WDM=str2double(get(handles.edit11,'String'));
ff=str2double(get(handles.edit14,'String'));
CMA=str2double(get(handles.edit15,'String'));

% input parameters
param=get_para;

param.spacing=str2double(get(handles.edit12,'String'));
param.spacing=param.spacing*1e9;
param.NRZdatarate=str2double(get(handles.edit13,'String'));
param.NRZdatarate=param.NRZdatarate*1e9;
param.BaudRate=str2double(get(handles.edit1,'String'));   % symbols per second
param.BaudRate=param.BaudRate*1e9;
param.noss=str2double(get(handles.edit2,'String')); % number of simulated symbols
param.oversampling=str2double(get(handles.edit3,'String')); % in optical domain
param.fiber_length=str2double(get(handles.edit4,'String'));
param.Pin=str2double(get(handles.edit5,'String'));

param.f_ocs=str2double(get(handles.edit7,'String'));
param.f_ocs=param.f_ocs*1e9;
param.f_LO=str2double(get(handles.edit8,'String'));
param.f_LO=param.f_LO*1e9;
temp_fwl=str2double(get(handles.edit9,'String'));
param.f_opt=3e8/(temp_fwl*1e-9);
param.line_width=str2double(get(handles.edit10,'String'));


TT_text=['Initializaiton' 10 '='];
for iii=1:20
TT_text=[TT_text '='];
set(handles.text17,'string',TT_text);
pause(0.1);
end


if(MMMM==1)
% QPSK
% some useful values
dt=1/(param.f_ocs*param.oversampling);
samples_symbol=round(1/(param.BaudRate*dt));
SymbolTime=1/param.BaudRate;
slen=floor((param.noss*samples_symbol));
t=(0:slen-1)'*dt; % time vector


% "Monte Carlo" optical carrier with phase noise (Lorenzian)
noise=rand(length(t),1)*pi;
linear_phase_error=cumsum(min(abs(tan(noise)),100000).*sign(tan(noise))) *dt*param.line_width*pi;
% linear_phase_error=0;
s_opt=sin(2*pi*param.f_ocs*t + linear_phase_error);
% s_opt=exp(1i*2*pi*param.f_ocs*t + linear_phase_error);


% Random data generation
data_I=round(rand(param.noss,1))-0.5;
data_Q=round(rand(param.noss,1))-0.5;
data_opt_I=zeros(slen,1);
data_opt_Q=zeros(slen,1);


% Spread data to sample_symbol sampling points
for i=1:slen
    data_opt_I(i)=data_I(1+floor((i-1)/samples_symbol));
    data_opt_Q(i)=data_Q(1+floor((i-1)/samples_symbol));
end


% Pulse shaping
% Generate a frequency vector for filtering
f_vec=-slen/2:slen/2-1;
f_vec=(f_vec/slen)/dt;
f_vec=fftshift(f_vec);


% Gaussian filter
H_Gauss=exp(-log(2)/2*(f_vec/(param.BaudRate)).^2);
Fs_I=fft(data_opt_I);
Fs_Q=fft(data_opt_Q);
Fs_I=Fs_I.*H_Gauss';
Fs_Q=Fs_Q.*H_Gauss';


% Apply filter on signal
data_opt_I=real(ifft(Fs_I));
data_opt_Q=real(ifft(Fs_Q));
data_opt_I=data_opt_I/max(data_opt_I); %normalize;
data_opt_Q=data_opt_Q/max(data_opt_Q); %normalize;


% create QAM signal on optical carrier
s_opt1=s_opt.*exp(1i*2*pi*param.Vpi*data_opt_I/2)-s_opt.*exp(-1i*2*pi*param.Vpi*data_opt_I/2);
s_opt2=s_opt.*exp(1i*2*pi*param.Vpi*data_opt_Q/2)-s_opt.*exp(-1i*2*pi*param.Vpi*data_opt_Q/2);
s_mod=s_opt1+1i*s_opt2;
% s_mod(round(slen/2))=1000;% 00arker:


% Adds 2 more modulated carriers to optical TX spectrum
if(WDM==1)
    more_carriers;
end

if(WDM==2)
    more_carriers_inten;
end

% Show optical TX spectrum
% figure(10001)
% OSA(s_mod,dt,1e9,param.f_ocs,'TX spectrum');

Dis_test=['Fiber Transmitting'];
set(handles.text17,'string',Dis_test);
pause(0.5);

% Dispersive nonlinear fiber
% s_mod=SSMF(s_mod,param.fiber_length,param.dispersion,dt,param.f_opt,param.Pin,param.alpha,param.gamma);
s_mod=SSMF_2(s_mod,param.fiber_length,param.dispersion,dt,param.f_opt,param.Pin,param.alpha,param.gamma);


% Show optical spectrum at fiber out
% figure(10002)
% OSA(s_mod,dt,1e9,param.f_ocs,'RX spectrum');
axes(handles.axes2)
% figure(10002)
% handles.axes1=OSA(s_mod,dt,1e8,param.f_ocs,'spectrum@RX');
% title('Spectrum@Rx GHz vs dB');
OSA(s_mod,dt,1e9,param.f_ocs,'Spectrum@RX');
title('Spectrum@Rx GHz vs dB');
% pause(1)


% Receiver with optical 90° Hybrid - e/o conversion in balanced RX
E_down_I=100*exp(2*pi*1j*(param.f_LO+param.IF_offset)*t); % optical LO 0°
E_down_Q=100*exp(2*pi*1j*(param.f_LO+param.IF_offset)*t+1j*pi/2); % optical LO 90°
s_I=abs((s_mod+E_down_I)).^2-abs((s_mod-E_down_I)).^2;
s_Q=abs((s_mod+E_down_Q)).^2-abs((s_mod-E_down_Q)).^2; % balanced RX


% Optional - Add the shot noise and thermal noise.
if(ST_noise==1)
    Rl=50;  %Load resistor/Unit omu
    Kb=1.380622e-23; %Boltzmann constant
    Tem=290;  %Temperature/Unit K
    Noise_BW=5e9;   %Effective noise bandwidth/Unit Hz
    q_EC=1.60218e-19;    %The elementary charge constant
    R=1;    %Response of PD
    Ip_I=R*mean(abs(s_I).^2);    %Average current/Unit A
    Ip_Q=R*mean(abs(s_Q).^2);
    
    sigma_short_I=sqrt(2*q_EC*Ip_I*Noise_BW);
    sigma_short_Q=sqrt(2*q_EC*Ip_Q*Noise_BW);
    sigma_thermal=sqrt(4*Kb*Tem*Noise_BW/Rl);
    
    sigma_t_I=sqrt(sigma_short_I^2+sigma_thermal^2);
    sigma_t_Q=sqrt(sigma_short_Q^2+sigma_thermal^2);
    
    s_I=s_I+sigma_t_I*randn(length(s_I),1)-sigma_t_I*randn(length(s_I),1);
    s_Q=s_Q+sigma_t_Q*randn(length(s_Q),1)-sigma_t_Q*randn(length(s_Q),1);
end


% FILTERING @ RX
% Bandwidth of PD << param.f_oc e. g. 2*f_ec
% Generate a frequency vector for filtering
f_vec=-slen/2:slen/2-1;
f_vec=(f_vec/slen)/dt;
f_vec=fftshift(f_vec);


% Filtering out optical carrier is modelled by simple rectangular filter
RectFilter=2/SymbolTime;
H_Rect=zeros(slen,1);
H_Rect(abs(f_vec)<=RectFilter)=1;
Fs_I=fft(s_I);
Fs_I=Fs_I.*H_Rect;
s_I=real(ifft(Fs_I));
s_I=s_I/max(s_I); %normalize; assuming s is offset-free
Fs_Q=fft(s_Q);
Fs_Q=Fs_Q.*H_Rect;
s_Q=real(ifft(Fs_Q));
s_Q=s_Q/max(s_Q); %normalize; assuming s is offset-free


% Low pass filtering down converted signals in electrical domain
% Raised cosine filter transfer function
atm=(1-param.roll_off)/(2*SymbolTime); %intermediate variable
atp=(1+param.roll_off)/(2*SymbolTime); %intermediate variable
pitro=(pi*SymbolTime)/(2*param.roll_off); %intermediate variable to simplyfy later expressions
hg=zeros(length(f_vec),1);   %transfer function
for i=1:length(f_vec)    %implement cos roll off filter
    absf=abs(f_vec(i));
    if(absf<=atm)
        hg(i)=1;
    end
    if ((atm<absf)&&(absf<atp))
        hg(i)=cos(pitro*(absf-atm));
    end
end


% Apply filter transfer function hg
fd_i=fft(s_I);
fd_q=fft(s_Q);
fd_i_filt=fd_i.*hg;
fd_q_filt=fd_q.*hg;
d_i=real((ifft(fd_i_filt)));
d_q=real((ifft(fd_q_filt)));
% End of filtering at RX


% Normalize the signal
d_i=d_i/max(d_i);
d_q=d_q/max(d_q);







%------------------------------Conventional DSP
% if (DSP==1)

Baudrate=param.BaudRate;
Samplingrate=samples_symbol*Baudrate;
rsps=2; % resample samples per symbol

% figure(301)
% plot(d_i,d_q,'+');
% axis([-1 1 -1 1]);
% grid on


T=1/Baudrate;
sps=samples_symbol; % samples per symbol
Tsamp=T/sps;
Tresamp=T/rsps;

Ix_rcvd=d_i;
Qx_rcvd=d_q;

% Resample the points
ts_Ix_rcvd=timeseries(Ix_rcvd, (1:1:length(Ix_rcvd))*Tsamp);
res_ts_Ix_rcvd=resample(ts_Ix_rcvd, (ts_Ix_rcvd.time(1):Tresamp:ts_Ix_rcvd.time(size(ts_Ix_rcvd.time, 1))));

ts_Qx_rcvd=timeseries(Qx_rcvd, (1:1:length(Qx_rcvd))*Tsamp);
res_ts_Qx_rcvd=resample(ts_Qx_rcvd, (ts_Qx_rcvd.time(1):Tresamp:ts_Qx_rcvd.time(size(ts_Qx_rcvd.time, 1))));

Ix_rcvd_rsamp = res_ts_Ix_rcvd.data;
Qx_rcvd_rsamp = res_ts_Qx_rcvd.data;

axes(handles.axes3)
% title('Constellation@Rx Before DSP');
plot(Ix_rcvd_rsamp,Qx_rcvd_rsamp,'+');
title('Constellation@Rx Before DSP');
axis([-1 1 -1 1]);
grid on

Ix_offset=1;
Qx_offset=1;
Ix_2sps=Ix_rcvd_rsamp(Ix_offset:rsps/2:end);
Qx_2sps=Qx_rcvd_rsamp(Qx_offset:rsps/2:end);

min_len_IQ=min([length(Ix_2sps), length(Qx_2sps)]);
Ix_2sps(min_len_IQ+1:end)=[];
Qx_2sps(min_len_IQ+1:end)=[];

Ix_2sps=Ix_2sps-mean(Ix_2sps);
Qx_2sps=Qx_2sps-mean(Qx_2sps);
Ix_normalized_2sps_1=Ix_2sps./(sqrt(mean(abs(Ix_2sps).^2)));
Qx_normalized_2sps_1=Qx_2sps./(sqrt(mean(abs(Qx_2sps).^2)));

Ix_normalized_2sps=Ix_normalized_2sps_1./(sqrt(mean(abs(Ix_normalized_2sps_1).^2+abs(Qx_normalized_2sps_1).^2)));
Qx_normalized_2sps=Qx_normalized_2sps_1./(sqrt(mean(abs(Ix_normalized_2sps_1).^2+abs(Qx_normalized_2sps_1).^2)));

% figure(303)
% plot(Ix_normalized_2sps,Qx_normalized_2sps,'+');
% axis([-2 2 -2 2]);
% grid on

data=Ix_normalized_2sps+1i*Qx_normalized_2sps;

% Test for frecomp
% [yout_fc, yout_pc, df]=frecomp(data', 1/(Tsamp/2), 1, 10);
% data = yout_fc;


% Pol_Demux (CMA filter)
if(CMA==1)
    
R=[1 0; 1 0];
mu=0.001;
taps=13;
% x=[data(1:length(data)/2) data(1:length(data)/2)];
x=[data data];
[yout, deth]=pol_demux(x,taps, mu, R);    
data=yout(:, 1);

% figure(303)
% plot(real(data),imag(data),'+');
% axis([-2 2 -2 2]);
% grid on

end

% ff_recovery
if(ff==1)

constell=[1+1i*1;-1+1i*1;1-1i*1;-1-1i*1;];
% constell=[1+1i*1;-1+1i*1;1-1i*1;-1-1i*1;3+1i*1;-3+1i*1;3-1i*1;-3-1i*1;1+1i*3;-1+1i*3;1-1i*3;-1-1i*3;3+1i*3;-3+1i*3;3-1i*3;-3-1i*3;];
constell=constell/sqrt(mean(abs(constell).^2));
data=data/sqrt(mean(abs(data).^2));
% Signal_x_Equalized=ff_carrier_recovery_mQAM(data, constell);
data=ff_carrier_recovery_mQAM(data, constell);

% Signal_x_Equalized=data;
% 
% figure(304)
% plot(real(Signal_x_Equalized),imag(Signal_x_Equalized),'+');
% axis([-2 2 -2 2]);
% grid on

end

axes(handles.axes4)
% title('Constellation@Rx After DSP');
plot(real(data),imag(data),'+');
title('Constellation@Rx After DSP');
axis([-2 2 -2 2]);
grid on

% Decision and Q calculation
constell=[1+1i*1;-1+1i*1;1-1i*1;-1-1i*1;];
% constell=[1+1i*1;-1+1i*1;1-1i*1;-1-1i*1;3+1i*1;-3+1i*1;3-1i*1;-3-1i*1;1+1i*3;-1+1i*3;1-1i*3;-1-1i*3;3+1i*3;-3+1i*3;3-1i*3;-3-1i*3;];
constell=constell/sqrt(mean(abs(constell).^2));
constell_size=length(constell);
Signal_x_Equalized_mat=kron(data.', ones(constell_size, 1));
Est_sym_x_rcvd_mat=kron(constell, ones(1, length(data)));

[Est_sym_x_val, Est_sym_x_ind]=min(abs(Est_sym_x_rcvd_mat-Signal_x_Equalized_mat).^2);
src_sym_x=constell(Est_sym_x_ind);

sigma2=mean((abs(data - src_sym_x)).^2);
sig_pwr=mean((abs(data).^2));
Qx=10*log10(sig_pwr/sigma2);
[Ix_offset Qx]

Dis_test=['DSP Done'];
set(handles.text17,'string',Dis_test);
pause(0.5);


% Text=get(handles.text17,'String');
Dis_test=['Simulation Done' 10 'Q factor=' num2str(Qx)];
set(handles.text17,'string',Dis_test);

  
end


%-----------------------------------------QAM

if(MMMM==2)
% QAM
% some useful values
dt=1/(param.f_ocs*param.oversampling);
samples_symbol=round(1/(param.BaudRate*dt));
SymbolTime=1/param.BaudRate;
slen=floor((param.noss*samples_symbol));
t=(0:slen-1)'*dt; % time vector


% "Monte Carlo" optical carrier with phase noise (Lorenzian)
noise=rand(length(t),1)*pi;
linear_phase_error=cumsum(min(abs(tan(noise)),100000).*sign(tan(noise))) *dt*param.line_width*pi;
% linear_phase_error=0;
s_opt=sin(2*pi*param.f_ocs*t + linear_phase_error);
% s_opt=exp(1i*2*pi*param.f_ocs*t + linear_phase_error);


% Random data generation
% data_I=round(rand(param.noss,1))-0.5;
% data_Q=round(rand(param.noss,1))-0.5;
data_opt_I=zeros(slen,1);
data_opt_Q=zeros(slen,1);
% For 16-QAM
data_I=rand(param.noss,1);
data_Q=rand(param.noss,1);
for i=1:1:param.noss
    if (data_I(i)<0.25)
        data_I(i)=0.2;
    elseif (data_I(i)>=0.25)&&(data_I(i)<0.5)
        data_I(i)=-0.2;
    elseif (data_I(i)>=0.5)&&(data_I(i)<0.75)
        data_I(i)=0.6;
    elseif (data_I(i)>=0.75)
        data_I(i)=-0.6;
    end
    
    if (data_Q(i)<0.25)
        data_Q(i)=0.2;
    elseif (data_Q(i)>=0.25)&&(data_Q(i)<0.5)
        data_Q(i)=-0.2;
    elseif (data_Q(i)>=0.5)&&(data_Q(i)<0.75)
        data_Q(i)=0.6;
    elseif (data_Q(i)>=0.75)
        data_Q(i)=-0.6;
    end
end


% Spread data to sample_symbol sampling points
for i=1:slen
    data_opt_I(i)=data_I(1+floor((i-1)/samples_symbol));
    data_opt_Q(i)=data_Q(1+floor((i-1)/samples_symbol));
end


% Pulse shaping
% Generate a frequency vector for filtering
f_vec=-slen/2:slen/2-1;
f_vec=(f_vec/slen)/dt;
f_vec=fftshift(f_vec);


% Gaussian filter
H_Gauss=exp(-log(2)/2*(f_vec/(param.BaudRate)).^2);
Fs_I=fft(data_opt_I);
Fs_Q=fft(data_opt_Q);
Fs_I=Fs_I.*H_Gauss';
Fs_Q=Fs_Q.*H_Gauss';


% Apply filter on signal
data_opt_I=real(ifft(Fs_I));
data_opt_Q=real(ifft(Fs_Q));
data_opt_I=data_opt_I/max(data_opt_I); %normalize;
data_opt_Q=data_opt_Q/max(data_opt_Q); %normalize;


% create QAM signal on optical carrier
s_opt1=s_opt.*exp(1i*2*pi*param.Vpi*data_opt_I/2)-s_opt.*exp(-1i*2*pi*param.Vpi*data_opt_I/2);
s_opt2=s_opt.*exp(1i*2*pi*param.Vpi*data_opt_Q/2)-s_opt.*exp(-1i*2*pi*param.Vpi*data_opt_Q/2);
s_mod=s_opt1+1i*s_opt2;
% s_mod(round(slen/2))=1000;% 00arker:


% Adds 2 more modulated carriers to optical TX spectrum
if(WDM==1)
    more_carriers;
end

if(WDM==2)
    more_carriers_inten;
end

Dis_test=['Fiber Transmitting'];
set(handles.text17,'string',Dis_test);
pause(0.5);

% Show optical TX spectrum
% figure(10001)
% OSA(s_mod,dt,1e8,param.f_ocs,'TX spectrum');


% Dispersive nonlinear fiber
% s_mod=SSMF(s_mod,param.fiber_length,param.dispersion,dt,param.f_opt,param.Pin,param.alpha,param.gamma);
s_mod=SSMF_2(s_mod,param.fiber_length,param.dispersion,dt,param.f_opt,param.Pin,param.alpha,param.gamma);


% Show optical spectrum at fiber out
axes(handles.axes2)
% figure(10002)
% handles.axes1=OSA(s_mod,dt,1e8,param.f_ocs,'spectrum@RX');
% title('Spectrum@Rx GHz vs dB');
OSA(s_mod,dt,1e9,param.f_ocs,'Spectrum@RX');
title('Spectrum@Rx GHz vs dB');
% pause(1)
% pause(1)


% Receiver with optical 90° Hybrid - e/o conversion in balanced RX
E_down_I=100*exp(2*pi*1j*(param.f_LO+param.IF_offset)*t); % optical LO 0°
E_down_Q=100*exp(2*pi*1j*(param.f_LO+param.IF_offset)*t+1j*pi/2); % optical LO 90°
s_I=abs((s_mod+E_down_I)).^2-abs((s_mod-E_down_I)).^2;
s_Q=abs((s_mod+E_down_Q)).^2-abs((s_mod-E_down_Q)).^2; % balanced RX


% Optional - Add the shot noise and thermal noise.
if(ST_noise==1)
    Rl=50;  %Load resistor/Unit omu
    Kb=1.380622e-23; %Boltzmann constant
    Tem=290;  %Temperature/Unit K
    Noise_BW=5e9;   %Effective noise bandwidth/Unit Hz
    q_EC=1.60218e-19;    %The elementary charge constant
    R=1;    %Response of PD
    Ip_I=R*mean(abs(s_I).^2);    %Average current/Unit A
    Ip_Q=R*mean(abs(s_Q).^2);
    
    sigma_short_I=sqrt(2*q_EC*Ip_I*Noise_BW);
    sigma_short_Q=sqrt(2*q_EC*Ip_Q*Noise_BW);
    sigma_thermal=sqrt(4*Kb*Tem*Noise_BW/Rl);
    
    sigma_t_I=sqrt(sigma_short_I^2+sigma_thermal^2);
    sigma_t_Q=sqrt(sigma_short_Q^2+sigma_thermal^2);
    
    s_I=s_I+sigma_t_I*randn(length(s_I),1)-sigma_t_I*randn(length(s_I),1);
    s_Q=s_Q+sigma_t_Q*randn(length(s_Q),1)-sigma_t_Q*randn(length(s_Q),1);
end


% FILTERING @ RX
% Bandwidth of PD << param.f_oc e. g. 2*f_ec
% Generate a frequency vector for filtering
f_vec=-slen/2:slen/2-1;
f_vec=(f_vec/slen)/dt;
f_vec=fftshift(f_vec);


% Filtering out optical carrier is modelled by simple rectangular filter
RectFilter=2/SymbolTime;
H_Rect=zeros(slen,1);
H_Rect(abs(f_vec)<=RectFilter)=1;
Fs_I=fft(s_I);
Fs_I=Fs_I.*H_Rect;
s_I=real(ifft(Fs_I));
s_I=s_I/max(s_I); %normalize; assuming s is offset-free
Fs_Q=fft(s_Q);
Fs_Q=Fs_Q.*H_Rect;
s_Q=real(ifft(Fs_Q));
s_Q=s_Q/max(s_Q); %normalize; assuming s is offset-free


% Low pass filtering down converted signals in electrical domain
% Raised cosine filter transfer function
atm=(1-param.roll_off)/(2*SymbolTime); %intermediate variable
atp=(1+param.roll_off)/(2*SymbolTime); %intermediate variable
pitro=(pi*SymbolTime)/(2*param.roll_off); %intermediate variable to simplyfy later expressions
hg=zeros(length(f_vec),1);   %transfer function
for i=1:length(f_vec)    %implement cos roll off filter
    absf=abs(f_vec(i));
    if(absf<=atm)
        hg(i)=1;
    end
    if ((atm<absf)&&(absf<atp))
        hg(i)=cos(pitro*(absf-atm));
    end
end


% Apply filter transfer function hg
fd_i=fft(s_I);
fd_q=fft(s_Q);
fd_i_filt=fd_i.*hg;
fd_q_filt=fd_q.*hg;
d_i=real((ifft(fd_i_filt)));
d_q=real((ifft(fd_q_filt)));
% End of filtering at RX


% Normalize the signal
d_i=d_i/max(d_i);
d_q=d_q/max(d_q);

% Re_Sam=round(samples_symbol/2):samples_symbol:round(param.noss*samples_symbol-samples_symbol/2);


% Constellation Diagram of direct I and Q 
% ff-carrier recovery
% Rin=d_i(Re_Sam)+1i*d_q(Re_Sam);
% Rin=Rin/sqrt(mean(abs(Rin).^2));
% Const_pr=[1+1i*1;-1+1i*1;1-1i*1;-1-1i*1;];
% Const_pr=Const_pr/sqrt(mean(abs(Const_pr).^2));
% Rout=ff_carrier_recovery_mQAM(Rin,Const_pr);
% Rout=Rout/sqrt(mean(abs(Rout).^2));
% figure(100)
% plot(real(Rout),imag(Rout),'+');
% axis([-1 1 -1 1]);
% grid on
% figure(101)
% plot(d_i(Re_Sam),d_q(Re_Sam),'+');


%-------------------------------------DSP for the mQAM 
Dis_test=['DSP Begins'];
set(handles.text17,'string',Dis_test);
pause(0.5);

Baudrate=param.BaudRate;
Samplingrate=samples_symbol*Baudrate;
rsps=2; % resample samples per symbol

% figure(301)
% plot(d_i,d_q,'+');
% axis([-1 1 -1 1]);
% grid on


T=1/Baudrate;
sps=samples_symbol; % samples per symbol
Tsamp=T/sps;
Tresamp=T/rsps;

Ix_rcvd=d_i;
Qx_rcvd=d_q;

% Resample the points
ts_Ix_rcvd=timeseries(Ix_rcvd, (1:1:length(Ix_rcvd))*Tsamp);
res_ts_Ix_rcvd=resample(ts_Ix_rcvd, (ts_Ix_rcvd.time(1):Tresamp:ts_Ix_rcvd.time(size(ts_Ix_rcvd.time, 1))));

ts_Qx_rcvd=timeseries(Qx_rcvd, (1:1:length(Qx_rcvd))*Tsamp);
res_ts_Qx_rcvd=resample(ts_Qx_rcvd, (ts_Qx_rcvd.time(1):Tresamp:ts_Qx_rcvd.time(size(ts_Qx_rcvd.time, 1))));

Ix_rcvd_rsamp = res_ts_Ix_rcvd.data;
Qx_rcvd_rsamp = res_ts_Qx_rcvd.data;

% figure(302)
% plot(Ix_rcvd_rsamp,Qx_rcvd_rsamp,'+');
% axis([-1 1 -1 1]);
% grid on

axes(handles.axes3)
% title('Constellation@Rx Before DSP');
plot(Ix_rcvd_rsamp,Qx_rcvd_rsamp,'+');
title('Constellation@Rx Before DSP');
axis([-1 1 -1 1]);
grid on

Ix_offset=1;
Qx_offset=1;
Ix_2sps=Ix_rcvd_rsamp(Ix_offset:rsps/2:end);
Qx_2sps=Qx_rcvd_rsamp(Qx_offset:rsps/2:end);

min_len_IQ=min([length(Ix_2sps), length(Qx_2sps)]);
Ix_2sps(min_len_IQ+1:end)=[];
Qx_2sps(min_len_IQ+1:end)=[];

Ix_2sps=Ix_2sps-mean(Ix_2sps);
Qx_2sps=Qx_2sps-mean(Qx_2sps);
Ix_normalized_2sps_1=Ix_2sps./(sqrt(mean(abs(Ix_2sps).^2)));
Qx_normalized_2sps_1=Qx_2sps./(sqrt(mean(abs(Qx_2sps).^2)));

Ix_normalized_2sps=Ix_normalized_2sps_1./(sqrt(mean(abs(Ix_normalized_2sps_1).^2+abs(Qx_normalized_2sps_1).^2)));
Qx_normalized_2sps=Qx_normalized_2sps_1./(sqrt(mean(abs(Ix_normalized_2sps_1).^2+abs(Qx_normalized_2sps_1).^2)));

% figure(303)
% plot(Ix_normalized_2sps,Qx_normalized_2sps,'+');
% axis([-2 2 -2 2]);
% grid on

data=Ix_normalized_2sps+1i*Qx_normalized_2sps;

% Test for frecomp
% [yout_fc, yout_pc, df]=frecomp(data', 1/(Tsamp/2), 1, 10);
% data = yout_fc;


% Pol_Demux (CMA filter)
if(CMA==1)

Dis_test=['CMA Begins'];
set(handles.text17,'string',Dis_test);
pause(0.5);

    
R=[1 0; 1 0];
mu=0.001;
taps=13;
% x=[data(1:length(data)/2) data(1:length(data)/2)];
x=[data data];
[yout, deth]=pol_demux(x,taps, mu, R);    
data=yout(:, 1);

% figure(303)
% plot(real(data),imag(data),'+');
% axis([-2 2 -2 2]);
% grid on

end

% ff_recovery
if(ff==1)

% constell=[1+1i*1;-1+1i*1;1-1i*1;-1-1i*1;];
constell=[1+1i*1;-1+1i*1;1-1i*1;-1-1i*1;3+1i*1;-3+1i*1;3-1i*1;-3-1i*1;1+1i*3;-1+1i*3;1-1i*3;-1-1i*3;3+1i*3;-3+1i*3;3-1i*3;-3-1i*3;];
constell=constell/sqrt(mean(abs(constell).^2));
data=data/sqrt(mean(abs(data).^2));
data=ff_carrier_recovery_mQAM(data, constell);

% Signal_x_Equalized=data;

% figure(304)
% plot(real(Signal_x_Equalized),imag(Signal_x_Equalized),'+');
% axis([-2 2 -2 2]);
% grid on

end

% figure(305)
% plot(real(Signal_x_Equalized),imag(Signal_x_Equalized),'+');
% axis([-2 2 -2 2]);
% grid on

axes(handles.axes4)
% title('Constellation@Rx After DSP');
plot(real(data),imag(data),'+');
title('Constellation@Rx After DSP');
axis([-2 2 -2 2]);
grid on

% Decision and Q calculation
% constell=[1+1i*1;-1+1i*1;1-1i*1;-1-1i*1;];
constell=[1+1i*1;-1+1i*1;1-1i*1;-1-1i*1;3+1i*1;-3+1i*1;3-1i*1;-3-1i*1;1+1i*3;-1+1i*3;1-1i*3;-1-1i*3;3+1i*3;-3+1i*3;3-1i*3;-3-1i*3;];
constell=constell/sqrt(mean(abs(constell).^2));
constell_size=length(constell);
Signal_x_Equalized_mat=kron(data.', ones(constell_size, 1));
Est_sym_x_rcvd_mat=kron(constell, ones(1, length(data)));

[Est_sym_x_val, Est_sym_x_ind]=min(abs(Est_sym_x_rcvd_mat-Signal_x_Equalized_mat).^2);
src_sym_x=constell(Est_sym_x_ind);

sigma2=mean((abs(data - src_sym_x)).^2);
sig_pwr=mean((abs(data).^2));
Qx=10*log10(sig_pwr/sigma2);
[Ix_offset Qx]

Dis_test=['DSP Done'];
set(handles.text17,'string',Dis_test);
pause(0.5);


% Text=get(handles.text7,'String');
Dis_test=['Simulation Done' 10 'Q factor=' num2str(Qx)];
set(handles.text17,'string',Dis_test);

end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
global MMMM;

contents = get(handles.popupmenu1,'String'); 
popupmenu4value = contents{get(handles.popupmenu1,'Value')};
switch popupmenu4value
   case 'Coherent_homodyne_QPSK'
        %function of A
%         I = imread('image_1.jpg');
%         axes(handles.axes4);
%         imshow(I);
        MMMM=1;
   case 'Coherent_homodyne_QAM'
        %function of B
%         I = imread('image_2.jpg');
%         axes(handles.axes4);
%         imshow(I);
        MMMM=2;
end





% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
