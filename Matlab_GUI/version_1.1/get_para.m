function param = get_para

param.spacing=50e9;
param.NRZdatarate=50e9;
param.BaudRate=28e9; % symbols per second (The symbol rate is 5G/s, bit rate is 10G/s with QPSK)
param.noss=2^10; % number of simulated symbols
param.oversampling=10; % in optical domain

% needed for demodulation & BER estimation 
param.CNR=72; % CNR e. g. from noise Tool  
param.ave_viterbi=8; % averaging over ... symbols for carrier phase recovery
param.IF_offset=10e6; % mismatch between optical IF and electrical IF
% param.IF_offset=0;
param.line_width=1e5; % laser linewidth

param.f_opt=195e12; % "real" optical carrier 
param.f_ocs=500e9; % frequency of "optical" carrier used in simulation
% param.f_ocs=195e12;
param.f_LO=500e9; % frequency optical LO in simulation
% param.f_LO=100.5e9;

param.fiber_length=50000; % m 
% param.dispersion=17*1e-6; % s/m^2 fiber dispersion for SSMF (D=17 ps/nm/km)
param.dispersion=16.972*1e-6;
% param.Pin=0.1e-3; % total power, Watt
param.Pin=1e-3;

param.Vpi=0.25; % modulation index for IQ modulator (nested MZM)
param.roll_off=0.5;  % raised cosine filtering at RX (in baseband)

% nonlinear effects, attenuation
% typical alfa: 5e-5 (Lnl=20 km)
% typical gamma: 3e-3   rad/W/m
param.alpha=5.06e-5;
param.gamma=2e-3;
end

