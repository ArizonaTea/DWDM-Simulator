function [yout deth]=pol_demux(x,taps, mu, R)
%%%%EXP:
%input:
% x=[datax datay]; %datax and datay are data from coherent receiver
% ntaps=7; % CMA tap number
% mu=1/2000; % CMA coefficiency
% R=[1 1]; % CMA convergence radius
%
%output:
%y is the demuxed data
%deth is error function ,can be used to analyze minimum number of taps
%%%%%
x(:,1)=en_nor(x(:,1),sqrt(1)); %normalization for QPSK X POL
x(:,2)=en_nor(x(:,2),sqrt(1)); %normalization for QPSK Y POL
halftaps= floor( taps / 2);          % half of (taps-1)
hzero = zeros( taps , 2, 2);
phi = 0;% initial phi
M = [ cos(phi) sin(phi); -sin(phi) cos(phi) ]; 
hzero( halftaps+1, :, :) = M;  % Initializing central taps
if halftaps
    extendedx = [ x(end-halftaps+1:end,:) ; x ; x(1:halftaps,:) ];
else
    extendedx = x;
end
h1 = squeeze( hzero(:, 1, :) );
h2 = squeeze( hzero(:, 2, :) );
if taps==1
    h1 = h1.';
    h2 = h2.';
end
convergence = false;
c=1; %iternation counter
% L=length(x);
repetitions = 50; %50*ceil(1./(L.*mu));
deth=[]; deth_temp=1000; %large initial
disp('Begin CMA iteration...')
while ~convergence && (c<(repetitions))

    
    h1_old = h1*1;
    h2_old = h2*1;
    [yout h1_new h2_new] = cmaadaptivefilter(extendedx, h1, h2, taps, mu, R, 2);
    
    if any(any(h1_new)) || any(any(h2_new))
        h1=h1_new;
        h2=h2_new;
    end
    deth1=max(max(abs([h1_old-h1 h2_old-h2])));
    
    deth=[deth deth1];
   
    if deth1 < 5e-7 || abs(deth1-deth_temp)<1e-8 % max(max(abs([(h1_old-h1)./(h1+h1_old)/2 (h2_old-h2)./(h1+h1_old)/2])))
        convergence = true;
    else
        deth_temp=deth1;
    end
    c=c+1;
end
disp ('Done')
end

function x=en_nor(x,e_n)
x=complex(real(x)-mean(real(x)),imag(x)-mean(imag(x)));
x=x*e_n./mean(abs(x));
end
