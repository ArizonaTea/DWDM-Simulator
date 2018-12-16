function [ Rout ] = ff_carrier_recovery_mQAM( Rx, constell )
%FF_CARRIER_RECOVERY_MQAM Summary of this function goes here
%   Detailed explanation goes here
% Rx: complex vector input 
% constell: normalized complex mQAM constellation map

    numberOfBits = length(Rx);
    B = 20;
    N = 30;   % 6~10

    b = -B+1:1:B-1;
    phi_b = exp(1i*b/B*1/4*pi);

    Y = [kron(Rx(1), ones(N, 1)); Rx; kron(Rx(numberOfBits), ones(N, 1))]; % Y = [Rx(numberOfBits-N+1:numberOfBits); Rx; Rx(1:N)];% received samples
    S = kron(Y, phi_b);

    constell_tmp = kron(constell, ones(1, length(Y)));

    for i = 1:length(phi_b)
        S_tmp = kron(ones(length(constell), 1), transpose(S(:, i)));
        S(:, i) = transpose(min((abs(S_tmp - constell_tmp)).^2));
    end

    phi_rec = zeros(size(Rx));   % phase noise recovery
    for i = N+1 : N+numberOfBits
        [val pos] = min(sum(S(i-N:i+N, :)));
        phi_rec(i-N) = phi_b(pos);
    end

    Rout = Rx.*phi_rec;
    %figure;
    %plot(Rout, '.')

end

