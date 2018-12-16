function out = SignalQuality(diffangle,BER_symbol,samples_symbol)

figure(99)
blt=1:20:2*samples_symbol;
ctt=1;
nhist=(1:360);

peye=zeros(length(blt),length(nhist));

symbol_0=[315:360 1:44]; % degree
symbol_90=[45:134];
symbol_180=[135:224];
symbol_270=[225:314];

%vary the sampling phase
BER_system=1;
position=1;

for start=(1:20:2*samples_symbol)+2*samples_symbol
    samplevec=start:samples_symbol:length(diffangle)-2*samples_symbol;  % a vector with the indices of the sampling point
    samplevec=round(floor(samplevec));  % and make sure its all integers!
    da=mod(diffangle(samplevec),2*pi)*180/pi; %wrap around
    n=hist(da,nhist); % do a histogram of the angle distribution
    BER_0=sum(BER_symbol.*n(symbol_0))/sum(n(symbol_0));
    BER_90=sum(BER_symbol.*n(symbol_90))/sum(n(symbol_90));
    BER_180=sum(BER_symbol.*n(symbol_180))/sum(n(symbol_180));
    BER_270=sum(BER_symbol.*n(symbol_270))/sum(n(symbol_270));
    BER=(BER_0+BER_90+BER_180+BER_270)/4;
    if BER_system>BER
        position=start;
        best_da=da;
        best_hist=nhist;
        best_samplevec=samplevec;
    end
    BER_system=min(BER_system,BER);
    peye(ctt,:)=n;  %save it into peye
    ctt=ctt+1;
    % if running in Octave, comment the following three lines, else the
    % run-time is too long
%     hist(da,nhist)
%     title(['BER:= ' num2str(BER,'%6.0e')],'fontsize',16)
%     pause(0.1); %in order to give time to plot the figure; without pause, only the last figure is plotted
end

    %Plot the best BER figure
    hist(best_da,best_hist)
    title(['BER:= ' num2str(BER_system,'%6.0e')],'fontsize',16)

out.samplevec=best_samplevec;
out.peye=peye;
out.ctt=ctt;
out.BER_system=BER_system;
out.sampling=position;
end