function [outvals] = BrayCurtisIndex(tempWaves,lead)

% Bray-Curtis similarity index (Lian et al., 2010: Signal Process 90:684-8)

wvs = tempWaves;

fname = strcat('Lead_',num2str(lead));

wb = waitbar(0,'Calculating...','name',fname);


outvals = zeros(length(tempWaves),1);
for i = 1:length(wvs)
    
    waveTemp = repmat(wvs(:,i),[1,length(wvs)]);
    
    tempInd = 1 - (sum(abs(waveTemp - wvs)) ./ sum(abs(waveTemp) + abs(wvs)));

    outvals(i) = mean(tempInd);
    
    waitbar(i/length(tempWaves),wb,sprintf('%d / %d',i,length(tempWaves)));

end

close(wb)


