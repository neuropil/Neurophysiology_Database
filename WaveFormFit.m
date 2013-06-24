function [WaveFitParams, WaveSum_DS] = WaveFormFit(disindex,Clust_Waves)

time = linspace(0,1,32);
time = round(time * 100)/100;
sampleindex = (1:1:32)';

WaveFitParams = struct;
WaveSum_DS = struct;
for wS = 1:length(disindex)
    
    tempWaves = squeeze(Clust_Waves(:,disindex(wS),:));

    wavemeanTemp = mean(tempWaves,2);
    wavestdTemp = std(tempWaves,0,2);
    wavestdPTemp = wavemeanTemp + wavestdTemp;
    wavestdMTemp = wavemeanTemp - wavestdTemp;
    
    waveAll = horzcat(wavemeanTemp, wavestdTemp, wavestdPTemp, wavestdMTemp);
    
    WaveSum_DS.(strcat('L',num2str(disindex(wS)))) = mat2dataset(waveAll, 'VarNames', {'Mean','StD','PStD','MStD'});

    uselead = wavemeanTemp;
    
    basesub = mean(uselead(30:32));
    
    uselead_base = uselead - basesub;
    
    [~, maxpoint] = max(uselead_base);
    [~, minpoint] = min(uselead_base(8:32));
    
    
    minpoint = minpoint + 7;
    
    if maxpoint == 1
        maxpoint = 8;
    end
    
    before_max = find(sampleindex < maxpoint);
    after_max = find(sampleindex >= maxpoint);
    
    if after_max == 32
        after_max = (10:32)';
    end
    
    first_pos_def = find(uselead_base(before_max) < 0,1,'last');
    second_pos_def = find(uselead_base(after_max) < 0,1,'first');
    second_pos_def = second_pos_def - 1;
    
    if isempty(first_pos_def) == 1;
        first_pos_def = 1;
    end
    
    if second_pos_def == 0;
        second_pos_def = 10;
    end
    
    pospeakstart = uselead_base(before_max(first_pos_def));
    pospeakend = uselead_base(after_max(second_pos_def));
    
    posps_loc = find(uselead_base == pospeakstart);
    pospe_loc = find(uselead_base == pospeakend);
    
    pospeak_width = time(pospe_loc) - time(posps_loc);
    
    % min calculation
    
    before_min = find(sampleindex < minpoint);
    first_neg_def = find(uselead_base(before_min) > 0,1,'last');
    
    if isempty(first_neg_def) == 1;
        first_neg_def = max(before_min);
    end
    
    negpeakstart = uselead_base(before_min(first_neg_def));
    negps_loc = find(uselead_base == negpeakstart);
    
    halfwidth = ceil((minpoint + negps_loc)/2);
    addhalf = minpoint - halfwidth;
    halfwidthpoint2 = minpoint + addhalf;
    
    if halfwidthpoint2 >= 32
        halfwidthpoint2 = 32;
    end
    
    negpeak_half_width = time(halfwidthpoint2) - time(halfwidth);
    
    baselinevolt = abs(basesub);
    spkamp = max(uselead_base) - baselinevolt;
    
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).jtkoyama_neg_width = negpeak_half_width;
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).jtkoyama_pos_width = pospeak_width;
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).amp = spkamp;
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).baseline = basesub;
    
    % Gidon Waveform Analysis
    
    last_ind = find(uselead(8:end) < basesub, 1, 'first') + 7;
    
    if isempty(last_ind) == 1
        last_ind = 24;
    end

    wpos = uselead(1:last_ind);
    
    xpos = [0.1:0.1:(0.1 * last_ind)];
    
    gain_init = max(wpos);
    center_init = xpos(find(wpos == max(wpos)));
    stdev_init = 0.3;
    dc_init = wpos(1);
    
    params_init = [gain_init, center_init, stdev_init, dc_init];
    
    [params_fit.pos, ~, ~, pmse] = RunGaussianFit(xpos, wpos', params_init);
    
    
    y_fit.pos = Gaussian(xpos, params_fit.pos(1), params_fit.pos(2), params_fit.pos(3), params_fit.pos(4));
    
    [pactual_re] = Rescale(wpos', [0 1]);
    [pfit_re] = Rescale(y_fit.pos, [0 1]);
    
    p_error = pfit_re - pactual_re;
    
    pos_rescale_mse = mean(p_error.^2);

    %% fit the negative component of the waveform
    
    wneg = uselead((last_ind - 1):end);
    xneg = (0.1 * (last_ind - 1)):0.1:3.2;
    %x.neg = [0.1:0.1:(length(w.neg) * 0.1)];
    
    gain_init = min(wneg);
    center_init = xneg(find(wneg == min(wneg)));
    stdev_init = 0.3;
    dc_init = wneg(end);
    
    params_init = [gain_init, center_init, stdev_init, dc_init];
    [params_fit.neg, ~, ~, nmse] = RunGaussianFit(xneg, wneg', params_init);
    
    y_fit.neg = Gaussian(xneg, params_fit.neg(1), params_fit.neg(2), params_fit.neg(3), params_fit.neg(4));
    
    [nactual_re] = Rescale(wneg', [0 1]);
    [nfit_re] = Rescale(y_fit.neg, [0 1]);
    
    n_error = nfit_re - nactual_re;
    
    neg_rescale_mse = mean(n_error.^2);
    
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).gauss_fit_neg_width = (round(params_fit.neg(3) * 100) / 100);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).gauss_fit_pos_width = (round(params_fit.pos(3) * 100) / 100);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).neg_mse = nmse;
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).ngain = params_fit.neg(1);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).ncenter = params_fit.neg(2);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).nstandev = params_fit.neg(3);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).ndc = params_fit.neg(4);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).nrmse = neg_rescale_mse;
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).nbase_cen = params_fit.neg(4) - params_fit.neg(2);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).pos_mse = pmse;
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).pgain = params_fit.pos(1);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).pcenter = params_fit.pos(2);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).pstandev = params_fit.pos(3);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).pdc = params_fit.pos(4);
    WaveFitParams.(strcat('L',num2str(disindex(wS)))).prmse = pos_rescale_mse;
    
end





