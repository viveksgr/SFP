 function features = extractSniffFeatures(fless_superset, feat_superset,fs)
    ntrials = size(fless_superset,1);
    samsize = size(fless_superset,2);
    
    % Initialize features vector
    features = zeros(ntrials, 14);

    for tt = 1:ntrials
        tin_off = feat_superset(tt,7)-feat_superset(tt,5);
        tin_peak = feat_superset(tt,1)-feat_superset(tt,5);
        tex_on = min(feat_superset(tt,6)-feat_superset(tt,5),samsize);
        tex_peak = min(feat_superset(tt,2)-feat_superset(tt,5),samsize);
        tex_off = min(feat_superset(tt,8)-feat_superset(tt,5),samsize);


        % Vol ratio
        features(tt,1) = feat_superset(tt,11)./feat_superset(tt,12);

        % Duration ratio
        features(tt,2) = feat_superset(tt,13)./feat_superset(tt,14);

        % Flow
        features(tt,3) = feat_superset(tt,3)./feat_superset(tt,4);

        % Times
        features(tt,4) = feat_superset(tt,2)-feat_superset(tt,1);


        features(tt,5) = calculateSmoothness(fless_superset(tt,1:tin_peak), fs);
        features(tt,6) = calculateSmoothness(fless_superset(tt,tin_peak: tex_peak), fs);
        features(tt,7) = calculateSmoothness(fless_superset(tt, tex_on:tex_off ), fs);

        features(tt,8) = calculateCurvature(fless_superset(tt,1:tin_peak), fs);
        features(tt,9) = calculateCurvature(fless_superset(tt,tin_peak: tex_peak), fs);
        features(tt,10) = calculateCurvature(fless_superset(tt, tex_on:tex_off), fs);

        phase1 = fless_superset(tt,1:tin_peak);
        phase2 = fless_superset(tt,tin_peak: tex_peak);
        phase3 = fless_superset(tt,tex_peak:end);

        slope1 = calculatePhase2Slope(phase1, fs);
        slope2 = calculatePhase2Slope(phase2, fs);
        slope3 = calculatePhase2Slope(phase3, fs);

        features(tt,11) =  slope2;
        features(tt,12) = max(diff(phase2)); % Max flow rate change
        features(tt,13) = calculateSmoothnessAroundPoint(phase2, tin_peak, fs);
        features(tt,14) = calculateSmoothnessAroundPoint(phase2,  tex_peak, fs);

        % % Slope ratios
        % features(tt,15) = abs(slope2/slope1);
        % features(tt,16) = abs(slope2/slope3);
        % features(tt,17) = abs(slope1/slope3);

        % Symmetry of trace
        % try
        %     features(tt,18) = calculateSymmetryAroundPoint(phase2, tin_peak);
        %     features(tt,19) = calculateSymmetryAroundPoint(phase2,  tex_peak);
        % catch
        %     save('dataset')
        % end

    end
end

% Define sub-functions for smoothness, curvature, phase 2 slope, etc.
function smoothness = calculateSmoothness(signal, fs)
    % Apply a Fast Fourier Transform
    % Y = fft(signal);
    % L = length(signal);
    % P2 = abs(Y/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);
    % 
    % f = fs*(0:(L/2))/L;
    % highFreqPower = sum(P1(f > fs/4)); % Power in the upper quarter of frequencies
    % totalPower = sum(P1);
    % 
    % % Smoothness can be inversely related to the proportion of high-frequency power
    % smoothness = 1 - (highFreqPower / totalPower);

    smooth_sig = movmean(signal,fs/10);
    smoothness = -mean(abs(smooth_sig-signal));
end

function curvature = calculateCurvature(signal, fs)
    % Calculate the second derivative
    secondDerivative = diff(signal, 2);

    % Curvature as the mean of the absolute second derivative
    curvature = mean(abs(secondDerivative));
    if isnan(curvature); curvature=0;end
end


function slope = calculatePhase2Slope(signal, fs)
    % Linear regression to find the slope
    x = (1:length(signal))';
    p = polyfit(x, signal, 1);
    slope = p(1); % The first coefficient is the slope
end

function smoothness = calculateSmoothnessAroundPoint(signal, point, fs)
    % Define a window size, e.g., 5% of the signal length
    windowSize = round(0.05 * length(signal));
    startIndex = max(1, point - windowSize);
    endIndex = min(length(signal), point + windowSize);

    % Calculate smoothness within this window
    smoothness = calculateSmoothness(signal(startIndex:endIndex), fs);
end


function symmetryScore_n = calculateSymmetryAroundPoint(signal, point)
    % Define a window size, e.g., 5% of the signal length
    windowSize = round(0.05 * length(signal));
    startIndex = max(1, point - windowSize);
    endIndex = min(length(signal), point + windowSize);

    % Calculate smoothness within this window
    if endIndex>startIndex+1
        symmetryScore_n =  computeTimeSymmetry(signal(startIndex:endIndex));
    else
        symmetryScore_n = nan;
    end
end

function symmetryScore = computeTimeSymmetry(signal)
    % Ensure the signal is a row vector
    if iscolumn(signal)
        signal = signal';
    end
    
    % Time-flip the signal
    flippedSignal = fliplr(signal);

    % Handle signals with an odd number of elements by removing the middle element
    if mod(length(signal), 2) == 1
        midIndex = ceil(length(signal) / 2);
        signal(midIndex) = [];
        flippedSignal(midIndex) = [];
    end

    % Compute the Pearson correlation coefficient between the original and flipped signals
    symmetryScore = corr(signal', flippedSignal');
end
