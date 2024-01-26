% Initialize variables
sampling_rate = 1000; % 1 kHz
duration = 1; % 1 second
t = 0:1/sampling_rate:duration-1/sampling_rate; % time vector
f = 5; % 5 Hz
phase_shifts = 0:pi/8:2*pi; % phase shifts from 0 to 2*pi
dtw_distances = zeros(size(phase_shifts)); % initialize DTW distances
pearsons = zeros(size(phase_shifts)); 

% Generate reference sine wave
reference_sine = sin(2*pi*f*t);

% Loop over each phase shift and compute DTW distance
for i = 1:length(phase_shifts)
    phase_shift = phase_shifts(i);
    shifted_sine = sin(2*pi*f*t + phase_shift);
    dtw_distances(i) = dtw(reference_sine, shifted_sine);
    pearsons(i) = 1-fastcorr(reference_sine, shifted_sine);
end

% Plot DTW distances as a function of phase shift
figure;
yyaxis left
plot(phase_shifts, dtw_distances, '-o');
ylabel('DTW Distance');
hold on
yyaxis right
plot(phase_shifts, pearsons, '-o');
xlabel('Phase Shift (radians)');
ylabel('r Distance');
title('DTW Distance vs Phase Shift');
grid on;

%% Frequency modulation
% Initialize variables
sampling_rate = 1000; % 1 kHz
duration = 1; % 1 second
t = 0:1/sampling_rate:duration-1/sampling_rate; % time vector
f = 5; % 5 Hz
freq_shifts = 0:0.1:5; % frequency shifts from 0 to 1 Hz
dtw_distances = zeros(size(freq_shifts)); % initialize DTW distances
pearsons = zeros(size(freq_shifts));

% Generate reference sine wave
reference_sine = sin(2*pi*f*t);

% Loop over each frequency shift and compute DTW distance
for i = 1:length(freq_shifts)
    freq_shift = freq_shifts(i);
    shifted_freq_sine = sin(2*pi*(f + freq_shift)*t);
    dtw_distances(i) = dtw(reference_sine, shifted_freq_sine);
    pearsons(i) = 1-fastcorr(reference_sine, shifted_freq_sine);
end

% Plot DTW distances as a function of frequency shift
figure;
yyaxis left
plot(freq_shifts, dtw_distances, '-o');
ylabel('DTW Distance');
hold on
yyaxis right
plot(freq_shifts, pearsons, '-o');
xlabel('Frequency Shift (Hz)');
ylabel('r Distance');
title('DTW Distance vs Frequency Shift');
grid on;

%% Amplitude modulation
% Initialize variables
sampling_rate = 1000; % 1 kHz
duration = 1; % 1 second
t = 0:1/sampling_rate:duration-1/sampling_rate; % time vector
f = 5; % 5 Hz
freq_shifts = 0:0.1:2; % frequency shifts from 0 to 1 Hz
dtw_distances = zeros(size(freq_shifts)); % initialize DTW distances
pearsons = zeros(size(freq_shifts));

% Generate reference sine wave
reference_sine = sin(2*pi*f*t);

% Loop over each frequency shift and compute DTW distance
for i = 1:length(freq_shifts)
    freq_shift = freq_shifts(i);
    shifted_freq_sine = sin(2*pi*(f)*t)*freq_shift;
    dtw_distances(i) = dtw(reference_sine, shifted_freq_sine);
    pearsons(i) = 1-fastcorr(reference_sine, shifted_freq_sine);
end

% Plot DTW distances as a function of frequency shift
figure;
yyaxis left
plot(freq_shifts, dtw_distances, '-o');
ylabel('DTW Distance');
hold on
yyaxis right
plot(freq_shifts, pearsons, '-o');
xlabel('Amplitude multiplier');
ylabel('r Distance');
if nanmean(pearsons)<(1.0e-5); ylim([0 2]); end
title('DTW Distance vs Frequency Shift');
grid on;
