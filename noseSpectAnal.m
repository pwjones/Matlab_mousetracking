function [f, Y] = noseSpectAnal(traj, dt)

min_df = .05;
n = numel(traj);
traj = traj(:); %make sure it's a column vector

T = dt;                     % Sample time (s)
Fs = 1/T;                % Sampling frequency (Hz)                   
L = numel(traj);               % Length of signal
t = (0:L-1)*T;                % Time vector
y = traj;
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
figure;
subplot(2,1, 1);
plot(t(1:1000),y(1:1000))
title('Nose position (1000 frames)')
xlabel('time (seconds)')

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
df = mean(diff(f));
samplef = floor(min_df/df);
if samplef > 1
    f = downsample(f, samplef)
    Y = downsample_summing(Y, samplef);
end
    
% Plot single-sided amplitude spectrum.
subplot(2,1,2);
%plot(f,2*abs(Y(1:NFFT/2+1)))
plot(f,2*abs(Y(1:numel(f))))
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')


function ret_vect = downsample_summing(vect, n)

ret_vect = zeros(floor(numel(vect)/n), 1);
for ii = 1:length(ret_vect)
    range = ((1:n)-1)*ii + 1;
    ret_vect(ii) = sum(vect(range));
end
