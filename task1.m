% Task 1
% Simulation of a M-PAM system with zero ISI using a Raised-cosine 
% filter Xrc(f).
% We assume that |Gr(f)| = |Gt(f)| = sqrt(Xrc(f)).
% We also use the following symbolism:
%  - k: number of bits in a symbol
%  - d: data that we transmit (uniform distribution)
%  - a: M-PAM data that we transmit
%  - n: white Gaussian noise with mean value 0 and variance N0/2
%  - y: data that we receive, y = a + n
%  - a_out: M-PAM data after demodulation/detection
%
% Author: Foteini Savvidou

% Clear the workspace
clc;
clear;
close all;

% Specify the levels of the PAM signal
k = 2;              % Change this value for a different PAM modulation
M = 2^k;
% Specify the number of symbols in a sequence
symbols = 1e4;      % Change this value for a different number of 
                    % symbols in a sequence
% Specify the variance of the Gaussian noise (var=N0/2)
N0 = 0.1;           % Change this value for a different variance
fprintf("Simulation of a %d-PAM system with zero ISI using a " + ...
    "sequence with %d symbols and Gaussian noise with variance " + ...
    "%.2f.\n", M, symbols,N0/2);

% Generate random data from the uniform distribution
d  = unidrnd(M,symbols,1) - 1;

% Convert the data to M-PAM symbols
% Create an array to save all the possible values for the PAM symbols
pam = -(M-1):2:(M-1);
% Find the PAM symbols
a = pam(d+1)';

% Generate Additive White Gaussian Noise
% n is noise after the Gr(f) filter
% The mean is 0 and the varince is N0/2
mu = 0;
sigma = sqrt(N0/2);
n = normrnd(mu,sigma,symbols,1);

% Calculate the received value
y = a + n;

% Detect the received symbol
a_out = detect(y,M);

% Display the received symbols, the ideal M-PAM symbols
% and the decision boundaries on a scatter plot
scatter(y,zeros(length(y),1),10,'black');
hold on;
scatter(pam,zeros(length(pam),1),20,'red','filled','s');
xline(pam(1:end-1)+1,'--b');
legend('Received symbols','PAM symbols','Decision boundaries');
title('Received constellation');

% Calculate the symbol error rate
n_of_errors = sum(a-a_out ~= 0);
ser = n_of_errors/symbols*100;
fprintf("The probability of a symbol being incorrectly detected " + ...
    "is %.6f%%.\n", ser);

% Calculate the baud rate for different values of the roll-of factor
rollof_factor = 0:0.1:1;
R = 2./(1+rollof_factor);
figure();
plot(rollof_factor, R);
xlabel('Roll-of factor');
ylabel('Baud rate (W symbols/sec)');
title('Baud rate for different values of the roll-of factor');

% Calculate the symbol error rate for different N0 values
% and plot the results
N0 = [0.01 0.03 0.05 0.08 0.1 0.15 0.2 0.3 0.5 0.8 1 1.5 2];
sigma = sqrt(N0/2);
lsigma = length(sigma);             % helper variable
n = zeros(symbols,lsigma);
y = zeros(symbols,lsigma);
a_out = zeros(symbols,lsigma);
ser = zeros(1,lsigma);
fprintf("Error Rate for different values of the noise standard " + ...
    "deviation s=sqrt(N0/2):\n");
for i = 1:lsigma
    % Generate Additive White Gaussian Noise
    n(:,i) = normrnd(mu,sigma(i),symbols,1);
    % Calculate the received value
    y(:,i) = a + n(:,i);
    % Detect the received symbol
    a_out(:,i) = detect(y(:,i),M);
    % Calculate the symbol error rate
    n_of_errors = sum(a-a_out(:,i) ~= 0);
    ser(1,i) = n_of_errors/symbols*100;
    fprintf("For s = %.3f: the probability of a symbol being " + ...
            "incorrectly detected is %.6f%%.\n", sigma(i), ser(1,i));
end
% Plot the results
figure();
plot(N0/2, ser);
xlabel('Noise Variance');
ylabel('Symbol Error Rate (%)');
title('Symbol Error Rate for different values of the AWGN variance');



function s = detect(y,M)
    % DETECT Local function that returns the detected PAM symbol based on the 
    % value of the received symbol and the level of the PAM modulation
    % Inputs:
    %   y: received data (with noise)
    %   M: level of the PAM modulation
    % Output:
    %   s: the detected PAM symbol

    % We want to find the odd number that is nearest to the received value
    % First, we calculate the floor of our value
    s = floor(y);
    % If s is even, then mod(y,2) < 1. Hence we have to add 1 to the 
    % number s to find the nearest PAM symbol
    % If s is odd, then mod(y,2) > 1 and we have found our symbol
    i = mod(y, 2) < 1;
    s(i) = s(i) + 1;
    
    % If the |s| > M-1 (the max possible value of a M-PAM symbol), then
    % set |s| = M-1
    if s > M-1
        s = M-1;
    elseif s < -(M-1)
        s = -(M-1);
    end
end