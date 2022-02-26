% Task 3
% Simulation of a M-PAM system with controlled ISI using a duobinary 
% signal X(f) and the Viterbi algorithm in the receiver.
% We assume that |Gr(f)| = |Gt(f)| = sqrt(X(f)).
% We also use the following symbolism:
%  - k: number of bits in a symbol
%  - d: data that we transmit (uniform distribution)
%  - a: M-PAM data that we transmit
%  - n: white Gaussian noise with mean value 0 and variance N0/2
%  - y: data that we receive
%  - b: data that we receive without noise
%  - a_out: data that we receive after detection
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
% Specify the variance of the Gaussian noise
N0 = 0.1;           % Change this value for a different variance
fprintf("Simulation of a %d-PAM system with controlled ISI using a " + ...
    "sequence with %d symbols and Gaussian noise with variance " + ...
    "N0/2=%.2f.\n", M, symbols,N0/2);

% Generate random data from the uniform distribution
d  = unidrnd(M,symbols,1) - 1;

% Convert the data to M-PAM symbols
% Create an array to save all the possible values for the PAM symbols
pam = -(M-1):2:(M-1);
% Set the a0 for t=0 to a0 = -(M-1) (the minimum possible value)
% Find the PAM symbols
a = zeros(symbols+1,1);
a(1) = -(M-1);
a(2:symbols+1) = pam(d+1)';

% Generate Additive White Gaussian Noise
% n is noise after the Gr(f) filter
% The mean is 0 and the varince is 2*N0/pi
mu = 0;
sigma = sqrt(2*N0/pi);
n = normrnd(mu,sigma,symbols,1);

% Calculate the received value, y(m) = a(m) + a(m-1) + n(m)
y = zeros(symbols,1);
for i = 2:symbols+1
    y(i-1) = a(i) + a(i-1) + n(i-1);
end

% Detect the received symbols using the Viterbi algorithm
a_out = detect_viterbi(y,M,symbols);

% Calculate the symbol error rate
n_of_errors = sum(a_out-a'~=0);
ser = n_of_errors/symbols*100;
fprintf("The probability of a symbol being incorrectly detected " + ...
    "is %.6f%%.\n", ser);

% Calculate the symbol error rate for different N0 values
% and plot the results
N0 = [0.01 0.03 0.05 0.08 0.1 0.15 0.2 0.3 0.5 0.8 1 1.5 2];
sigma = sqrt(2*N0/pi);
lsigma = length(sigma);             % helper variable
n = zeros(symbols,lsigma);
y = zeros(symbols,lsigma);
a_out = zeros(symbols+1,lsigma);
ser = zeros(1,lsigma);
fprintf("Error Rate for different values of the noise standard " + ...
    "deviation s=sqrt(2N0/pi):\n");
for i = 1:lsigma
    % Generate Additive White Gaussian Noise
    n(:,i) = normrnd(mu,sigma(i),symbols,1);
    % Calculate the received value
    for j = 2:symbols+1
        y(j-1,i) = a(j) + a(j-1) + n(j-1,i);
    end
    % Detect the received symbols using the Viterbi algorithm
    tmp = detect_viterbi(y(:,i),M,symbols);
    a_out(:,i) = tmp';
    % Calculate the symbol error rate
    n_of_errors = sum(a_out(:,i)-a~=0);
    ser(1,i) = n_of_errors/symbols*100;
    fprintf("For s = %.3f: the probability of a symbol being " + ...
            "incorrectly detected is %.6f%%.\n", sigma(i), ser(1,i));
end
% Plot the results
figure();
plot(N0/2, ser);
xlabel('Noise Variance (N0/2)');
ylabel('Symbol Error Rate (%)');
title('Symbol Error Rate for different values of the AWGN variance');

function s = detect_viterbi(y,M,symbols)
    % DETECT_VITERBI Local function that returns the detected PAM symbols 
    % using the Viterbi algorithm
    % Inputs:
    %   y: received data (with noise)
    %   M: level of the PAM modulation
    %   symbols: number of symbols in the transmitted sequence
    % Output:
    %   s: the detected data

    % All the possible values of the PAM symbols
    S = -(M-1):2:(M-1);
    % Calculate the received value (without noise) for all the possible
    % combinations of the PAM symbols
    b = repmat(S,M,1)+repmat(S,M,1)';
    % Set the μ0 value for t=0
    m = Inf(1,M);
    m(1) = 0;               % a0 = -(M-1)
    % Create an array to save the surviving paths and the detected path
    surv_paths = zeros(M,symbols);
    path = zeros(symbols+1,1);
    % Loop over every situation for t=1 to t=symbols
    for i = 1:symbols
        tmp = m'+(b-y(i)).^2;
        [m,p] = min(tmp);
        surv_paths(:,i) = p';
    end
    [~,path(symbols+1)] = min(m);
    for i = symbols:-1:1
        path(i) = surv_paths(path(i+1),i);
    end
    s = S(path);
end