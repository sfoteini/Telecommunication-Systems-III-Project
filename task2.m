% Task 2
% Simulation of a M-PAM system with controlled ISI using a duobinary 
% signal X(f) and precoding in the transmitter.
% We assume that |Gr(f)| = |Gt(f)| = sqrt(X(f)).
% We also use the following symbolism:
%  - k: number of bits in a symbol
%  - d: data that we transmit (uniform distribution)
%  - p: precoded data
%  - a: M-PAM data that we transmit
%  - n: white Gaussian noise with mean value 0 and variance N0/2
%  - y: data that we receive
%  - b: data that we receive without noise
%  - d_out: data that we receive after detection
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

% Precoding: Pm = (Dm - Pm-1)(modM), P0 = 0
p = zeros(symbols+1,1);
for i=2:symbols+1
    p(i) = mod(d(i-1)-p(i-1),M);
end

% Convert the precoded data to M-PAM symbols
% am = (2pm - (M-1))*d, where we d = 1
a = 2*p - (M-1);

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

% Detect the received symbol, d = (b/(2d)+M-1)(modM),
% where d = 1 and b is the received symbol (without noise)
d_out = detect(y,M);

% Display the received y symbols, the ideal b symbols (without noise)
% and the decision boundaries on a scatter plot
% i = p(m)+p(m-1), i = {0,1,...,2(M-1)}
i = 0:2*(M-1);
% b = 2*(p(m)+p(m-1)-(M-1)) = 2(i-(M-1)), b = {-2(M-1),...,0,...,2(M-1)}
b = -2*(M-1):2:2*(M-1);
% d = (b/2+M-1)modM = i mod M = i-M, if M <= k <= 2(M-1)
%                               i  , if 0 <= k <= M-1
db = zeros(1,length(i));
db(1,1:M) = i(1,1:M);
db(1,M+1:end) = i(1,M+1:end) - M;

scatter(y,zeros(length(y),1),10,'black');
hold on;
scatter(b,zeros(length(b),1),20,'red','filled','s');
% Display the decision boundaries
xline(b(1:end-1)+1,'--b');
% Display the detected symbol for each decision area
txt = strings(1,length(db));
for j = 1:length(db)
    txt(j) = ['d=' num2str(db(j))];
end
text(b-M/10,0.2*ones(1,length(b)),txt,'FontSize',8);
legend('Received symbols','Ideal symbols','Decision boundaries');
title('Received constellation');

% Calculate the symbol error rate
n_of_errors = sum(d-d_out ~= 0);
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
d_out = zeros(symbols,lsigma);
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
    % Detect the received symbol
    d_out(:,i) = detect(y(:,i),M);
    % Calculate the symbol error rate
    n_of_errors = sum(d-d_out(:,i) ~= 0);
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



function d = detect(y,M)
    % DETECT Local function that returns the detected data based on the
    % value of the received symbol and the level of the PAM modulation
    % Inputs:
    %   y: received data (with noise)
    %   M: level of the PAM modulation
    % Output:
    %   s: the detected data
    
    % We want to find the even number that is nearest to the received value
    % First, we calculate the floor of our value
    b = floor(y);
    % If s is even, then mod(y,2) < 1.
    % If s is odd, then mod(y,2) > 1. Hence we have to add 1 to the 
    % number s to find the nearest even number.
    i = mod(y, 2) > 1;
    b(i) = b(i) + 1;
    
    % If the |b| > 2*(M-1) (the max possible value), then
    % set |b| = 2*(M-1)
    if b > 2*(M-1)
        b = 2*(M-1);
    elseif b < -2*(M-1)
        b = -2*(M-1);
    end

    % Calculate the received data using the quantized value b
    d = mod(b/2+M-1,M);
end