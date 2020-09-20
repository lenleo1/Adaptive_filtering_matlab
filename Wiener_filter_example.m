%% Wiener Filter Example
t = (0:0.005:10)';
fs = 1/0.005;
dn = 10*sin(2*pi*5*t)+15*sin(2*pi*6*t);  % expected signal
xn = dn + 5*randn(length(t),1);          % Observed signal

M  = floor(length(xn)*0.3);              % order of filter
itr = floor(length(xn)*1);               % Break point coordinates
rho_max = max(eig(xn(1:M)*xn(1:M)'));    % Maximum eigenvalue of input signal correlation matrix
miu = 1/sqrt(2)*(1/rho_max);             % 0 < mu < 1/rho

[yn,Wopt,en] = classic_lms(xn,dn,M,miu,itr);

subplot(411);plot(xn);          ylabel('Observed signal X');axis tight;
subplot(412);plot(dn);          ylabel('expected signal d');axis tight;
subplot(413);plot(yn);          ylabel('recovered signal y');axis tight;
subplot(414);plot(abs(en),'r'); ylabel('error signal e');axis tight;
%% function
function [yn,Wopt,en] = classic_lms(xn,dn,M,miu,itr)
% -------------------------------------------------
% function [yn,Wopt,en] = classiclms(xn,dn,M,mu,itr)
% use the method of steepest descent
% Inputs:
%     xn   Observation signal sequence   (column vector)
%     dn   expected signal sequence      (column vector)
%     M    filter order (>1)             (scalar)
%     miu  convergence factor(step size) (scalar)
%          greater than 0 and less than the reciprocal 
%          of the maximum eigenvalue of the correlation matrix   
%     itr  Iterative end point index,    (scalar)      
%          Defaults to the length of xn,
%          M < itr <= length(xn)
% Outputs:
%    Wopt  filter weight vector (column vector)
%     en   error sequence (N x 1) (column vector)  
%     yn   actual output sequence (column vector)
% ---------------------------------------------------
N = length(dn);                 % total number of samples
                                % The number of parameters must be 4 or 5
if nargin == 4                  % 4, the number of recursive iterations defaults to the length of xn
    itr = length(xn);
elseif nargin == 5              % 5,M < itr =< length(xn)
    if itr>N || itr<M
        error('The number of iterations is too large or too small!');
    end
else
    error('Please check the number of input parameters!');
end

if M < 2
    error('The filter order must be greater than one!');
end

% Initialization parameters
en = NaN*zeros(N,1);    % error sequence, en(k) represents the error between the expected output d and the actual output y at the kth iteration
W  = zeros(M,1);        % weighting parameter, the initial elements are 0                           

% Iteration
for k = M:itr                  % k-th iteration
    x = xn(k:-1:k-(M-1));      % filter input with M taps
    y = W'* x;                 % k-th output
    en(k) = dn(k) - y ;        % k-th error                            
    W = W + 2*miu * en(k) * x; % weight updating 
end
Wopt = W;                      % Optimal weight

% The output sequence of the optimal filter
yn = NaN * ones(size(xn));
for k = M:N
    x = xn(k:-1:k-M+1);
    yn(k) = W'* x;
end
end
