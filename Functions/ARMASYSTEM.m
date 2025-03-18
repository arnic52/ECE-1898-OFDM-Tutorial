function y = ARMASYSTEM(x, b, a)
% ARSYSTEM Implement ARMA system
%   y = ARSYSTEM(x) implements the autoregressive moving average system
%   y[n] = b_0 x[n] + b_1 x[n-1] + ... + b_M x[n-M]
%          - a_1 y[n-1] - ... - a_N x[n-N]

    N = length(a);                             % order of AR system
    M = length(b)-1;                           % order of MA system
    nmax = length(x)-1;                        % max time for x[n] & y[n]
    x = [zeros(M,1); x];                        % set x[n] = 0 for n < 0
    y = [zeros(N,1); zeros(nmax+1,1)];          % init. y to all zeros
    
    for n = 0 : nmax                           % loop over 0 <= n <= nmax
        iX = n + M + 1;                        % index for x array
        iY = n + N + 1;                        % index for y array
        
        xPast = x(iX : -1 : iX-M);                              % pres/past inputs, reversed
        yPast = y(iY-1 : -1 : iY-N+1);                          % past outputs, reversed
        y(iY) = (sum(b.*xPast) - sum(a(2:end).*yPast))/a(1);    % compute y[n] by diff eqn
    end

    y = y(N+1+extra : end);                          % remove y[-N] ... y[-1]

%     Nm1 = length(a);                        % number of a_k coeffs
%     M = length(b);                          % number of b_m coeffs
%     nmax = length(x)-1;                     % max time for x[n] & y[n]
% 
%     x = [zeros(M-1) x];                     % set x[n] = 0 for n < 0
%     y = [zeros(1,Nm1) zeros(1,nmax+1)];     % init. y to all zeros
%     
%     for n = 0 : nmax                        % loop over 0 <= n <= nmax
%         indexX = n + M;                     % index of x array
%         indexY = n + Nm1 + 1;               % index for y array
% 
%         % past inputs & outputs in reverse order
%         xPast = x(indexX-1 : -1 : indexX-M+1);
%         yPast = y(indexY-1 : -1 : indexY-Nm1);
%         % compute y[n] by difference equation
%         y(indexY) = sum(b .* xPast) - sum(a .* yPast);  
%     end
% 
%     y = y(Nm1+1 : end);                     % remove y[-N+1] ... y[-1]
end
