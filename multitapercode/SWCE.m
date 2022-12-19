function [h,s] = SWCE(N, NumWindows)

% Sinusoidal Windows Cepstrum Estimator (SWCE)
% Different window functions as columns in h (N X NumWindows)
% Different weights in the vector s (NumWindows X 1)
%
% The method is based on the following reference:
%
%   M. Hansson-Sandsten and J. Sandberg, "Optimal cepstrum estimation
%   using multiple windows," in Proc. ICASSP 2009, Taipei, Taiwan, April
%   2009, pp. 3077�3080.

M = fix(N/NumWindows);
for i=1:NumWindows
    h(:,i)=sqrt(2/(N+1))*sin((pi*i*[1:N]')/(N+1));
end
s=((cos(2*pi*[0:fix(N/M)-1]'*M/N/2))+1)./sum(cos(2*pi*[0:fix(N/M)-1]'*M/N/2)+1);