function [multipeak,s] = multipeakwind(N,now);

% [multipeak,s]=multipeakwind(N,now);
% Calculates the 'now' (default 8) windows
%
% now : A number between 2 and 8 (default 8)
% N   : Window length
% multipeak: The windows in a matrix as columns (with weights).
%
% The program creates a file and a variable
% called 'multipeak.mat' and 'multipeak'.
%
% The method is based on the following reference:
%
%   M. Hansson and G. Salomonsson, "A Multiple Window Method for Estimation of Peaked Spectra",
%   IEEE T. on Sign. Proc. 45(3), pp. 778--781, March 1997.

if nargin==1
    now=8;
end

file='multipeak';
loadfile=[file int2str(N) '_' int2str(now) '.mat'];

c=exist(loadfile);

if c==2
    eval(['load ' loadfile])
    if length(multipeak(:,1))~=N | length(multipeak(1,:))~=now
        c=0;
    end
end

if c==0 | c==4 | c==3 | c==5

    K1=20;         % Peak in dB
    K2=30;         % Penalty value in dB
    
    B=(now+2)/N ;   % Resolution in spectrum

    loge=log10(exp(1));
    C=2*K1/10/B/loge;
    l=[1:N-1]';
    r0=2/C*(1-exp(-C*B/2));
    r=(2*C-exp(-C*B/2)*(2*C*cos(pi*B.*l)-4*pi.*l.*sin(pi*B.*l)))./(C^2+(2*pi.*l).^2);

    rpeak=[r0;r];  % Covariance function peaked spectrum

    r=2*sin(pi*B.*l)./(2*pi.*l);
    rbox=[B;r];

    rpen=10^(K2/10)*[1;zeros(N-1,1)]-(10^(K2/10)-1)*rbox; % Covariance function penalty function

    Ry=toeplitz(rpeak);
    Rx=toeplitz(rpen);
    [RR]=chol(Rx);
    C=inv(RR')*Ry*inv(RR);
    [Q,T]=schur(C);
    F=inv(RR)*Q;
    RD=F'*Ry*F;
    RD=diag(RD);
    [RDN,h]=sort(RD);
    for i=1:length(RD)
        FN(:,i)=F(:,h(i));
        FN(:,i)=FN(:,i)/sqrt(FN(:,i)'*FN(:,i));
    end
    RDN=RDN(length(RD):-1:1);
    FN=FN(:,length(RD):-1:1);

    s=RDN(1:now)/sum(RDN(1:now))  ;
    multipeak=FN(:,1:now);
    eval(['save ' loadfile ' ' 'multipeak s'])
end