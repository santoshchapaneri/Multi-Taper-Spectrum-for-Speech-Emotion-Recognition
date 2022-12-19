function [thomson,s]=thomson(N,now);

% [thomson,s]=thomson(N,now);

% Calculates the 'now' (default 8) windows 
% 
% now : A number between 2 and 8 (default 8)
% N   : Window length
% thomson: The windows in a matrix as columns (with weights). 
%
% The program creates a file and a variable 
% called 'thomson.mat' and 'thomson'.
%
% The method is based on the following reference:
%
%   D.J. Thomson, "Spectrum Estimation and Harmonic Analysis", Proc. IEEE,
%   70(9), pp. 1055--1096, September 1982.
%
% (C) Copyright 2010 
% Mathematical Statistics, Centre for Mathematical Sciences, Lund University, Lund, Sweden
% School of Computing, University of Eastern Finland, Joensuu, Finland

if nargin==1
  now=8;
end

file='thomson';
loadfile=[file int2str(N) '_' int2str(now) '.mat'];

c=exist(loadfile);

if c==2 
  eval(['load ' loadfile])
  if length(thomson(:,1))~=N | length(thomson(1,:))~=now
    c=0;
  end
end

if c==0 | c==4 | c==3 | c==5

  B=(now+2)/N ;   % Resolution in spectrum

  l=[1:N-1]';

  r=2*sin(pi*B.*l)./(2*pi.*l);
  rbox=[B;r];

  Ry=toeplitz(rbox);
  [FN,s,v]=svd(Ry); 

  thomson=FN(:,1:now);
  
  s=ones(now,1)/now;

  eval(['save ' loadfile ' ' 'thomson s'])
end