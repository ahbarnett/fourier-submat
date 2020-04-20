function s = sinc(x)
% SINC
%
% fixed up for vector arguments. Barnett 3/30/20
s = sin(x)./x;
s(x==0) = 1;
