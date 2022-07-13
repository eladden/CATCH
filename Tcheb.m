function [T] = Tcheb(x,j,xmax,xmin)
%% [T] = Tcheb(k,N)
%calculate the monomial of the chebychev polynomail

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Elad Denenberg
% Version: 1.0
% Date: 20.12.2019
% License:
% CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 2 && nargin == 4
    x = (2*x  - (xmax+xmin))/(xmax-xmin);
elseif nargin ~= 2
    error("wrong input");
end
T = cos(j*acos(x));