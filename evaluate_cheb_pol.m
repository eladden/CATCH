function y = evaluate_cheb_pol(x,a,xmax,xmin)
%% function y = evaluate_cheb_pol(x,xmax,xmin,a)
% evaluates the function with the Chebyshev polynomail

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Elad Denenberg
% Version: 1.0
% Date: 20.12.2019
% License:
% CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = 0;
if nargin >2 && nargin ==4
    for i = 1:length(a)
        y = y + a(i)*Tcheb(x,i-1,xmax,xmin); 
    end
elseif nargin == 2
    for i = 1:length(a)
        y = y + a(i)*Tcheb(x,i-1); 
    end
else
    error("wrong input");
end