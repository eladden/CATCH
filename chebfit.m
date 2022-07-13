function [a] = chebfit(g)
%% [a] = chebfit(g)
%fit a chebichev polynomial to the values of g given. The resulting
%coeeficients are given in a.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Elad Denenberg
% Version: 1.0
% Date: 20.12.2019
% License:
% CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This is Matlab, we do the calculation using VECTORS
N = length(g)-1;
I = 2/N*cos((pi()/N)*[0:N].*[0:N]');
I(1,1) = I(1,1)/2;
I(1,N+1) = I(1,N+1)/2;
I(N+1,1) =I(N+1,1)/2;
I(N+1,N+1) =I(N+1,N+1)/2;
I(1,:) = I(1,:)./2;
I(2:N,1) = I(2:N,1)./2;
I(N+1,:) =I(N+1,:)./2;
I(2:N,N+1) =I(2:N,N+1)./2;

a = I*g;

%% The loop way of doing it - should be much slower. Here for documentation
% for j = 0:N 
%     pj = pjcalc(j,N);
%     Ij = zeros(1,N+1);
%     for k = 0:N
%         pk = pjcalc(k,N);
%         
%         Ij(k+1) = 2/(pj*pk*(N))*cos(j*pi()*k/(N));
%     end
%     a(j+1) = sum(Ij.*g'); %we use j+1 here because the loop starts with 0 like this was C++ 
% end
% 
% function pjout = pjcalc(j,N)
% 
% if (j==0) || (j==N)
%     pjout = 2;
% else 
%     pjout = 1;
% end
        
