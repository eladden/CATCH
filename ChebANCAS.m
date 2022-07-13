function [tmin,Dmin,Sampled,Built] = ChebANCAS (satrec1_0,satrec2_0,dt,max_test_time,consts,min_t,N,delta)
%% [tmin,Dmin,multimin] = ChebANCAS (satrec1_0,satrec2_0,dt,max_test_time,consts,min_t,N,delta)
% use Chebyshev Polynomials in the same way Alfano and Negron Conjuction 
% Analysis Software uses polynomials to find the times of closest approach
% The algorithm and the useful values for the parameters min_t N and delta
% are given in:
%   Elad Denenberg, "Satellite closest approach calculation through Chebyshev
%   Proxy Polynomials," Acta Astronautica, Volume 170, 2020, Pages 55-65, ISSN 0094-5765, 
%   https://doi.org/10.1016/j.actaastro.2020.01.020 
%
% inputs:
%    satrec1_0, satrec2_0 - satellites in the SGP4 implementation format
%    dt                   - the time difference between the epoch of the
%                           first and second satellites
%    max_test_time        - The horizon of the test (Using SGP4 longer than
%                           two weeks is discouraged
%    consts               - Model Constants
%    min_t                - Minimal time difference between samples -
%                           required accuracy
%    N                    - Maximal order
%    delta                - The number used to devide the minimal orbital
%                           period
%In the paper it was found that good values are N=16 and delta=2, But this
%can be fine tuned for better accuray at the price of longer computation
%time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Elad Denenberg
% Version: 1.0
% Date: 20.12.2019
% License:
% CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T1 = 2*pi()/satrec1_0.no;
T2 = 2*pi()/satrec2_0.no;
if nargin < 6
    min_t = 0;
    delta = 2;
    N = 16;
end

T = min([T1 T2])/delta; %This is \Gamma in the paper

p = min_t+T/2:T:max_test_time;
p = repmat(p',[1,N+1]);
k = 0:1:N;
xk = repelem(T/2*cos(pi()*k/N),size(p,1),1);

p = p + xk; %this now contains all the timepoints to be sampled 


%% evaluate each point

%This is the epoch time difference between the two objects
if dt > 0
    dts = 0;
    dtd = dt;
elseif dt < 0
    dts = -dt;
    dtd = 0;
else
    dts = 0; dtd =0;
end

%% evaluate
%slice p for the parfor loop
a = zeros(size(p));
bx = zeros(size(p));
by = zeros(size(p));
bz = zeros(size(p));
fp_dot = zeros(size(p,2),1);
rdx = zeros(size(p,2),1);
rdy = zeros(size(p,2),1);
rdz = zeros(size(p,2),1);

for i = 1:size(p,1)
    for j = 1:size(p,2)
          if ~((i > 1) && (j ==size(p,2))) 
            [~, r1_teme,v1_teme] = sgp4(satrec1_0,p(i,j)+dts,consts);
            [~, r2_teme,v2_teme] = sgp4(satrec2_0,p(i,j)+dtd,consts);
         
        
  
            rd      = r1_teme - r2_teme;
            rd_dot  = v1_teme - v2_teme;
        
            rdx(j) = rd(1);
            rdy(j) = rd(2);
            rdz(j) = rd(3);
            fp_dot(j)  = 2*dot(rd_dot,rd);
    
         else
            rdx(j) = rdxprev;
            rdy(j) = rdyprev;
            rdz(j) = rdzprev;
            fp_dot(j)  = fpprev;             
         end
    end
    rdxprev =  rdx(1);
    rdyprev = rdy(1);
    rdzprev = rdz(1);
    fpprev = fp_dot(1);
    bx(i,:) = chebfit(rdx); %fit polinomial curves
    by(i,:) = chebfit(rdy);
    bz(i,:) = chebfit(rdz);
    a(i,:) = chebfit(fp_dot);
end
Sampled = size(p,1)*size(p,2);
Built = size(p,1); %for profiling

Dminv = ones(size(p,1),1)*inf;
tminv = ones(size(p,1),1)*(-1);
%% spline and find roots
for i = 1:size(p,1)
    aloop = a(i,:);
    Nloop = N;
    while (aloop(end) ==0)
        aloop = aloop(1:end-1);
        Nloop = Nloop -1;
    end
    A = zeros(Nloop,Nloop);
    A(1,2) = 1; % first row
    
    A(Nloop,:) = -aloop(1:Nloop)/(2*aloop(Nloop+1)); % last row
    A(Nloop, Nloop - 1) = A(Nloop , Nloop - 1) + 0.5;
    for j = 2:Nloop-1 %every other row
        A(j,j-1) = 0.5;
        A(j,j+1) = 0.5;
    end
    
    tcand = eig(A);
    tcand = tcand(imag(tcand) == 0); %this is between 1 and -1
    tcand = tcand(tcand < 1 & tcand > -1); %only solutions in the interval
    if ~isempty(tcand)
        rx = evaluate_cheb_pol(tcand,bx(i,:));
        ry = evaluate_cheb_pol(tcand,by(i,:));
        rz = evaluate_cheb_pol(tcand,bz(i,:));
        
        Dcand = sqrt(rx.^2 + ry.^2 + rz.^2);
        tcandReal = tcand*T/2 + (p(i,end)+p(i,1))/2;
       
        
        [Dcand,idx] = min(Dcand);
    
        %if (Dmin > Dcand)
           Dminv(i) = Dcand;
           tminv(i) = tcandReal(idx);
        %end
    end
end
[Dmin,indx] = min(Dminv);
tmin = tminv(indx);
