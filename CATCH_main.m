%% Constants
fprintf("setting up the problem...");
xpdotp = 1440.0 / (2.0*pi); 
mu = 398600.5;
Re = 6378.137;
J2 = 0.00108262998905;
%j2 = J2;
whichconsts = 72;
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc(whichconsts);
consts = [tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2];

%% Example of how to run CATCH
% SKYSAT-1
s1 = '1 39418U 13066C   19042.20180787 -.00000231  00000-0 -14135-4 0  9996';
s2 = '2 39418  97.6613 126.3667 0023103  57.2848 303.0584 14.98788403285542';
% FENGYUN 1C DEB
d1 = '1 36258U 99025DWK 19042.03235178  .00000157  00000-0  33857-4 0  9997';
d2 = '2 36258  99.0353 139.3348 0205784 231.9458 126.3017 14.62433061544647';

%% load the TLE

[satrecse, ~, ~, ~] = twoline2rv(whichconsts, s1, s2, 'c');
[satrecde, ~, ~, ~] = twoline2rv(whichconsts, d1, d2, 'c');

dt = (satrecse.jdsatepoch-satrecde.jdsatepoch)*24*60; %epoch difference in minutes

%% compute initial state and time
%Uncomment the next line to download the latest data
%websave('finals2000A.txt','https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt'); %load the information about dT1 and dAT
datfileh = fopen('finals2000A.txt');
fseek(datfileh, 0,'eof');
filelength = ftell(datfileh);
fclose(datfileh);
finals2000A = importfile('finals2000A.txt', 1, filelength);
finals = finals2000A;
save finals finals
dAT = 35;

if dt > 0
    dts = 0;
    dtd = dt;
    [year,mon,day,hr,mi,sec] = invjday ( satrecse.jdsatepoch );
    ep = satrecse.jdsatepoch;
    index = find(ismember(finals(:,1:3),[year-2000 mon day],'rows'));
    dUT1 = finals(index,11);
    [~, ~, ~, ~, ~, ~, ttt] =...
        convtime (year,mon,day,hr,mi,sec,0, dUT1, dAT );
elseif dt < 0
    dtd = 0;
    dts = -dt;
    [year,mon,day,hr,mi,sec] = invjday ( satrecde.jdsatepoch );
    ep = satrecse.jdsatepoch;
    index = find(ismember(finals(:,1:3),[year-2000 mon day],'rows'));
    dUT1 = finals(index,11);
    [~, ~, ~, ~, ~, ~, ttt] =...
        convtime (year,mon,day,hr,mi,sec,0, dUT1, dAT );  
else
    dts = 0;
    dtd = 0;
    [year,mon,day,hr,mi,sec] = invjday ( satrecse.jdsatepoch );
    index = find(ismember(finals(:,1:3),[year-2000 mon day],'rows'));
    dUT1 = finals(index,11);
    [~, ~, ~, ~, ~, ~, ttt] =...
        convtime (year,mon,day,hr,mi,sec,0, dUT1, dAT );

end

[~,rs_0_teme,vs_0_teme] = sgp4(satrecse,dts,consts);
[~,rd_0_teme,vd_0_teme] = sgp4(satrecde,dtd,consts);

order = 106;
eqeterms = 2;
opt = 'a';

[rs_0,vs_0] = teme2eci(rs_0_teme',vs_0_teme',[],ttt,order,eqeterms,opt );
[rd_0,vd_0] = teme2eci(rd_0_teme',vd_0_teme',[],ttt,order,eqeterms,opt );
rs_0 = rs_0';
vs_0 = vs_0';
rd_0 = rd_0';
vd_0 = vd_0';

max_test_time = 14*24*60;
fprintf("done!\n");

%% solve
fprintf("Finding TCA...");
runtime_ = tic;
[TCA,Dmin,samp_,built_] = ChebANCAS (satrecse,satrecde,dt,max_test_time,consts);
runtime = toc(runtime_);
fprintf("done in %g seconds\n",runtime);

%% plot results
fprintf("plotting...");
times_for_plot = 0:0.5:max_test_time;
d_for_plot = zeros(1,length(times_for_plot));
v2_for_plot = zeros(1,length(times_for_plot));
fp_dot=zeros(1,length(times_for_plot));
for i = 1:length(d_for_plot) 
    [~, r1_teme,v1_teme] = sgp4(satrecse,times_for_plot(i)+dts,consts);
    [~, r2_teme,v2_teme] = sgp4(satrecde,times_for_plot(i)+dtd,consts);
    
    rd = r1_teme-r2_teme;
    rd_dot = v1_teme-v2_teme;
    d_for_plot(i) = norm(r1_teme-r2_teme); 
    v_for_plot(i) = norm(v1_teme-v2_teme);
    fp_dot(i) = 2*dot(rd_dot,rd);
    v2_for_plot(i) = v_for_plot(i)^2;
end
%%
figure
plot(times_for_plot,d_for_plot,'LineWidth', 1.5);
hold;
scatter(TCA,Dmin,'sc');
xlabel("time from epoch (mins)");
ylabel('distance(km)')
set(findall(gcf,'type','axes'),'fontsize',12)
set(findall(gcf,'type','text'),'fontSize',12)
xlim([0 times_for_plot(end)])
grid on

fprintf("Done! look at those beautiful results!\n");