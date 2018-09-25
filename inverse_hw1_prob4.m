% Paula Burgi
% Homework 1, problem 4

db = 40e-3; %depth of well, (km) 
st = 2; %slowness at top of well (s/km)
sb = 5; %slowness at bottom of well (s/km)
dd = 0.5e-3; %interval of slowness measurements (km)

n = 80;  % number of data points
dpt = linspace(dd,db,n)'; % depth intervals
m = 80; % number of model params

s = linspace(st,sb,n)'; % linear slowness increase linearly
% PROB 4A
    s((3/8)*n:(5/8)*n) = 2; % modify s.t. n=15-25 = 2 s/m
% PROB 4B
    G = tril(ones(m,n).*0.5); 
% PROB 4C
    d1 = G*s; % data generated for n=80
% PROB 4D
    ns   = rand(m,1).*0.1; % random noise with std=0.1
    dns1 = d1+ns; % add noise to data
% PROB 4E 
    slsq1 = inv(G'*G)*G'*dns1; % lsq solution
    sdd1  = diff(dns1)/0.5; % take difference and normalize to dd


% Sample 160 points
n = 160;  % number of data points
dpt2 = linspace(dd,db,n)'; % depth intervals

% PROB 4B
    G = tril(ones(m,m).*0.5); 
    G = reshape([G(:) G(:)]',2*size(G,1), []); % interleaf elements of G to reflect over sampling
    for i = 1:2:length(G) 
        G(i,(i+1)/2) = .25; % set elements at 1/2 steps equal to half the depth interval 
    end
% PROB 4C
    d = G*s; % data generated for n=160
% PROB 4D
    ns = rand(n,1).*0.1; % random noise with std=0.1
    dns = d+ns; % add noise to data
    
% PROB 4E 
    slsq = inv(G'*G)*G'*dns; % lsq solution
    sdd = diff(dns)/.25; % take difference and normalize to dd

% Plot 
figure('units','normalized','outerposition',[.2 .2 .7 .5]);  hold on; 
subplot(1,2,1); hold on; box on; 
plot(dpt*1e3, s, '-k', 'linewidth', 2); 
plot(dpt*1e3, slsq1, '-r.'); 
plot(dpt*1e3, [s(1); sdd1], '--b.'); 
xlabel('depth (m)'); 
ylabel('slowness (s/km)'); 
title('N=80 data points'); 

subplot(1,2,2); hold on; box on; 
plot(dpt*1e3, s, '-k', 'linewidth', 2); 
plot(dpt2*1e3, [s(1); sdd], '-b.'); 
plot(dpt*1e3, slsq, '-r.'); 
xlabel('depth (m)'); 
ylabel('slowness (s/km)'); 
title('N=160 data points'); 
legend('True slowness', 'Est Slowness (LSQ)', 'Est. Slowness (diff)', ...
    'location', 'northwest'); 





