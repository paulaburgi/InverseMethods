% close; 
% t = [-5:.5:99.5]'; 
% %t = [0:.5:100]'; 
% n = length(t); 
% 
% t0 = 12; 
% g0 = exp(1)/t0;
% g = zeros(n,1);
% for i = 1:n
%     if i > find(t == t0)
%         g(i)  = g0.*t(i).*exp(-t(i)/t0); 
%     end
% end
% %plot(t,g)
% G = zeros(n,n);
% 
% for i = 1:size(G,1)
%     for j = 1:size(G,2)
%         if i >= j 
%             G(i,j) = g0*(t(i)-t(j)).*exp(-(t(i)-t(j))./t0).*0.5;
%         end
%     end
% end
% 
% idx = find(t == 12); 
% s = zeros(n,1); s(1) = 1;
% %s = exp((-(t-8).^2)/8)+0.5*exp((-(t-25).^2/(2*8))); 
%     s1 = normpdf(t,7,2.5);
%     s1 = s1/max(s1); 
%     s2 = normpdf(t,20,2.5); 
%     s2 = (s2/max(s2))*.38; 
%     s  = [s1(1:39); s2(40:end)];
% % G = G(2:end,1:end-1);
% % s = s(2:end);
% % t = t(2:end); 
% test = G*s;
% plot(t,s); hold on; %axis([t(1) t(end) 0 6]);
% plot(t,test);
% %keyboard
% 
% [u, s, v] = svd(G);
% sd = diag(s); 





%% problem 2
close; 

n = 160; 
m=160; 
d = linspace(0,40,n'); 
st = 2; %slowness at top of well (s/km)
sb = 5; %slowness at bottom of well (s/km)
s = linspace(st,sb,n)'; % linear slowness increase linearly
s((3/8)*n:(5/8)*n) = 2; % modify s.t. n=15-25 = 2 s/m

G = tril(ones(m,n).*0.5); 
d1 = G*s; 

ns   = rand(m,1).*.1; % random noise with std=0.1
dns1 = d1+ns; % add noise to data
 
figure; hold on; 
plot(d, s); 
mest = inv(G'*G)*G'*dns1; 
plot(d, mest); 


[Usvd,Ssvd,Vsvd] = svd(G); 
r = rank(G);
% r = 100; 
Vsvd = Vsvd(:,1:r); 
Ssvd = Ssvd(1:r,1:r); 
Usvd = Usvd(:,1:r); 
Gsvd    = Vsvd*inv(Ssvd)*Usvd'; 
 mest3    = Gsvd*dns1;

plot(d,mest3, 'b'); 


sd = diag(Ssvd); 
a2 = 1; 
Fd = (sd.^2)./((sd.^2)+(a2.^2)); 
Fd = eye(r,r).*Fd; 

mest4 = Vsvd*Fd*inv(Ssvd)*Usvd'*dns1; 

plot(d,mest4,'r'); 


L = diff(diff(eye(n,n))); 
L = [[-1 1 zeros(1,n-2)]; L]; 
L = [L; [zeros(1,n-2) -1 1]]; 

[U,V,X,C,S] = gsvd(G,L); 
C = C(1:r,1:r); 
S = S(1:r,1:r); 
U = U(:,1:r); 
X = X(:,1:r); 


g = sqrt(diag(C'*C)./diag(S'*S));
a = 1; %1e20; %1e17; 
F = (g.^2)./((g.^2)+(a.^2)); 
F = eye(r,r).*F; 

Gi = inv(X)'*F*inv(C)*U';  

mest2 = (Gi)*dns1; 

plot(d,mest2, 'k'); 





%% problem 3

close all; 

% set r to be variable
r = sym('r');
% mass and inertia equations
Md = (4.*pi.*r.^2); 
Id = ((8/3).*pi.*r.^4); 
G  = [Md; Id]; 
M  = int(Md,r); 
I  = int(Id,r); 

%EARTH
% parameters for earth (data)
Re = 6.3708e6; 
Me = 5.972e24;
Ie = 8.034e37; 
d  = [Me; Ie];

% calculate q
q1 = eval(int(Md,r, 0, Re)); 
q2 = eval(int(Id,r, 0, Re)); 
q  = [q1; q2];

% set up to interate over all radii 
Ri = [0:50000:Re]; 
mest_all = []; 
k_all = [];

% interate over radii
for i=1:length(Ri)
    % radius
    ri = Ri(i); 
    % calculate H
    H0 = G*G'*((r-ri).^2); 
    H = eval(int(H0, 0, Re)); 
    %calculate c
    c = (inv(H)*q)./(q'*inv(H)*q);  
    % calculate mest
    mest = c'*d; 
    mest_all = [mest_all; mest]; 
    % calculate k
    k(r) = c'*G; 
    k_all = [k_all; eval(k(Ri))]; 
end

figure('units','normalized','outerposition',[0 0 .7 .5]); 
subplot(1,2,1)
plot(mest_all, Ri./1e3); 
title('Earth Density'); 
xlabel('Density (g/cm^3)'); 
ylabel('Radial Distance from center (km)'); 
xlim([3800 7800]); 
subplot(1,2,2)
pcolor(flipud(k_all)); shading flat;
colormap gray
title('Kernel K(r,rh)'); 

% MARS 
% parameters for mars (data)
Re = 3.389e6; 
Me = 0.642e24;
Ie = 2.709e36; 
d  = [Me; Ie];

% calculate q
q1 = eval(int(Md,r, 0, Re)); 
q2 = eval(int(Id,r, 0, Re)); 
q  = [q1; q2];

% set up to interate over all radii 
Ri = [0:50000:Re]; 
mest_all = []; 
k_all = [];

% interate over radii
for i=1:length(Ri)
    % radius
    ri = Ri(i); 
    % calculate H
    H0 = G*G'*((r-ri).^2); 
    H = eval(int(H0, 0, Re)); 
    %calculate c
    c = (inv(H)*q)./(q'*inv(H)*q);  
    % calculate mest
    mest = c'*d; 
    mest_all = [mest_all; mest]; 
    % calculate k
    k(r) = c'*G; 
    k_all = [k_all; eval(k(Ri))]; 
end

figure('units','normalized','outerposition',[0 0 .7 .5]); 
subplot(1,2,1)
plot(mest_all, Ri./1e3); 
title('Mars Density'); 
xlabel('Density (g/cm^3)'); 
ylabel('Radial Distance from center (km)'); 
xlim([3300 4700]); 
subplot(1,2,2)
pcolor(flipud(k_all)); shading flat;
colormap gray
title('Kernel K(r,rh)'); 








































