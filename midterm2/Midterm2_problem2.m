% Paula Burgi
% Midterm 2, Problem 2

% model parameters
n  = 40; 
m  = 40; 
L  = 40; 
%x = linspace(1,L,n); 
s = linspace(1,L,m)'; 
ss = sym('ss'); 
x = sym('x'); 
%h = sym('h'); 
xm = linspace(1,L,n); 

% calculate true height model
h1 = (2.*pi.*(x-(L./4)))./L;
h2 = (2.*pi.*(x-((3.*L)./4)))./L;
h  = 15 - (5./2).*tanh(h1) + (5/2).*tanh(h2); 

d1 = h./(((ss-x).^2 + h.^2).^(3/2)); 
% d1 = ((5*tanh((pi*(x - 30))/20))/2 - (5*tanh((pi*(x - 10))/20))/2 + 15)/(((s-x)^2 + ((5*tanh((pi*(x - 30))/20))/2 - (5*tanh((pi*(x - 10))/20))/2 + 15)^2)^(3/2))

% % calc dn
%     si = linspace(1,L,m); 
%     xi = linspace(1,L,n); 
%     dx = xi(2)-xi(1); 
%     da = [];
%     ii = []; 
%     for i = 1:m
%         ss = si(i);
%         ii = []; 
%         for j = 1:n
%             x = xi(j); 
%             ii = [ii; eval(subs(d1))]; 
%         end
%         da = [da; sum(ii).*dx]; 
%         disp(i);
%     end
%     st  = max(da).*0.01; 
%     dn  = da + randn(40,1).*st; 

% calculate function & jacobian
d2 = h./(((s-x).^2 + h.^2).^(3/2)); 
clear dii
clear di
clear J
for i = 1:m
    xii = ['x' num2str(i)]; 
    syms(xii);
    di = subs(d2,xii); 
    if i == 1
        dii = di;
        J   = repmat(di, 1, m); 
    else 
        dii  = dii+di; 
        J(:,i) = subs(J(:,i), xii); 
    end
end
hdd = dii - dn; 

h1 = (2.*pi.*(xm-(L./4)))./L;
h2 = (2.*pi.*(xm-((3.*L)./4)))./L;
hm  = 15 - (5./2).*tanh(h1) + (5/2).*tanh(h2); 

Cd = eye(n);

%var0 = hm - 0.5; 
var0  = (1.1:40.1)'; 

[x, k, Cm, X2] = LMLSQ(hdd, Cd, var0, J);

 % calc dest
xf = x; 
h1 = (2.*pi.*(xf-(L./4)))./L;
h2 = (2.*pi.*(xf-((3.*L)./4)))./L;
hf  = 15 - (5./2).*tanh(h1) + (5/2).*tanh(h2);  

figure; hold on; 
plot(s,x,'b-'); 
plot(s,var0,'r-'); 
plot(s,hm,'k-'); 

keyboard
 
    


% da = [];
% ii = []; 
% for i = 1:m
%     s = si(i);
%     ii = []; 
%     for j = 1:n
%         x = xi(j); 
%         ii = [ii; sum(eval(subs(d1)))]; 
%     end
%     da = [da; sum(ii).*dx]; 
%     disp(i);
% end
% 
% st  = max(da).*0.01; 
% dn  = da + randn(40,1).*st; 
% figure; 
% plot(xm, dn, 'k-'); 

keyboard








da = [];
ii = []; 
for i = 1:m
    s = si(i);
    ii = []; 
    for j = 1:n
        x = xi(j); 
        ii = [ii; eval(subs(d1))]; 
    end
    da = [da; sum(ii).*dx]; 
    disp(i);
end

st  = max(da).*0.01; 
dn  = da + randn(40,1).*st; 

% plot true height model
% plot(x, h, '.k-'); 
figure; 
plot(xm, dn, 'k-'); 
