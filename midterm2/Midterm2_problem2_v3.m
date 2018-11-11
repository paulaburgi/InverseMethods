% Paula Burgi
% Midterm 2, Problem 2

clear

n = 40; 
m = 40; 
L  = 40; 

sn = linspace(1,40,n)'; 
xn = linspace(1,40,n)';
dx = xn(2)-xn(1); 
hs = sym('hs'); 
xs = sym('xs'); 

h1 = (2.*pi.*(xn-(L./4)))./L;
h2 = (2.*pi.*(xn-((3.*L)./4)))./L;
h  = 15 - (5./2).*tanh(h1) + (5/2).*tanh(h2); 

df  = hs./(((sn-xs).^2 + hs.^2).^(3/2)); 

% manually sum integral in gravity function
d_all = []; 
for i = 1:m
    if i < 10
        hsi = ['hs0' num2str(i)]; 
    else
        hsi = ['hs' num2str(i)]; 
    end
    sym(hsi); 
    di = subs(df, hs, hsi);
    di = subs(di, xs, i); 
    if i == 1
        d_all = di.*dx; 
    else
        d_all = d_all+di; 
    end
end

% create fake data & Jacobian
J = [];
for i = 1:m
    if i < 10
        hsi = ['hs0' num2str(i)]; 
    else
        hsi = ['hs' num2str(i)]; 
    end
    evalc([hsi '=' num2str(h(i))]);
    Ji  = diff(d_all, hsi);
    J   = [J Ji]; 
end

dnn = eval(subs(d_all)); 
st  = max(dnn).*0.01; 
dn  = dnn + randn(40,1).*st; 
% figure; plot(xn, dn);

f    = d_all - dn; 


var0 = h - 0.5; 
%var0 = ones(n,1).*12.5; 
Cd   = eye(m); 

[varf, k, Cm, X2] = LMLSQ(f, Cd, var0, J); 
varf = varf';
vstr = char(symvar(f)); 


figure; hold on; 
plot(xn, h, 'k'); 
plot(xn, var0, 'r'); 
plot(xn, varf, 'b'); 

% forward model
for i = 1:m
    if i < 10
        hsi = ['hs0' num2str(i)]; 
    else
        hsi = ['hs' num2str(i)]; 
    end
    evalc([hsi '=' num2str(varf(i))]);
end
dest = eval(subs(d_all)); 
for i = 1:m
    if i < 10
        hsi = ['hs0' num2str(i)]; 
    else
        hsi = ['hs' num2str(i)]; 
    end
    evalc([hsi '=' num2str(var0(i))]);
end
d0  = eval(subs(d_all)); 



%% regularization 
a  = 1; 
L  = [-2 1 zeros(1, n-2); diff(diff(eye(n))); zeros(1, n-2) -2 1]; 
fa = a.*[f; L*d0]; 
K =  [J; a.*L]; 

[varfr, k, Cm, X2] = LMLSQ(fa, Cd, var0, K); 

% forward model
for i = 1:m
    if i < 10
        hsi = ['hs0' num2str(i)]; 
    else
        hsi = ['hs' num2str(i)]; 
    end
    evalc([hsi '=' num2str(varfr(i))]);
end
destr  = eval(subs(d_all)); 

figure; hold on; 
plot(xn, dn, 'k'); 
plot(xn, d0, 'r'); 
plot(xn, dest, 'b'); 
plot(xn, dest, 'g'); 


























