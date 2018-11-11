% Paula Burgi
% Midterm 2, Problem 1

clear
close all

t = [5 10 20 30 40 50]'; % from aster
hm = [0.72 0.49 0.30 0.20 0.16 0.12]'; % from aster
%h = [0.80 0.78 0.54 0.41 0.33 0.27]'; % for midterm

n = length(t); 
m = 2; 

T = sym('T'); % m^2/hr 
S = sym('S'); % dimensionless
Q = 50;       % m^3
d = 60;       % m
s = 0.01; 

h1 = Q./(4.*pi.*T.*t); 
h2 = exp((-(d.^2).*S)./(4.*T.*t)); 
hd  = h1.*h2;
h = hd-hm; 

fiS = diff(h, S); 
fiT = diff(h, T); 
J = [fiS fiT];

Cd = eye(n)*(s.^2); 

T0    = .7; %0.584; 
S0    = 1e-3; %2.1e-3;
l0    = 1000; 
nints = 90; 
v     = 0.1; 
T = T0; 
S = S0; 
var0 = [S0; T0];
h0 = eval(subs(hd));

figure(1); hold on; 
plot(t, hm, '*'); 
T = .585; 
S = .00207; 
hb = eval(subs(hd));
plot(t, hb, 'k.', 'markersize', 35); 
plot(t, h0, 'r.', 'markersize', 40); 
cmap = jet(nints); 

ep = 1e-10; 
c_all = [];
e_all = []; 

lmvars = symvar(h); 
nv     = length(lmvars); 

for k = 1:nints
    if k == 1
        for j = 1:nv
            varstr = char(lmvars(j)); 
            evalc([varstr '=' num2str(var0(j))]); 
        end
        T = T0; 
        S = S0; 
        l = l0; 
        %h1 = (eval(subs(h))-hm); 
        %ri = h1'*h1; 
        fi = sqrt(Cd)*eval(subs(h)); 
        ri = (fi'*fi); 
        ci = 0;
    end
    
    Ji = eval(subs(J)); 
     Gmi = eval(subs(hd)); 
    fi = eval(subs(h)); 

    %fi = (Gmi-hm); 
    
    % convergence
    c  = 2.*Ji'*fi; 
    dc = mean(abs(c-ci)); 
    de = mean(ep.*(1+abs(c))); 
    c_all = [c_all; (dc)];
    e_all = [e_all; (de)]; 
    if dc < de
        break % stops the for-loop based on convergence criteria
    end
    
        
%     figure(1); 
%     plot(t, Gmi, '.', 'color', cmap(k,:), 'markersize', 15); 
    
    Hi  = (Ji'*Ji)+(eye(m).*l);
    dm = -inv(Hi)*Ji'*fi; 

    
    rn = (fi'*fi); 
    rs = ri-rn; 
    if rs < 0 % if current guess is worse
        l = l/v; 
    elseif rs >= 0 % if current guess is better
        l = l*v;
    end
    
    for j = 1:nv
        varstr = char(lmvars(j)); 
        evalc([varstr '=' varstr '+' num2str(dm(j))]); 
    end
    
%     T  = T+dm(1); 
%     S  = S+dm(2); 
    ri = (fi'*fi); 
    ci = c; 
    
    
end

plot(t, Gmi, 'g.', 'markersize', 40); 
axis([4 51 0 2]); 

x = eval(lmvars); 
    


test = [c_all e_all];

keyboard














%% test
Ti = [.1:.1:2]'; 
Si = [3e-4:3e-4:6e-3]';
X2 = []; 
h_all = [];
idx = [];
for i=1:length(Ti)
    for j = 1:length(Si)
        T = Ti(i); 
        S = Si(j); 
        hi = eval(subs(h)); 
        h_all = [h_all hi]; 
        r  = (hm-hi)./s; 
        xi = r'*r; 
        X2 = [X2; xi]; 
        idx = [idx; i j];
    end
end

close all; 
figure; hold on; 
plot(t, hm, '*'); 
T = .585; 
S = .00207; 
hb = eval(subs(h));
plot(t, hb, 'ko'); 
plot(t, h_all(:,107), 'xg'); 
        
        
X2r = reshape(X2, length(Si), length(Ti))'; 

close all; 
figure(2); hold on; 
[col,row] = meshgrid(Si,Ti);
[M,c] = contour(col, row, X2r, [10 100 1000 3000], 'ShowText','on', 'color', 'k');
c.LineWidth = 2;
xlabel('Storage coef (S)'); 
ylabel('Transmissivity (T)'); 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        



