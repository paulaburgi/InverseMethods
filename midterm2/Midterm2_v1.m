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
h = hm-hd; 

fiS = int(h, S); 
fiT = int(h, T); 
J = [fiT fiS];

Cd = eye(n)*s; 


T0    = .7; %0.584; 
S0    = 1e-3; %2.1e-3;
l0    = 1000; 
nints = 100; 
v     = 0.1; 
T = T0; 
S = S0; 
h0 = eval(subs(hd));

figure(1); hold on; 
plot(t, hm, '*'); 
T = .585; 
S = .00207; 
hb = eval(subs(hd));
plot(t, hb, 'k.', 'markersize', 35); 
plot(t, h0, 'r.', 'markersize', 40); 
cmap = jet(nints); 

for k = 1:nints
    if k == 1
        T = T0; 
        S = S0; 
        l = l0; 
        h1 = (eval(subs(h))-hm); 
        ri = h1'*h1; 
        uJ = 1; 
    end
    
    %if uJ == 1
        Ji = eval(subs(J)); 
    %end
    Gmi = -eval(subs(h))-hm; 
 

    fi = (Gmi-hm); 
    
    
    xi = fi'*fi; 
%     figure(2); 
%     plot(S,T, 'rx', 'markersize', 2); 
        
    figure(1); 
    plot(t, Gmi, '.', 'color', cmap(k,:), 'markersize', 15); 
    
    Hi  = (Ji'*Ji)+(eye(m).*l);
    dm = -inv(Hi)*Ji'*fi; 

    
    rn = (fi'*fi); 
    rs = ri-rn; 
    if rs < 0 % if current guess is worse
        l = l/v; 
        uJ = 0; 
    elseif rs >= 0 % if current guess is better
        l = l*v;
        uJ = 1; 
        
    end
    T = T+dm(1); 
    S = S+dm(2); 
    ri = (fi'*fi); 
    
    %%%%%%% TRY PLOTTING EACH GUESS ON CHI2 PLOT %%%%%%%%%
    
end

plot(t, Gmi, 'g.', 'markersize', 40); 
axis([4 51 0 2]); 

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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        



