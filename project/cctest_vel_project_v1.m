%% Load data, build basic G

close all
clear 

ccdir = '/data/pmb229/other/clearcuttingTStest/'; 

% load baseline data & int data
    load([ccdir 'baselines.mat']); 
    load([ccdir 'meancor_bl_dates_area2_HH.mat']); 
    l = 0.236; sr = 8.5e5; los = 38.7; 

% extract useful info
    dn_all = baselines2.dn_all; 
    bl_all = baselines2.bl_all; 
    d      = meancor_bl_dates; 
    gidx   = d.good_cor_idx; 
    dc     = d.dateCombos; 
    bl     = abs(d.bl); % have to make bl positive for inversion to work
    dta    = diff(dn_all); 


% extra ints to connect with network
    exints = [733236 733604; 733742 733972; 733788 733972];
    [a,b]  = dsearchn(dc, exints);
     %gidx   = [gidx; a]; 

% select which ints to work with
    bl     = bl(gidx,:); 
    dc     = dc(gidx,:); 
    dcv    = dc(:);
    du     = sort(unique(dcv));
    nvels  = length(du)-1;
    nints  = length(dc); 

% build G matrix
    G  = zeros(length(dc), length(du)-1);
    xa = [];
    ya = [];
    for i = 1:length(G)
        xi     = dc(i,1); 
        yi     = dc(i,2); 
        x = find(du == xi); 
        y = find(du == yi); 
        
        dtj   = dta(x:y-1);
        for j = 1:length(dtj)
            G(i,x+j-1) = dtj(j); 
        end

        xa = [xa; x];
        ya = [ya; y];
    end
    d2y = 0.00274; % convert days to years
    Ga = G.*d2y; 
    
    % deal with unconstrained dates 
    rdi = find(sum(G)' == 0);
    nrdi = length(rdi);
    G  = [Ga; ones(nrdi, nvels)*(1/(nvels-1))];
    rng = 1:nrdi; 
    for i = rng
        idx = i-1; 
        G(end-idx, rdi(i)) = -1; 
    end


%% Compare data to inv with and without DEM Error
% noise weight
nw  = 0.1; 
n   = [(rand(length(dc), 1)-0.5).*nw; zeros(nrdi, 1)];

% define time steps and "real" deformation
vel  = ones(nvels, 1).*3; 
def  = Ga*vel; 

% create ints
intsr = [def; zeros(nrdi, 1)];
ints  = intsr+n;

% add baseline term to G matrix
blg    = [(4.*pi.*bl)./(l.*sr.*sind(los)); zeros(nrdi, 1)]; 
%blg    = [bl; zeros(nrdi, 1)]; 
Gbl    = [G abs(blg)]; 
dz     = 30; 
velz   = [vel; dz];
intsrz = intsr+abs(blg.*dz); 
intsz  = intsrz+n; 

% calc Gsvd with baseline term 
r2      = rank(Gbl); 
[U,S,V] = svd(Gbl); 
Gsvd    = V(:,1:r2)*inv(S(1:r2,1:r2))*U(:,1:r2)'; 
mzsvd   = Gsvd*intsz; 
mvel0   = mzsvd(1:end-1); 
mz0     = mzsvd(end); 
    

%% Add it DEM error with each time step 

dz     = 30;           % topo error
iz     = abs(blg.*dz); % addition to phase due to topo error
intsiz = intsr;        % initiate ints where only some have additional topo error phase
c      = 10;           % date index after which topo error is introduced

for j = 1:length(dc)
    % add phase to ints that include dates greater than j
    if xa(j) > c
        intsiz(j) = intsiz(j)+iz(j); 
    end
end

intsiz = intsiz + n; 

mest     = Gsvd*intsiz;
mbl      = mest(end); 
mvel     = mest(1:end-1); 
mz       = mest(end); 


%% nonlinear inversion 

% make symbolic variables
z = sym('z'); 
v = [];
for i = 1:nvels
    if i < 10
        vi = ['v0' num2str(i)]; 
    else
        vi = ['v' num2str(i)]; 
    end
    vi = sym(vi); 
    v = [v; vi]; 
end
svar  = [v; z];

% create residual function
G_m_ = G*v;
for j = 1:length(dc)
    if xa(j) > c
        G_m_(j) = G_m_(j)+blg(j)*z; 
    end
end
% covariance matrix
Cd   = eye(nints+nrdi).*(nw^2); 
Cdi  = inv(sqrt(Cd)); 
% residual function
f = Cdi*(G_m_ - intsiz); 

% create Jacobian 
J = [];
for i = 1:nvels+1
    if i == nvels+1
        Ji  = diff(f, z); 
    else
        vii = v(i); 
        Ji  = diff(f, vii);
    end
    J   = [J Ji]; 
end

% eps
%var0  = zeros(nvels+1, 1); 
var0  = 1:nvels+1; 
ep    = 1e-3; 

% LM method
varf   = LMLSQ_project(f, var0, J, ep, svar); 
mnlvel = varf(1:end-1); 
mnlz   = varf(end); 



%% plot

figure; hold on; box on; 
dd     = diff(du); 
vdates = du(1:end-1)+dd;  
plot(vdates, zeros(nvels, 1), 'k', 'HandleVisibility','off'); 
plot(vdates, [mvel0], 'linewidth', 2);
plot(vdates, [mvel], 'linewidth', 2);
plot(vdates, [mnlvel], 'linewidth', 2); 
yl = ylim; 
plot([vdates(c); vdates(c)], [-100 100], 'k--'); 
plot(vdates(rdi), varf(rdi), '*'); 
xlabel('Date'); 
ylabel('Deformation'); 
datetick('x'); 
xlim([vdates(1) vdates(end)]);
ylim(yl); 
p = 3; % precision
legend(['no error (lin), z=' num2str(mz0, p) 'm'], ...
       ['dem error (lin), z=' num2str(mz, p) 'm'], ...
       ['dem error (non lin), z=' num2str(mnlz, p) 'm'], ... 
       ['dem error start, z=' num2str(dz, p) 'm'], ... 
       'missing dates', ...
       'location', 'southeast'); 




























    
    
    
    
    
    
    
    
