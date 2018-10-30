% load image
close all
pics = imread('cheese.jpg');
pics2 = flipud(double(pics(:,:,3)));
pics2 = pics2(1:5:end,1:5:end); 
subplot(2,2,1); pcolor(pics2); shading flat; colormap gray; caxis([0 255]); 

[n,m] = size(pics2); 
npts = n*m; 
idim = n; 


sig = 1.5; 
wid = 3; 
scale = 0.01.*(max(max(pics2))); 
a = [exp(-([0:wid-1].^2)./(2*sig.^2)), zeros(1,idim-wid)]; 
g = toeplitz(a); 
g = sparse(g); 
g = (1/(2.*pi.*sig.^2)).*kron(g,g); 
gs = sparse(g); 

d = gs*reshape(pics2,npts,1)+(scale.*randn(npts,1)); 
di = reshape(d, n,m); 

subplot(2,2,2); pcolor(di); shading flat; colormap gray; caxis([0 255]); 

m0 = zeros(npts,1); 
%m0 = sparse(m0); 
s0 = d-(gs*m0); 
r0 = gs'*s0; 
p0 = r0; 
q0 = gs*p0;

nints = 50; 
mk_all = [];
m_all = [];
d_all = [];
for k=1:nints
    if k == 1
        rk = r0;
        qk = q0;
        mk = m0;
        pk = p0;
        sk = s0;
    end
    
    
    ak1 = (rk'*rk)./(qk'*qk);
    mk1 = mk+(ak1*pk);
    sk1 = sk-(ak1*qk);
    rk1 = (gs'*sk1); 
    bk1 = (rk1'*rk1)./(rk'*rk);
    pk1 = rk1+(bk1'*pk); 
    qk1 = (gs*pk1); 
    
    rk = rk1;
    qk = qk1;
    mk = mk1;
    pk = pk1;
    qk = qk1;
    sk = sk1;
    
    % model length
    m_all = [m_all; sqrt(sum(mk).^2)];
    % data residual
    gm = gs*mk;
    d_all = [d_all; sqrt(sum((rk).^2))];
    % all models
    mk_all = [mk_all mk]; 
    
    
%     figure(2); 
%     pcolor(reshape(rk1,n,m)); shading flat; colormap gray; %caxis([0 255]); 
%     keyboard
    
end

subplot(2,2,3); 
plot(d_all, m_all, 'o-');
% plot(log10(d_all), log10(mk_all), 'o-');

idx = 10; 
mkp = mk_all(:,idx); 


subplot(2,2,4); pcolor(reshape(mkp,n,m)); shading flat; colormap gray; caxis([0 255]); 


        
 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
        