% syms x
% syms f(x,y)
% syms s
% syms t
% 
% Rdd = dirac((x.*cos(t) + y.*sin(t)) - s);
% Rd  = int(Rdd,'x')
% R   = int(Rd,'y')

% img = imread('phantom_test.png'); 
% img = double(img(:,:,1)); 

% rim = img; 
% % da = floor(size(img, 1)/64); 
% % dr = floor(size(img, 2)/120); 
% % rim = rim(1:da:end, 1:dr:end); 
% R = []; 
% for i = 1:180
%     rim = imrotate(img,i,'bicubic','crop'); 
%     %rim = [rim(end, :); rim(1:end-1,:)];
%     
%     R(:,i) = sum(rim); 
% end
% %close all; 
% pcolor(R); shading flat; colormap gray;

t = 1; 
nx = size(img,1); 
ny = size(img,2);  
ox = floor(nx/2); 
oy = floor(ny/2); 

[x,y] = meshgrid(1:nx, 1:ny); 
x = x - ox; 
y = y - oy; 
r = sqrt(x.^2 + y.^2); 
t = atan2(y,x); 
ri = linspace(min(min(r)), max(max(r)), nx); 
ti = linspace(min(min(t)), max(max(t)), ny); 
[tg, rg] = meshgrid(ti, ri); 
xi = rg.*cos(tg) + ox; 
yi = rg.*sin(tg) + oy; 


cds = [xi(:) yi(:)]; 

imP = ImToPolar (img, 0, 1, 200, 180);


R = []; 
rim = imP; 
for i = 1:size(imP,2)
    rim = [rim(:,end) rim(:,1:end-1)];
    rimi = reshape(rim, l, 1); 
    
    
    
    R(i,:) = sum(rim); 
end
pcolor(R); shading flat; colormap gray;











