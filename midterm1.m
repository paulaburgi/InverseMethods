close; 
t = [-5:.5:99.5]'; 
%t = [0:.5:100]'; 
n = length(t); 

t0 = 12; 
g0 = exp(1)/t0;
g = zeros(n,1);
for i = 1:n
    if i > find(t == t0)
        g(i)  = g0.*t(i).*exp(-t(i)/t0); 
    end
end
%plot(t,g)
G = zeros(n,n);

for i = 1:size(G,1)
    for j = 1:size(G,2)
        if i >= j 
            G(i,j) = g0*(t(i)-t(j)).*exp(-(t(i)-t(j))./t0).*0.5;
        end
    end
end

idx = find(t == 12); 
s = zeros(n,1); s(1) = 1;
%s = exp((-(t-8).^2)/8)+0.5*exp((-(t-25).^2/(2*8))); 
    s1 = normpdf(t,7,2.5);
    s1 = s1/max(s1); 
    s2 = normpdf(t,20,2.5); 
    s2 = (s2/max(s2))*.38; 
    s  = [s1(1:39); s2(40:end)];
% G = G(2:end,1:end-1);
% s = s(2:end);
% t = t(2:end); 
test = G*s;
plot(t,s); hold on; %axis([t(1) t(end) 0 6]);
plot(t,test);
%keyboard

[u, s, v] = svd(G);
sd = diag(s); 