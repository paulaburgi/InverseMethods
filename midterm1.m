% close; 
% t = [-5:.5:99.5]'; 
% %t = [0:.5:100]'; 
% n = length(t); 
% 
% t0 = 10; 
% g0 = t0.*exp(-1)*.01; 
% g = zeros(n,1);
% for i = 1:n
%     if i > find(t == t0)
%         g  = g0.*t.*exp(-t/t0); 
%     end
% end
% 
% G = zeros(n,n);
% 
% for i = 2:size(G,1)
%     for j = 1:size(G,2)
%         if j > i 
%             G(i,j) = g0.*(t(j)-t(i)).*exp(-(t(j)-t(i))./t0).*0.5;
%         end
%     end
% end
% %G = G(2:end,2:end);
% % G = fliplr(flipud(G)); 
% idx = find(t == 12); 
% s = zeros(n,1); s(1) = 1;
% s = exp((-(t-8).^2)/8)+0.5*exp((-(t-25).^2/(2*8))); 
% test = G*s;
% plot(t,s); hold on; 
% plot(t,test);
% %keyboard

