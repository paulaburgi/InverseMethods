figure; hold on; box on;
plot(xn, h, 'k', 'linewidth', 3); 

for i=1:size(varfra, 2)
    plot(xn, varfra(:,i), 'b'); 
    keyboard
    oops
end
