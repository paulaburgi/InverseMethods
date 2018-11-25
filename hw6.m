% Paula Burgi
% HW 6, #2 

function [lamda, fi] = loadedDice(xt)


% loop through possible values of lamda
li = -10:.001:10; 
dc = 1:6;

% start loops
xi = [];
for i = 1:length(li)
    l = li(i); 
    x1 = [];
    x2 = [];
    % create vectors for x to sum
    for j=dc
        x1 = [x1; j.*exp(-l.*j)]; 
        x2 = [x2; exp(-l.*j)]; 
    end
    x = sum(x1)./sum(x2); 
    xi = [xi; x]; 
end

% find closest value of lamda to match desired x
[k, d] = dsearchn(xi, xt); 
lamda = li(k); 

% find Z
Zi = []; 
for j=dc;
    Zi = [Zi; exp(lt.*j)]; 
end
Z = sum(Zi); 

% find fi
fi = exp(-lt*dc)./Z'; 

















