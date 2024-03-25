% <k> := k_n 
function k_em = k_1n(i,k,r) % i is power of k,  n is 100 or even more
k_em = sum([1:k].^i.*probabilityC(k,r))  ;       % k^i*p(k), if i = 1 then k*p(k).  if i = 2 then (k^2)*p(k)   
    
