function p_m = probabilityC(n,r)
p_C = 0;
for k=1:n
    p_i = k^-r ;
    p_C = p_i+p_C;
end
 C=1/p_C;
% fprintf('C=''%s''.\n',C)
 
p_m=C*[1:n].^(-r);

%p(k)=k^{-gamma}=1/C要輸入多少n得多少C,得C=0.8020,gamma=2.8