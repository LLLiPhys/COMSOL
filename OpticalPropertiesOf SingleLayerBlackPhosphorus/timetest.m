function timetest
clear all;
N=1e9;
u = rand(N,1);
v = zeros(N,1);
tic
v = u + 1;
toc
tic
for jj = 1:N
     v(jj) = u(jj)+1;
end
toc
tic
 for ii = 1:N
     v(ii) = u(ii)+1;
end
toc
end