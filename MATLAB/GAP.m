%Greedy Analysis Pursuit (GAP) algorithm for non-linear sparse recovery of
%This function returns a column vector recovered using non-linear GAP algorithm

function x = GAP(y,f,D,N)
[p, d] = size(D);
r=@(x)(y-f(x));
fun=@(x)[r(x);-D*x];
options=optimset('Display','iter-detailed'); 
o = rand(p,1);
x = lsqnonlin(fun,o,[],[],options);%initial estimate of x
n1 = norm(y-f(x));
supp = (1:p); %initialize co-support of x

j=1;

while(j<=N && n1>10^(-3))
    c = abs(D*x);
    [r, i]= max(c); %find the largest entry of Dx
    supp(i)=[]; %update co-support by eliminating the index corresponding to largest entry
    D(supp,:); %update analysis operator
    fun=@(x)[r(x);-D*x];%function handle for ||y-f(x)||
                        %                    || -Dx  ||
    x=lsqnonlin(fun,x,[],[],options);%solve non-linear least squares problem to update x
    j=j+1;
    n1=norm(y-f(x));
end