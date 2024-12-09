function c = diff_steady_all(xbar,x,k,m)
%The program requires x to be input in the form of a row vector
%To express the k-order derivative coefficient at xbar using the function value at x, with m-order convergence accuracy, it is required that the length of x must be greater than or equal to k+m,
%When the length of x is k+m, return to the normal method
%Output format is the same as x
L=length(x);
if L<=k+m
    c = fdcoeffF(k,xbar,x);
else
    N=k+m;%The number of constraints
    h=x-xbar;
    factor=mean(abs(diff(h)));
    h=h/factor;%By performing approximate normalization, the singularity of the matrix is reduced, making the calculation more accurate
    H=zeros(1,2*N-1);
    for i=1:2*N-1
        H(i)=sum(h.^(i-1));
    end
    A=zeros(N,N);
    j=1:N;
    for i=1:N
        A(i,:)=H(1,i:i+N-1)./factorial(j-1)/factorial(i-1);
    end
    B=zeros(N,1);B(k+1)=1;
    D=A\B;
    M=zeros(N,L);
    for i=1:N
        M(i,:)=h.^(i-1)/factorial(i-1);
    end
    c=D'*M;
    c=c/factor^k;
end
end

