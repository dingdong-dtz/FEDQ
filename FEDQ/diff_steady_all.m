function c = diff_steady_all(xbar,x,k,m)
%程序要求x以行向量形式输入
%求用x处函数值表达xbar处k阶导数系数,且具有m阶收敛精度，要求x的长度必须大于等于k+m，
%当x的长度为k+m时，退回到普通方法
%输出格式与x相同
L=length(x);
%m=min(m,L-k);%当输入的阶数要求过大时，即退回到普通方法
if L<=k+m
    c = fdcoeffF(k,xbar,x);
else
    N=k+m;%约束条件的数目
    h=x-xbar;
    factor=mean(abs(diff(h)));
    h=h/factor;%通过近似归一化的操作，使得矩阵的奇异性降低，使得计算更加准确
    H=zeros(1,2*N-1);
    for i=1:2*N-1
        H(i)=sum(h.^(i-1));%程序中将0^0计算为1
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

