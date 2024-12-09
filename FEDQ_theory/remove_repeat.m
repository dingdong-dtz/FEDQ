function  [jx,jy]= remove_repeat(jx,jy)
%[jx,jy]= remove_repeat(jx,jy)
%   jx,jy表示一条曲线，可能是闭合的也可能不是闭合的，本函数去除相邻的重复值和首位重复值
index_repeat=~or(abs(diff([jx jx(1)]))>1e-10,abs(diff([jy jy(1)]))>1e-10);
jx(index_repeat)=[];jy(index_repeat)=[];
end