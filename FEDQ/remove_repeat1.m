function  [jtheta,index_repeat]= remove_repeat1(jtheta)
%[jtheta,index_repeat]= remove_repeat1(jtheta)
%   jtheta表示一条曲线theta坐标，可能是闭合的也可能不是闭合的，本函数去除相邻的重复值和首位重复值
index_repeat=~abs(diff([jtheta jtheta(1)]))>1e-10;
jtheta(index_repeat)=[];
end