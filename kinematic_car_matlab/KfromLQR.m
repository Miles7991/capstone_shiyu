clc;clear all;
A=zeros(4,4);
A(2,3)=1;
A(3,4)=1;
B=zeros(4,2);
B(1,1)=1;
B(4,2)=1;
K=-lqr(A,B,eye(4),eye(2))



