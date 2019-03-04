function [A]=make_lin_graph(M)

A=zeros(M);
for i=1:M-1
    A(i,i+1)=1;
    A(i+1,i)=1;
end