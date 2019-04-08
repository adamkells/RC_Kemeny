function [A]=erdosrenyi_2(nodes,N1,W)
%
% nodes: number of nodes in the 2 cluster graph
% N1: number of nodes in cluster 1 
% W: is the probability of inter/intracluster connections
%
p1=N1/nodes;
p2=1-p1;
p=[p1,p2];
x(1:nodes)=2;
x(1:N1)=1;
A=zeros(nodes);


for i=1:nodes
    for j=i:nodes
        P(i,j) = (4)*W(x(i),x(j))/(nodes*(p(x(i))*p(x(j))));
        r = rand(1);
        if r<P(i,j)
            A(i,j)=1;
        end
    end
end
for i=1:nodes
    for j=1:i
        A(i,j)=A(j,i);
    end
end

% this bit adds symmetry
%A=(A+A')/2;
%A=ceil(A);


rem=[];
for i=1:nodes
    if max(A(i,:))==0 || max(A(:,i))==0
        rem=[rem i];
    end
end
A(rem,:)=[];
A(:,rem)=[];
rem
end