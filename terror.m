function [K,Adj,corr] = terror()

load('examples/terrorist.mat')
Adj = A;
corr = correct_clustering;

nodes = size(A,1);

for i = 1:nodes
    for j = 1:nodes
        K(i,j)=Adj(i,j)/sum(Adj(:,j));
    end
end
for i = 1:nodes
    K(i,i)=0;
    K(i,i)=-sum(K(i,:));
end

end