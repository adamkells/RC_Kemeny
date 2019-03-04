function [A]=make_graph(M)

n_edges = 0;
A = [0,1;1,0];
while n_edges < M
    a = randi(size(A,1),1);
    b = randi(size(A,1),1);
    if A(a,b)~=1
        if a~=b
            A(a,b) = 1;
            A(b,a) = 1;
            n_edges = n_edges + 2;
        end
    elseif A(a,b)==1
        r_var=rand(1);
        if r_var<0.5
            A(size(A,1)+1,b)=1;
            A(b,size(A,1))=1;
            %A(:,size(A,1))=zeros([size(A,1),1]);
        else
            A(a,size(A,1)+1)=1;
            A(size(A,2),a)=1;
            %A(size(A,2),:)=zeros([1,size(A,2)]);
        end
        n_edges = n_edges + 2;
    end
end