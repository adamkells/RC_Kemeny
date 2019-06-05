A=zeros(length(unique(data)));
val=1;
x=unique(data);
for i=1:length(x)
    100*i/length(x)
    if max(max(ismember(data,x(i))))==1
        data(data==x(i))=val;
        val=val+1;
    end
end

for i=1:size(data,1)
    A(data(i,1),data(i,2))=1;
end
A=sparse(A);

% for i=1:length(A)
%     if A(i,i)==1
%         if sum(A(i,:))==1
%             A(i,:)=[];
%             A(:,i)=[];
%         end
%     end
% end