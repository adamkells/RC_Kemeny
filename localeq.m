function R = localeq(K,eq,A)

eqmat=zeros(length(eq));
for i=1:length(eq)
   eqmat(:,i)=eq; 
end
R=(A'*(K.*eqmat)*A)./(eq'*A);

end