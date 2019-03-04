function [committor]=compute_commit(K,end_points)

T = expm(K); % This line defines the Markovian transition matrx
N=size(T,1);
L = T - eye(N);
L(end_points,:) = 0;
for i = end_points
    L(i,i) = 1;
end
Z = zeros(N,1);
Z(end_points(2)) = 1;
committor = linsolve(L,Z);

end