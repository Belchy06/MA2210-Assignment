function V = gramSchmidt_errorcheck(U)

% Check that the input is sensible.
if ischar(U)
	error('gramSchmidt cannot work with strings!')
end

N = size(U,1);
V = zeros(N,N);

for i = 1:N
	V(:,i) = U(:,i);
	
	for j = 1:i-1
		V(:,i) = V(:,i) - (V(:,i)'*V(:,j)) * V(:,j);
	end

	% Check for lack of linear independence:
	if norm(V(:,i)) < 1e-10
		error('Vectors in U are not independent!')
	end	

	V(:,i) = V(:,i) / norm(V(:,i));
end

end
