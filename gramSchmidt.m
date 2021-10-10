function V = gramSchmidt(U)

N = size(U,1);

% Initialise with an empty set of vectors.
V = zeros(N,N);
% Run over all of the vectors. For each i we find the properly orthogonalised vector v.
for i = 1:N
	% v starts off being equal to u
	V(:,i) = U(:,i);
	% Now we run over all vectors j which have already been orthogonalised. These are vectors with indices from 1 up to (but not including) i.
	for j = 1:i-1
		% Subtract off the component in the direction of the v_j vector.
		V(:,i) = V(:,i) - (V(:,i)'*V(:,j)) * V(:,j);
	end
	% Finally, we have an orthogonal vector. We then normalise it (necessary for the formula above).
	V(:,i) = V(:,i) / norm(V(:,i));
end

end
