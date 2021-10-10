clear all
clc

%% Setup
% Etching circuit boards. Do tests on different parameters which result in the variables:
% quality, thickness, duration, acid concentration, temperature
% The study looks at how the quality of the board is related to the other variables.
%
% Perhaps the analysis results in an eigenvalue problem for a matrix of the form:

A = [0.75 0 1.75 0 ; -0.25 1 1.75 0;5.25 0 4.25 0;0 0 0 -5]

[U,D] = eig(A);
evals = diag(D)

U

%%
% Consider large eigenvalues to be the significant ones (either large +ve or large -ve).
% Rank these by descending absolute order:

[dummy,ind] = sort(abs(evals), 'descend')
evals = evals(ind);
U = U(:,ind);

evals
U

%%
%
% evals =
%      6.0000   -5.0000    1.0000   -1.0000
%
% U =
%     -0.3015         0         0   -0.5774   <-- thickness
%     -0.3015         0    1.0000   -0.5774   <-- duration
%     -0.9045         0         0    0.5774   <-- acid concentration
%           0    1.0000         0         0   <-- temperature
%
% We see there is a determental contribution from temperature (-6 and evec is purely in direction of T).
% The largest contribution is a positive effect of a combination of thickness, duration and acid concentration.
% Similarly, there is a positive effect from the duration alone.
% Finally, there is a small negative effect from a different combination of thickness, duration, and concentration.

% We multiply some eigenvectors by a negative as we are focusing on positive variables.
% Note: this is purely cosmetic.
U(:,[1,4]) = -U(:,[1,4]);

U

%%
% Although these vectors tell us lots of information, they are not orthogonal:

disp('Overlap of eigenvectors 1 and 3')
U(:,1)' * U(:,3)

%%
disp('Overlap of all eigenvectors as a matrix')
U'*U

%%
% So we can't use them as a coordinate system. For example, if we wanted to fix 1 lot of the 3rd vector (for its positive contribution), we would ask:

vec = 1*U(:,3)

%%
% This has set the duration. And if we added 2 lots of the 1st vector (for its positive contribution) we would find:

vec = vec + 2*U(:,1)

% But this new vec has changed the duration required! Of course, this is what we want practically, but it means that the original statement of "1 lot of the 3rd eigenvector" is no longer true. Hence, we have a poor description of the vector.

%% Constructing an orthogonal space
% So we should come up with a new basis that is better suited for representing vectors.
% The Gram-Schmidt procedure does this. It works by taking each vector in turn and making it orthogonal to all other vector already considered.
% We will call this orthogonal set "V"

% The first vector is just the same as the original:

V(:,1) = U(:,1)

%%
% The second vector, in this case, happens to already be orthogonal to the first, so we take it as is:

disp('Overlap:')
V(:,1)' * U(:,2)

V(:,2) = U(:,2)

%%
% The third vector does have overlap with the first. So we simply subtract this away:

V(:,3) = U(:,3) - (U(:,3)'*V(:,1)) * V(:,1);
V(:,3) = V(:,3) / norm(V(:,3));

V

%%
% This looks like we have taken a simple eigenvector [0 1 0] and made it
% more complicated. This is motivated by the eigenvalues which indicate
% quality - in words "Duration alone does not have a significant effect on
% quality, but the combination of thickness, acid concentration AND
% duration does".

% V is orthogonal so far:

V'*V

%%
% Finally, the 4th vector has some overlap with the 1st and 3rd which must
% be removed:

V(:,4) = U(:,4) - (U(:,4)'*V(:,1)) * V(:,1) - (U(:,4)'*V(:,3)) * V(:,3);
V(:,4) = V(:,4) / norm(V(:,4));

V

%%
% V is fully orthogonal
V'*V

% Note that the original ordering of the vectors changes the outcome of the
% Gram-Schmidt procedure. In this case, the ordering is chosen to favour
% the more significant eigenvectors.


%% Loop constructs

% Although the above is simple to code up for a 4x4 matrix, it is something
% we would like to automate for arbitrarily large matrices. Consider of two
% diagonals filled with ones.

N = 10;
U = eye(N) + diag(ones(N-1,1), -1)

%%
% We first would like to normalise all of these vectors. If we were to
% normalise the first three vectors, we could write:

U(:,1) = U(:,1) / norm(U(:,1));
U(:,2) = U(:,2) / norm(U(:,2));
U(:,3) = U(:,3) / norm(U(:,3));

U

%%

% These types of sequences, where a single command is repeated many times
% with only one number changing, are simple to write using a "for loop".

for i = 1:10
	U(:,i) = U(:,i) / norm(U(:,i));
end

%%
% The for loop will run its "body" for each value of i that is given. It is similar to writing:
i = 1;
U(:,i) = U(:,i) / norm(U(:,i));
i = 2;
U(:,i) = U(:,i) / norm(U(:,i));
i = 3;
U(:,i) = U(:,i) / norm(U(:,i));
% All the way to i=10.

%%
% Of course, we can use a variable instead of a fixed limit:
for i = 1:N
	U(:,i) = U(:,i) / norm(U(:,i));
end

U

%%
% The vectors in U are not orthogonal:
U'*U

%%
% The orthogonalised vectors of V can be built up using a for loop. We will
% start with just copying them:
V = zeros(N,N);
for i = 1:N
	V(:,i) = U(:,i);
end
V
%%
% Next we will subtract any overlap with the first vector from the 2nd
% vectors onwards:
for i = 2:N
	V(:,i) = V(:,i) - (V(:,i)'*V(:,1)) * V(:,1);
end

V

%%
% This has sorted out the overlap with the first vector - all other vectors
% are now orthogonal to the first vector.
V'*V

%%
% We now need to repeat this process with the 2nd, 3rd, ... vectors to
% orthogonalise the remaining vectors. This can also be done with a for
% loop. We present the entire segment of code with comments interspersed:

% Initialise with an empty set of vectors.
V = zeros(N,N);
% Run over all of the vectors. For each i we find the properly
% orthogonalised vector v.
for i = 1:N
	% v starts off being equal to u
	V(:,i) = U(:,i);
	% Now we run over all vectors j which have already been orthogonalised.
	% These are vectors with indices from 1 up to (but not including) i.
	for j = 1:i-1
		% Subtract off the component in the direction of the v_j vector.
		V(:,i) = V(:,i) - (V(:,i)'*V(:,j)) * V(:,j);
	end
	% Finally, we have an orthogonal vector. We then normalise it
	% (necessary for the formula above).
	V(:,i) = V(:,i) / norm(V(:,i));
end

V

%%
% There are no overlaps now
V'*V

%%
% The result has some interesting structure. If we put the above
% Gram-Schmidt code into a function we can easily see how this structure
% changes with increasing matrix size.

N = 5;
U = eye(N) + diag(ones(N-1,1), -1);

V5 = gramSchmidt(U)

N = 100;
U = eye(N) + diag(ones(N-1,1), -1);

V100 = gramSchmidt(U);

% The first 4 vectors are the same for this matrix:

V5
V100(1:5,1:5)

%%
% We can also use for loops to test out a set of different parameters:

for N = [5,10,50,100,3,2]
	U = eye(N) + diag(ones(N-1,1), -1);
	V = gramSchmidt(U);

	% Here we use the "string concatenation" ability of [].
	disp(['For N = ', num2str(N), ' we find the top right element is ', num2str(V(1,end))])
end

% This demonstrates that the for loop can take any set of values, not just
% values that are in sequence. They do not even have to be in order.

%% Conditional statements
% Finally, we often need to "branch" based on a parameter. That is, do one
% of two things depending on a value. For example:

val = -0.1
if val < 0
	disp('Negative')
else
	disp('Positive')
end

% Here the "if" statement checks whether what it is given is true. If it is
% true, it executes the next set of commands. If it is instead false, it
% will move to an "else" block (which does not have to be included) and
% execute those commands.

%%
% A more practical example is error checking. If we gave the gramSchmidt
% procedure a string, it should report a meaningful error. Currently it
% doesn't:

gramSchmidt('asdf')

%%
% Updating the function with some if statements solves that problem.
gramSchmidt_errorcheck(rand(2))
gramSchmidt_errorcheck('asdf')

%%
% The gramSchmidt procedure would also not work if the vectors were
% linearly dependent. We can see that this is sometimes missed:

V = gramSchmidt(ones(2))

V'*V

%%
% Obviously these are the same vector! The issue is that one vector became
% very close to zero. We can check this with an if statement:

V = gramSchmidt_errorcheck(ones(2))


