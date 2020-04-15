function q2 = update4(q1)

global dt mass massMat maxIter Pvec v ep nodes

q2 = q1;
q_free = q1(3:2*nodes-1);

error = 10 * ep;
iter = 0;

while (error > ep) && (iter < maxIter)

[Fs, Js] = stretchForce2(q2);
[Fb, Jb] = bendForce2(q2);

% first simulate all of the points as free to move
f = mass.*(q2-q1)/dt^2 - mass.*v/dt - (Fs+Fb+Pvec);
J = massMat/dt^2 - (Js+Jb);

%isolate the free to move portions
f_free = f(3: 2*nodes-1);
J_free = J(3:2*nodes-1, 3:2*nodes-1);

q_free = q_free - J_free\f_free;

%combine the free and fixed vectors
q2(1:2) = 0;
q2(2*nodes) = 0;
q2(3:2*nodes-1) = q_free;


error = abs(sum(f_free));
iter = iter + 1;
end

end