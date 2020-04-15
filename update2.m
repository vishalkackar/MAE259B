function q2 = update2(q1)

    global dt ep maxIter Weight mass v massMat
    
    q2 = q1;
    
    error = 10 * ep;
    iter = 1;
    while ((error > ep) && (iter <= maxIter))
        [Fv, Jv] = velForce2(q2);
        [Fs, Js] = stretchForce2(q2);
        [Fb, Jb] = bendForce2(q2);
        
                
        f = mass.*(q2-q1)/dt^2 - mass.*v/dt - (Fv+Fs+Fb+Weight);
        
        J = massMat/dt^2 - (Jv+Js+Jb);
        
        
        q2 = q2 - J\f;
        error = abs(sum(f));
        
        iter = iter + 1;
    end
    
    
end