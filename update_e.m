function q2 = update_e(q1)
    
    global dt mass Weight v

    Fv = Fv_e();
    Fs = Fs_e(q1);
    Fb = Fb_e(q1);
    
    F = Fv+Fs+Fb+Weight;
    
    q2 = q1 + dt^2 ./ mass .*(mass.*v/dt + F);
    

end