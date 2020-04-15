function [Fv, Jv] = velForce2(q)

    global q0 dt mu nodes middle R Rm
    
    Fv = -6*pi*mu*(q - q0) / dt;
    Jv = -6*pi*mu/dt*eye(length(q), length(q));
    
    for i=1:1:nodes
        if(i==middle)
            Fv(2*i-1:2*i) = Fv(2*i-1:2*i)* Rm;
                        
            Jv(2*i-1:2*i,2*i-1:2*i) = Jv(2*i-1:2*i,2*i-1:2*i) * Rm;
        else
            Fv(2*i-1:2*i) = Fv(2*i-1:2*i)* R;
                        
            Jv(2*i-1:2*i,2*i-1:2*i) = Jv(2*i-1:2*i,2*i-1:2*i) * R;
        end

    end
end