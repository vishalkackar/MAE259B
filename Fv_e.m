function Fv = Fv_e()

    global mu nodes middle R Rm v
    
    Fv = -6 * pi * mu * v;
    
    for i=1:1:nodes
        if(i == middle)
            Fv(2*i-1:2*i) = Fv(2*i-1:2*i)* Rm;
                        
        else
            Fv(2*i-1:2*i) = Fv(2*i-1:2*i)* R;
                        
        end

    end
end