function [Fs, Js] = stretchForce2(q)

    global EA dL nodes
        
    Fs = zeros(length(q),1);
    Js = zeros(length(q), length(q));
    
    for i=1:1:nodes-1
        x1 = q(2*i-1); y1 = q(2*i);
        x2 = q(2*i+1); y2 = q(2*i+2);
        
        g = gradEs(x1,y1,x2,y2,dL,EA);
        h = hessEs(x1,y1,x2,y2,dL,EA);
        
        Fs(2*i-1:2*i+2) = Fs(2*i-1:2*i+2) - g;
        Js(2*i-1:2*i+2,2*i-1:2*i+2) = Js(2*i-1:2*i+2,2*i-1:2*i+2) - h;
    end

    
end