function [Fb, Jb] = bendForce2(q)
    
    global EI dL nodes
    
    Fb = zeros(numel(q),1);
    Jb = zeros(numel(q),numel(q));
    
    for i=1:1:nodes-2
        x1 = q(2*i-1); y1 = q(2*i);
        x2 = q(2*i+1); y2 = q(2*i+2); 
        x3 = q(2*i+3); y3 = q(2*i+4);
        
        g = gradEb(x1, y1, x2, y2, x3, y3, dL, EI);
        h = hessEb(x1, y1, x2, y2, x3, y3, dL, EI);
        
        Fb(2*i-1:2*i+4) = Fb(2*i-1:2*i+4) - g;
        Jb(2*i-1:2*i+4,2*i-1:2*i+4) = Jb(2*i-1:2*i+4,2*i-1:2*i+4) - h;
    end 
    
end