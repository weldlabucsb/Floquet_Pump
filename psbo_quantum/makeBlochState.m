function [wf] = makeBlochState(v)
    wf=@(x) v(1)*sqrt(1/pi); 
    for nn=2:2:(length(v)-1)
       wf=@(x) wf(x) + ...
           v(nn)*sqrt(2/pi)*sin(nn.*x)+...
           v(nn+1)*sqrt(2/pi)*cos(nn.*x);    
    end
end

