function q = qrotate(v,q)    
    q = qmult(qmult(q,v),qconj(q));
end