function theta = qangle(q)
  theta = atan2(sqrt(q(2)^2 + q(3)^2 + q(4)^2), q(1));
 end
