function A = A_q(q)
q_skew = [0 -q(3) q(2);
          q(3), 0, -q(1);
          -q(2), q(1), 0];
A = (q(4)^2-norm(q(1:3))^2)*eye(3) - 2*q(4)*q_skew + 2*q(1:3)*q(1:3)';
end