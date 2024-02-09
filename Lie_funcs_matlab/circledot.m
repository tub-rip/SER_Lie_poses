function M = circledot(p)
M = [p(4)*eye(3), -hat_so3(p(1:3)); zeros(1,6)];