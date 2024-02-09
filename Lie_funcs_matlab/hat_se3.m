% Auxiliary "hat" (lift) function, from the vector space to the Lie algebra

function V = hat_se3(xi)
V = [hat_so3(xi(4:6)), xi(1:3); 0 0 0 0];
