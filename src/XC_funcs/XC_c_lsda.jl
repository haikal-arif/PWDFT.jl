# https://sci-hub.se/10.1103/PhysRevB.45.13244
# Perdew, J. P., & Wang, Y. (1992). 
# Accurate and simple analytic representation of 
# the electron-gas correlation energy. 
# Physical Review B, 45(23), 13244–13249. 
function XC_c_lsda(n::Float64)
  # We assume there is no spin polarization so we only have εc(rs, 0) term

  rs = (3 / (4 * π * n))^(1 / 3)
  DrsDn = (1 / 3) * rs / n

  A = 0.031091
  alpha = 0.21370
  β1 = 7.5957
  β2 = 3.5876
  β3 = 1.6382
  β4 = 0.49294

  Q0 = -2 * A * (1 + alpha * rs)
  Q1 = 2 * A * (β1 * rs^0.5 + β2 * rs + β3 * rs^1.5 + β4 * rs * rs)
  Q1prime = A * (β1 * rs^-0.5 + 2 * β2 + 3 * β3 * rs^0.5 + 4 * β4 * rs)

  ec = Q0 * log(1 + 1 / Q1)
  vc = (-2 * A * alpha * log(1 + 1 / Q1) - (Q0 * Q1prime) / (Q1 * Q1 + Q1)) * DrsDn


  return ec, vc
end
