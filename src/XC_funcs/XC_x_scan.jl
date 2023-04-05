using ForwardDiff: derivative
using ForwardDiff

# Strongly Constrained and Appropriately Normed Semilocal Density Functional
# Sun, J., Ruzsinzky, A., Perdew, J.P
# 10.1103/PhysRevLett.115.036402
function XC_x_scan(ρ::Float64, norm∇ρ::Float64, τ::Float64)
  τunif(n) = (0.3) * (3 * π^2)^(2 / 3) * n^(5 / 3)
  τW(n, norm∇n) = norm∇n^2 / 8 * n

  α(n, norm∇n, tau) = (tau - τW(n, norm∇n)) / τunif(n)
  DαDρ(n, norm∇n, tau) = derivative(r -> α(r, norm∇n, tau), r)
  ∂α∂ρ = (norm∇ρ^2 - 5 * ρ * τ) / (2.7 * π * π * ρ^(11 / 3))
  ∂α∂norm∇ρ = -(5 * norm∇ρ) / (6 * (3 * π * π)^(2 / 3) * ρ^(8 / 3))
  ∂α∂τ = 1 / (0.3 * (3 * π * π)^(2 / 3) * ρ^(5 / 3))

  s = norm∇ρ / (2 * (3 * π^2)^(1 / 3) * ρ^(4 / 3))
  ∂s∂ρ = -2 * norm∇ρ / (3 * (3 * π * π)^(1 / 3) * ρ^(7 / 3))
  ∂s∂norm∇ρ = 1 / (2 * (3 * π * π)^(1 / 3) * ρ^(4 / 3))
  ∂s∂τ = 0

  μAK = 10 / 81
  b2 = sqrt(5913 / 405000)
  b1 = (511 / 13500) / (2 * b2)
  b3 = 0.5

  k1 = 0.065
  b4 = μAK^2 / k1 - (1606 / 18225) - b1^2

  x = μAK * s^2 * (1 + (b4 * (s^2) / μAK) * exp(-abs(b4) * s^2 / μAK)) + (b1 * s^2 + b2 * (1 - α) * exp(-b3 * (1 - α)^2))^2
  ∂x∂ρ = 2(b1 * s^2 + b2 * (1 - α) * exp(-b3 * (1 - α) * (1 - α))) * (2 * b1 * s * ∂s∂ρ + (2 * b3 * (1 - α)^2 - 1) * b2 * exp(-b3 * (1 - α)^2) * ∂α∂ρ) + 2(μAK + b4 * s * s * exp(-b4 * s * s / μAK) * (2 - b4 * s * s / μAK)) * s * ∂s∂ρ
  ∂x∂norm∇ρ = 2(b1 * s^2 + b2 * (1 - α) * exp(-b3 * (1 - α) * (1 - α))) * (2 * b1 * s * ∂s∂norm∇ρ + (2 * b3 * (1 - α)^2 - 1) * b2 * exp(-b3 * (1 - α)^2) * ∂α∂norm∇ρ) + 2(μAK + b4 * s * s * exp(-b4 * s * s / μAK) * (2 - b4 * s * s / μAK)) * s * ∂s∂norm∇ρ
  ∂x∂τ = 2(b1 * s^2 + b2 * (1 - α) * exp(-b3 * (1 - α) * (1 - α))) * (2 * b1 * s * ∂s∂τ + (2 * b3 * (1 - α)^2 - 1) * b2 * exp(-b3 * (1 - α)^2) * ∂α∂τ) + 2(μAK + b4 * s * s * exp(-b4 * s * s / μAK) * (2 - b4 * s * s / μAK)) * s * ∂s∂τ

  h1x = 1 + k1 * (1 - k1 / (k1 + x))
  ∂h1x∂ρ = 1 / ((1 + x / k1) * (1 + x / k1)) * ∂x∂ρ
  ∂h1x∂norm∇ρ = 1 / ((1 + x / k1) * (1 + x / k1)) * ∂x∂norm∇ρ
  ∂h1x∂τ = 1 / ((1 + x / k1) * (1 + x / k1)) * ∂x∂τ

  c1x = 0.667
  c2x = 0.8
  dx = 1.24

  fx = if ((1 - α) > 0)
    exp(-c1x * α / (1 - α))
  else
    (-dx * exp(c2x / (1 - α)))
  end

  ∂fx∂ρ = if ((1 - α) > 0)
    fx * (-c1x / (1 - α)) * ∂α∂ρ
  else
    fx * (c2x / (1 - α)^2) * ∂α∂ρ
  end

  ∂fx∂norm∇ρ = if ((1 - α) > 0)
    fx * (-c1x / (1 - α)) * ∂α∂norm∇ρ
  else
    fx * (c2x / (1 - α)^2) * ∂α∂norm∇ρ
  end

  ∂fx∂τ = if ((1 - α) > 0)
    fx * (-c1x / (1 - α)) * ∂α∂τ
  else
    fx * (c2x / (1 - α)^2) * ∂α∂τ
  end

  h0x = 1.174
  a1 = 4.9479

  gx = 1 - exp(-a1 / sqrt(s))
  ∂gx∂n = -0.5 * a1 * s^(-1.5) * (1 - gx) * ∂s∂ρ
  ∂gx∂norm∇ρ = -0.5 * a1 * s^(-1.5) * (1 - gx) * ∂s∂norm∇ρ
  ∂gx∂τ = -0.5 * a1 * s^(-1.5) * (1 - gx) * ∂s∂τ

  Fxa = (h1x + fx * (h0x - h1x))
  ∂Fxa∂n = ∂fx∂ρ * (h0x - h1x) + (1 - fx) * ∂h1x∂ρ
  ∂Fxa∂norm∇ρ = ∂fx∂norm∇ρ * (h0x - h1x) + (1 - fx) * ∂h1x∂norm∇ρ
  ∂Fxa∂τ = ∂fx∂τ * (h0x - h1x) + (1 - fx) * ∂h1x∂τ

  Fx = Fxa * gx
  ∂Fx∂n = Fxa * ∂gx∂n + gx * ∂Fxa∂n
  ∂Fx∂norm∇ρ = Fxa * ∂gx∂norm∇ρ + gx * ∂Fxa∂norm∇ρ
  ∂Fx∂τ = Fxa * ∂gx∂τ + gx * ∂Fxa∂τ

  nexunif = -(3 / (4 * π)) * (3 * π^2 * ρ)^(1 / 3)
  ∂nexunif∂n = -1 * (3 * ρ / π)^(1 / 3)

  sx = ρ * nexunif * Fx
  v1x = (Fx * ∂nexunif∂n) + (ρ * nexunif * (∂Fx∂n)) # ∂e/∂n
  v2x = nexunif * ∂Fx∂norm∇ρ / norm∇ρ # ∂e/∂n|∇ρ| * 1/|ρ|
  v3x = nexunif * ∂Fx∂τ

  return sx, v1x, v2x, v3x
end

function XC_x_scan(ρ::Vector{Float64}, norm∇ρ::Vector{Float64}, τ::Vector{Float64})
  sx = zeros(size(ρ))
  v1x = zeros(size(ρ))
  v2x = zeros(size(ρ))
  v3x = zeros(size(ρ))

  for index in 1:length(ρ)
    sx[index], v1x[index], v2x[index], v3x[index] = XC_x_scan(ρ[index], norm∇ρ[index], τ[index])
  end

  return sx, v1x, v2x, v3x
end
