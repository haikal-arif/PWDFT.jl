using ForwardDiff: derivative

# Strongly Constrained and Appropriately Normed Semilocal Density Functional
# Sun, J., Ruzsinzky, A., Perdew, J.P
# 10.1103/PhysRevLett.115.036402
function XC_x_scan(ρ::Float64, norm∇ρ::Float64, τ::Float64)
  τunif(n) = (0.3) * (3 * π^2)^(2 / 3) * n^(5 / 3)
  τW(n, norm∇n) = norm∇n^2 / 8 * n
  α(n, norm∇n, tau) = (tau - τW(n, norm∇n)) / τunif(n)
  s(n, norm∇n) = norm∇n / (2 * (3 * π^2)^(1 / 3) * n^(4 / 3))

  μAK = 10 / 81
  b2 = sqrt(5913 / 405000)
  b1 = (511 / 13500) / (2 * b2)
  b3 = 0.5

  k1 = 0.065
  b4 = μAK^2 / k1 - (1606 / 18225) - b1^2

  x(n, norm∇n, tau) = μAK * s(n, norm∇n)^2 * (1 + (b4 * (s(n, norm∇n)^2) / μAK) * exp(-abs(b4) * s(n, norm∇n)^2 / μAK)) + (b1 * s(n, norm∇n)^2 + b2 * (1 - α(n, norm∇n, τ)) * exp(-b3 * (1 - α(n, norm∇n, tau))^2))^2

  h1x(n, norm∇n, tau) = 1 + k1 * (1 - k1 / (k1 + x(n, norm∇n, tau)))

  c1x = 0.667
  c2x = 0.8
  dx = 1.24

  heaviside(t) = 0.5 * (sign(t) + 1)
  fx(n, norm∇n, tau) = exp(-c1x * α(n, norm∇n, tau) / (1 - α(n, norm∇n, tau))) * heaviside(1 - α(n, norm∇n, tau)) + (-dx * exp(c2x / (1 - α(n, norm∇n, tau)))) * (1 - heaviside(1 - α(n, norm∇n, tau)))

  h0x = 1.174
  a1 = 4.9479

  gx(n, norm∇n) = 1 - exp(-a1 / sqrt(s(n, norm∇n)))
  Fxa(n, norm∇n, tau) = h1x(n, norm∇n, tau) + fx(n, norm∇n, tau) * (h0x - h1x(n, norm∇n, tau))
  Fx(n, norm∇n, tau) = Fxa(n, norm∇n, tau) * gx(n, norm∇n)

  nexunif(n) = -(3 / (4 * π)) * (3 * π^2 * n)^(1 / 3)

  sx(n, norm∇n, tau) = nexunif(n) * Fx(n, norm∇n, tau)
  v1x(n, norm∇n, tau) = derivative(r -> sx(r, norm∇n, tau), n)
  v2x(n, norm∇n, tau) = derivative(r -> sx(n, r, tau), norm∇n)
  v3x(n, norm∇n, tau) = derivative(r -> sx(n, norm∇n, r), tau)

  return sx(ρ, norm∇ρ, τ), v1x(ρ, norm∇ρ, τ), v2x(ρ, norm∇ρ, τ) / norm∇ρ, v3x(ρ, norm∇ρ, τ)
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
