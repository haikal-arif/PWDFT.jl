using ForwardDiff

function XC_c_scan(ρ::Float64, norm∇ρ::Float64, τ::Float64)
  # Strongly Constrained and Appropriately Normed Semilocal Density Functional
  # Sun, J., Ruzsinzky, A., Perdew, J.P
  # 10.1103/PhysRevLett.115.036402
  # This equation isn't included in the main paper
  # You can find it in its supplement
  #
  # TODO! we need to implement this according to paper
  # This implementation is under assumption that the system is non spin-polarized
  ζ = 0
  ϕ = 1
  γ = 0.031091
  τunif = (0.3) * (3 * π^2)^(2 / 3) * ρ^(5 / 3)
  ∂τunif∂ρ = (0.5) * (3 * π^2)^(2 / 3) * ρ^(2 / 3)

  τW = norm∇ρ^2 / (8 * ρ)
  ∂τW∂ρ = -norm∇ρ^2 / (8 * ρ^2)
  ∂τW∂norm∇ρ = 2 * norm∇ρ / (8 * ρ)

  α = (τ - τW) / τunif

  ∂α∂ρ = (-τunif * ∂τW∂ρ - (τ - τW) * ∂τunif∂ρ) / (τunif * τunif)
  ∂α∂norm∇ρ = -(∂τW∂norm∇ρ) / τunif
  ∂α∂τ = 1 / τunif

  s = norm∇ρ / (2 * (3 * π^2)^(1 / 3) * ρ^(4 / 3))
  ∂s∂ρ = -2 * norm∇ρ / (3 * (3 * π * π)^(1 / 3) * ρ^(7 / 3))
  ∂s∂norm∇ρ = 1 / (2 * (3 * π * π)^(1 / 3) * ρ^(4 / 3))
  ∂s∂τ = 0

  r = (3 / (4 * π * ρ))^(1 / 3) # steinitz radius
  t = (3 * π * π / 16)^(1 / 3) * s / (ϕ * sqrt(r))
  heaviside(t) = 0.5 * (sign(t) + 1)
  β = 0.066725 * (1 + 0.1 * r) / (1 + 0.1778 * r)
  c1c = 0.64
  dc = 0.7
  c2c = 1.5
  b3c = 0.125541
  b2c = 0.0889
  b1c = 0.0285764

  εLDA0c = -b1c / (1 + b2c * sqrt(r) + b3c * r)
  w0 = exp(-εLDA0c / b1c) - 1
  χ = 0.128026
  g = 1 / (1 + 4 * χ * s * s)^0.25
  H0 = b1c * log(1 + w0 * (1 - g))
  ε0c = (εLDA0c + H0)

  eLSDAc, vLSDAc = XC_c_lsda(ρ)
  w1 = exp(-eLSDAc / γ) - 1
  t = (3 * π * π / 16)^(1 / 3) * s / (ϕ * sqrt(r))
  A = β / (γ * w1)
  gAt2 = 1 / (1 + 4 * A * t * t)^(1 / 4)
  H1 = γ * ϕ * log(1 + w1 * (1 - gAt2))
  ε1c = eLSDAc + H1


  fc = exp(-c1c * α / (1 - α)) * heaviside(1 - α) - dc * exp(c2c / (1 - α)) * heaviside(α - 1)
  sc = (ε1c + fc * (ε0c - ε1c))
  v1c = 0
  v2c = 0
  v3c = 0
  return sc, v1c, v2c, v3c
end


function XC_c_scan(ρ::Vector{Float64}, norm∇ρ::Vector{Float64}, τ::Vector{Float64})
  sc = zeros(size(ρ))
  v1c = zeros(size(ρ))
  v2c = zeros(size(ρ))
  v3c = zeros(size(ρ))

  for index in 1:length(ρ)
    sc[index], v1c[index], v2c[index], v3c[index] = XC_c_scan(ρ[index], norm∇ρ[index], τ[index])
  end

  return sc, v1c, v2c, v3c
end
