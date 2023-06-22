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
  ∂r∂ρ = -(1 / 3) * (3 / (4 * π))^(1 / 3) * ρ^(-4 / 3)
  ∂r∂normρ = 0
  ∂r∂τ = 0

  t = (3 * π * π / 16)^(1 / 3) * s / (ϕ * sqrt(r))
  ∂t∂ρ = (3 * π * π / 16)^(1 / 3) * (sqrt(r) * ∂s∂ρ - s * 0.5 * ∂r∂ρ / sqrt(r)) / r
  ∂t∂norm∇ρ = (3 * π * π / 16)^(1 / 3) * ∂s∂norm∇ρ / (ϕ * sqrt(r))
  ∂t∂τ = 0

  β = 0.066725 * (1 + 0.1 * r) / (1 + 0.1778 * r)
  ∂β∂ρ = 0.066725 * ((1 + 0.1778 * r) * 0.1 * ∂r∂ρ - (1 + 0.1 * r) * 0.1778 * ∂r∂ρ) / (1 + 0.1778 * r)^2
  ∂β∂normρ = 0
  ∂β∂τ = 0

  c1c = 0.64
  dc = 0.7
  c2c = 1.5
  b3c = 0.125541
  b2c = 0.0889
  b1c = 0.0285764

  εLDA0c = -b1c / (1 + b2c * sqrt(r) + b3c * r)
  vLDA0c = b1c * (1 + b2c * sqrt(r) + b3c * r)^(-2) * (0.5 * b2c / sqrt(r) + b3c) * ∂r∂ρ

  w0 = exp(-εLDA0c / b1c) - 1
  ∂w0∂ρ = -1 * vLDA0c * (w0 + 1) / b1c
  ∂w0∂normρ = 0
  ∂w0∂τ = 0

  χ = 0.128026
  g = 1 / (1 + 4 * χ * s * s)^0.25
  ∂g∂ρ = -2 * (1 + 4 * χ * s * s)^(-5 / 4) * χ * s * ∂s∂ρ
  ∂g∂normρ = -2 * (1 + 4 * χ * s * s)^(-5 / 4) * χ * s * ∂s∂norm∇ρ

  H0 = b1c * log(1 + w0 * (1 - g))
  ∂H0∂ρ = b1c * (∂w0∂ρ * (1 - g) - w0 * ∂g∂ρ) / (1 + w0 * (1 - g))
  ∂H0∂normρ = -b1c * w0 * ∂g∂normρ / (1 + w0 * (1 - g))

  ε0c = (εLDA0c + H0)
  ∂ε0c∂ρ = vLDA0c + ∂H0∂ρ
  ∂ε0c∂normρ = ∂H0∂normρ

  eLSDAc, vLSDAc = XC_c_lsda(ρ)

  w1 = exp(-eLSDAc / γ) - 1
  ∂w1∂ρ = -1 * vLSDAc * (w1 + 1) / γ
  ∂w1∂nomrρ = 0
  ∂w1∂τ = 0

  A = β / (γ * w1)
  ∂A∂ρ = (w1 * ∂β∂ρ - β * ∂w1∂ρ) / (γ * w1 * w1)
  ∂A∂normρ = 0
  ∂A∂τ = 0

  gAt2 = 1 / (1 + 4 * A * t * t)^(0.25)
  ∂gAt2∂ρ = -(0.25) * (1 + 4 * A * t * t)^(-5 / 4) * (4 * (A * 2 * t * ∂t∂ρ + ∂A∂ρ * t * t)) # ask why A doesn't derived
  ∂gAt2∂normρ = -2 * gAt2 * A * t * ∂t∂norm∇ρ / (1 + 4 * A * t^2)
  ∂gAt2∂τ = 0

  H1 = γ * ϕ^3 * log(1 + w1 * (1 - gAt2))
  ∂H1∂ρ = γ * ϕ^3 * (∂w1∂ρ * (1 - gAt2) - w1 * ∂gAt2∂ρ) / (1 + w1 * (1 - gAt2))
  ∂H1∂normρ = γ * ϕ^3 * (-w1 * ∂gAt2∂normρ) / (1 + w1 * (1 - gAt2))

  ε1c = eLSDAc + H1
  ∂ε1c∂ρ = vLSDAc + ∂H1∂ρ
  ∂ε1c∂normρ = ∂H1∂normρ

  fc = if (1 - α) > 0
    exp(-c1c * α / (1 - α))
  else
    -dc * exp(c2c / (1 - α))
  end

  ∂fc∂ρ = if (1 - α) > 0
    (-c1c * ∂α∂ρ / (1 - α)^2) * exp(-c1c * α / (1 - α))
  else
    -(c2c * ∂α∂ρ / (1 - α)^2) * dc * exp(c2c / (1 - α))
  end
  ∂fc∂normρ = if (1 - α) > 0
    (-c1c * ∂α∂norm∇ρ / (1 - α)^2) * exp(-c1c * α / (1 - α))
  else
    -(c2c * ∂α∂norm∇ρ / (1 - α)^2) * dc * exp(c2c / (1 - α))
  end
  ∂fc∂τ = if (1 - α) > 0
    (-c1c * ∂α∂τ / (1 - α)^2) * exp(-c1c * α / (1 - α))
  else
    -(c2c * ∂α∂τ / (1 - α)^2) * dc * exp(c2c / (1 - α))
  end


  sc = (ε1c + fc * (ε0c - ε1c))
  v1c = sc + (∂ε1c∂ρ + fc * (∂ε0c∂ρ - ∂ε1c∂ρ) + ∂fc∂ρ * (ε0c - ε1c)) * ρ
  v2c = (∂ε1c∂normρ + fc * (∂ε0c∂normρ - ∂ε1c∂normρ) + (ε0c - ε1c) * ∂fc∂normρ) * ρ
  v3c = ∂fc∂τ * (ε0c - ε1c) * ρ
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
