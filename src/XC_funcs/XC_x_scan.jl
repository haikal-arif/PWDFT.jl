
# Strongly Constrained and Appropriately Normed Semilocal Density Functional
# Sun, J., Ruzsinzky, A., Perdew, J.P
# 10.1103/PhysRevLett.115.036402

function XC_x_scan(ρ, norm∇ρ, τ) # (ρ, |∇ρ|, τ)
  τunif = (0.3) * (3 * π^2)^(2 / 3) * ρ .^ (5 / 3)
  ∂τunif∂ρ = (0.5) * (3 * π^2)^(2 / 3) * ρ .^ (2 / 3)

  τW = norm∇ρ .^ 2 ./ (8 * ρ)
  ∂τW∂ρ = -norm∇ρ .^ 2 ./ (8 * ρ .^ 2)
  ∂τW∂norm∇ρ = 2 * norm∇ρ ./ (8 * ρ)

  α = (τ - τW) ./ τunif

  ∂α∂ρ = (-τunif .* ∂τW∂ρ - (τ - τW) .* ∂τunif∂ρ) ./ (τunif .* τunif)
  ∂α∂norm∇ρ = -(∂τW∂norm∇ρ) ./ τunif
  ∂α∂τ = 1 ./ τunif

  s = norm∇ρ ./ (2 * (3 * π^2)^(1 / 3) * ρ .^ (4 / 3))
  ∂s∂ρ = -2 * norm∇ρ ./ (3 * (3 * π * π)^(1 / 3) * ρ .^ (7 / 3))
  ∂s∂norm∇ρ = 1 ./ (2 * (3 * π * π)^(1 / 3) * ρ .^ (4 / 3))
  ∂s∂τ = 0

  μAK = 10 / 81
  b2 = sqrt(5913 / 405000)
  b1 = (511 / 13500) / (2 * b2)
  b3 = 0.5

  k1 = 0.065
  b4 = μAK^2 / k1 - (1606 / 18225) - b1^2



  x = μAK * s .^ 2 .* (1 .+ (b4 * (s .^ 2) / μAK) .* ℯ .^ (-abs(b4) * s .^ 2 / μAK)) +
      (b1 * s .^ 2 + b2 * (1 .- α) .* ℯ .^ (-b3 * (1 .- α) .^ 2)) .^ 2
  ∂x∂ρ = 2 * (b1 * s .^ 2 + b2 * (1 .- α) .* ℯ .^ (-b3 * (1 .- α) .* (1 .- α))) .* (2 * b1 * s .* ∂s∂ρ + b2 * (2 * b3 * (1 .- α) .^ 2 .- 1) .* ℯ .^ (-b3 * (1 .- α) .^ 2) .* ∂α∂ρ) +
         2 * (μAK .+ b4 * s .* s .* ℯ .^ (-b4 * s .* s / μAK) .* (2 .- b4 * s .* s / μAK)) .* s .* ∂s∂ρ
  ∂x∂norm∇ρ = 2 * (b1 * s .^ 2 + b2 * (1 .- α) .* ℯ .^ (-b3 * (1 .- α) .* (1 .- α))) .* (2 * b1 * s .* ∂s∂norm∇ρ + b2 * (2 * b3 * (1 .- α) .^ 2 .- 1) .* ℯ .^ (-b3 * (1 .- α) .^ 2) .* ∂α∂norm∇ρ) +
              2 * (μAK .+ b4 * s .* s .* ℯ .^ (-b4 * s .* s / μAK) .* (2 .- b4 * s .* s / μAK)) .* s .* ∂s∂norm∇ρ
  ∂x∂τ = 2 * (b1 * s .^ 2 + b2 * (1 .- α) .* ℯ .^ (-b3 * (1 .- α) .* (1 .- α))) .* (2 * b1 * s .* ∂s∂τ + b2 * (2 * b3 * (1 .- α) .^ 2 .- 1) .* ℯ .^ (-b3 * (1 .- α) .^ 2) .* ∂α∂τ) +
         2 * (μAK .+ b4 * s .* s .* ℯ .^ (-b4 * s .* s / μAK) .* (2 .- b4 * s .* s / μAK)) .* s .* ∂s∂τ

  h1x = 1 .+ k1 * (1 .- k1 ./ (k1 .+ x))
  ∂h1x∂ρ = 1 ./ ((1 .+ x / k1) .* (1 .+ x / k1)) .* ∂x∂ρ
  ∂h1x∂norm∇ρ = 1 ./ ((1 .+ x / k1) .* (1 .+ x / k1)) .* ∂x∂norm∇ρ
  ∂h1x∂τ = 1 ./ ((1 .+ x / k1) .* (1 .+ x / k1)) .* ∂x∂τ

  c1x = 0.667
  c2x = 0.8
  dx = 1.24

  fx = map(
    x -> (1 - x) > 0
         ? ℯ^(-c1x * x / (1 - x))
         : (-dx * ℯ^(c2x / (1 - x))),
    α)

  ∂fx∂ρ = map(
    (x, y, z) -> (1 - x) > 0
                 ? z * (-c1x / (1 - x)^2) * y
                 : z * (c2x / (1 - x)^2) * y,
    α, ∂α∂ρ, fx)

  ∂fx∂norm∇ρ = map(
    (x, y, z) -> (1 - x) > 0
                 ? z * (-c1x / (1 - x)^2) * y
                 : z * (c2x / (1 - x)^2) * y,
    α, ∂α∂norm∇ρ, fx)

  ∂fx∂τ = map(
    (x, y, z) -> (1 - x) > 0
                 ? z * (-c1x / (1 - x)^2) * y
                 : z * (c2x / (1 - x)^2) * y,
    α, ∂α∂τ, fx)

  h0x = 1.174
  a1 = 4.9479

  gx = 1 .- ℯ .^ (-a1 ./ (s .^ 0.5))
  ∂gx∂n = -0.5 * a1 * s .^ (-1.5) .* (1 .- gx) .* ∂s∂ρ
  ∂gx∂norm∇ρ = -0.5 * a1 * s .^ (-1.5) .* (1 .- gx) .* ∂s∂norm∇ρ
  ∂gx∂τ = -0.5 * a1 * s .^ (-1.5) .* (1 .- gx) .* ∂s∂τ

  Fxa = (h1x + fx .* (h0x .- h1x))
  ∂Fxa∂n = ∂fx∂ρ .* (h0x .- h1x) + (1 .- fx) .* ∂h1x∂ρ
  ∂Fxa∂norm∇ρ = ∂fx∂norm∇ρ .* (h0x .- h1x) + (1 .- fx) .* ∂h1x∂norm∇ρ
  ∂Fxa∂τ = ∂fx∂τ .* (h0x .- h1x) + (1 .- fx) .* ∂h1x∂τ

  Fx = Fxa .* gx
  ∂Fx∂n = Fxa .* ∂gx∂n + gx .* ∂Fxa∂n
  ∂Fx∂norm∇ρ = Fxa .* ∂gx∂norm∇ρ + gx .* ∂Fxa∂norm∇ρ
  ∂Fx∂τ = Fxa .* ∂gx∂τ + gx .* ∂Fxa∂τ

  nexunif = -(3 ./ (4 * π)) * (3 * π^2 * ρ) .^ (1 / 3)
  ∂nexunif∂n = -1 * (3 * ρ / π) .^ (1 / 3)

  sx = nexunif .* Fx
  v1x = (Fx .* ∂nexunif∂n) + (ρ .* nexunif .* ∂Fx∂n) # ∂e/∂n
  v2x = ρ .* nexunif .* ∂Fx∂norm∇ρ 
  v3x = ρ .* nexunif .* ∂Fx∂τ

  return sx, v1x, v2x, v3x
end

