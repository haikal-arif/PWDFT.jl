# Strongly Constrained and Appropriately Normed Semilocal Density Functional
# Sun, J., Ruzsinzky, A., Perdew, J.P
# 10.1103/PhysRevLett.115.036402
function XC_x_scan(rho, grho, psi, τ)
    τunif = (0.3) * (3 * π^2)^(2 / 3) * rho^(5 / 3)
    τW = grho^2 / 8 * rho
    α = (τ - τW) / τunif
    s = abs(grho) / (2 * (3 * π^2)^(1 / 3) * rho^(4 / 3))
    μAK = 10 / 81
    b2 = sqrt(5913 / 405000)
    b1 = (511 / 13500) / (2 * b2)
    b3 = 0.5

    k1 = 0.065
    b4 = μAK^2 / k1 - (1606 / 18225) - b1^2
    x = μAK * s^2 * (1 + (b4 * (s^2) / μAK) * exp(-abs(b4) * s^2 / μAK)) + (b1 * s^2 + b2 * (1 - α) * exp(-b3 * (1 - α)^2))^2
    h1x = 1 + k1 - k1 / (1 + x / k1)

    c1x = 0.667
    c2x = 0.8
    dx = 1.24
    fx = ((1 - α) > 0) ? exp(-c1x * α / (1 - α)) : (-dx * exp(c2x / (1 - α)))

    h0x = 1.174
    a1 = 4.9479
    gx = 1 - exp(-a1 / sqrt(s))
    Fx = (h1x + fx * (h0x - h1x)) * gx
    exunif = -(3 / (4 * π)) * (3 * π^2 * rho)^(1 / 3)
    sx = rho * exunif * Fx
    return sx
end
