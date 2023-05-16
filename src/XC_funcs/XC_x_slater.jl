function XC_x_slater(Rhoe::Float64)

    third = 1.0 / 3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    rs = pi34 / Rhoe^third

    f = -0.687247939924714
    alpha = 2.0 / 3.0

    ex = f * alpha / rs
    vx = 4.0 / 3.0 * f * alpha / rs
    return ex, vx
end

function XC_x_slater(Rhoe::Vector{Float64})
    ex = zeros(length(Rhoe))
    vx = zeros(length(Rhoe))
    for i in 1:length(Rhoe)
        ex[i], vx[i] = XC_x_slater(Rhoe[i])
    end

    return ex, vx
end
