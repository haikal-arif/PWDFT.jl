# A wrapper struct for Nspin wavefunction with Gamma-point trick
mutable struct BlochWavefuncGamma
    data::Vector{Array{ComplexF64,2}}
end

# Constructor, from an array to BlochWavefuncGamma (Nspin=1)
function BlochWavefuncGamma(psi::Array{ComplexF64,2})
    return BlochWavefuncGamma([psi])
end


import Base: getindex
function getindex(X::BlochWavefuncGamma, i)
    return X.data[i]
end

import Base: setindex!
function setindex!(X::BlochWavefuncGamma, v::Array{ComplexF64}, i)
    X.data[i] = v
    return
end


# Also these?
# firstindex(X)   The first index, used in X[begin]
# lastindex(X)    The last index, used in X[end]


import LinearAlgebra: dot
function dot(v1::BlochWavefuncGamma, v2::BlochWavefuncGamma)
    c = dot(v1.data, v2.data)
    s = c + conj(c)
    Nspin = length(v1)
    Nstates = size(v1.data[1], 2)
    ds = 0.0 + im * 0.0
    for ispin in 1:Nspin, ist in 1:Nstates
        ds = ds + conj(v1.data[ispin][1, ist]) * v2.data[ispin][1, ist]
    end
    return s - ds
end


function dot_gamma(v1::Array{ComplexF64}, v2::Array{ComplexF64})
    c = dot(v1, v2)
    s = c + conj(c)
    Nstates = size(v1, 2)
    ds = 0.0 + im * 0.0
    for ist in 1:Nstates
        ds = ds + conj(v1[1, ist]) * v2[1, ist]
    end
    return s - ds
end


import Base: length
function length(v::BlochWavefuncGamma)
    return length(v.data)
end

#
# Operators
#

import Base: +
function +(v1::BlochWavefuncGamma, v2::BlochWavefuncGamma)
    return BlochWavefuncGamma(v1.data + v2.data)
end

import Base: -
function -(v1::BlochWavefuncGamma, v2::BlochWavefuncGamma)
    return BlochWavefuncGamma(v1.data - v2.data)
end



function zeros_BlochWavefuncGamma(Ham::HamiltonianGamma)
    return zeros_BlochWavefuncGamma(Ham.pw.gvecw.Ngw, Ham.electrons.Nstates, Nspin=Ham.electrons.Nspin)
end


function zeros_BlochWavefuncGamma(Nbasis::Int64, Nstates::Int64; Nspin=1)
    data = Vector{Array{ComplexF64,2}}(undef, Nspin)
    for ispin in 1:Nspin
        data[ispin] = zeros(ComplexF64, Nbasis, Nstates)
    end
    return BlochWavefuncGamma(data)
end




# random BlochWavefuncGamma from randn

function randn_BlochWavefuncGamma(Ham::HamiltonianGamma; zero_dc=false)
    return randn_BlochWavefuncGamma(
        Ham.pw.gvecw.Ngw, Ham.electrons.Nstates,
        Nspin=Ham.electrons.Nspin, zero_dc=zero_dc
    )
end

function randn_BlochWavefuncGamma(Nbasis::Int64, Nstates::Int64; Nspin=1, zero_dc=false)
    #
    data = Vector{Array{ComplexF64,2}}(undef, Nspin)
    for ispin in 1:Nspin
        data[ispin] = randn(ComplexF64, Nbasis, Nstates)
        if zero_dc
            # Set DC component to zero
            for ist in 1:Nstates
                data[ispin][1, ist] = 0.0
            end
        else
            # Set DC component (ig=1) to real number
            for ist in 1:Nstates
                data[ispin][1, ist] = data[ispin][1, ist] + conj(data[ispin][1, ist])
            end
        end
        ortho_sqrt_gamma!(data[ispin])
    end

    return BlochWavefuncGamma(data)
end

# 
function rand_BlochWavefuncGamma(Ham::HamiltonianGamma; zero_dc=false)
    return rand_BlochWavefuncGamma(
        Ham.pw.gvecw.Ngw, Ham.electrons.Nstates,
        Nspin=Ham.electrons.Nspin, zero_dc=zero_dc
    )
end

function rand_BlochWavefuncGamma(Nbasis::Int64, Nstates::Int64; Nspin=1, zero_dc=false)
    #
    data = Vector{Array{ComplexF64,2}}(undef, Nspin)
    for ispin in 1:Nspin
        data[ispin] = rand(ComplexF64, Nbasis, Nstates)
        if zero_dc
            # Set DC component to zero
            for ist in 1:Nstates
                data[ispin][1, ist] = 0.0
            end
        else
            # Set DC component (ig=1) to real number
            for ist in 1:Nstates
                data[ispin][1, ist] = data[ispin][1, ist] + conj(data[ispin][1, ist])
            end
        end
        ortho_sqrt_gamma!(data[ispin])
    end

    return BlochWavefuncGamma(data)
end



function ortho_check(psis::BlochWavefuncGamma)
    Nspin = length(psis)
    for ispin in 1:Nspin
        ortho_check_gamma(psis.data[ispin])
    end
    return
end

function ortho_check_gamma(psis::BlochWavefuncGamma)
    Nspin = length(psis)
    for ispin in 1:Nspin
        ortho_check_gamma(psis.data[ispin])
    end
    return
end

function ortho_check_gamma(psi::Array{ComplexF64,2})
    Nstates = size(psi)[2]
    @printf("\nNorm check:\n")
    for ist = 1:Nstates
        c = dot(psi[:, ist], psi[:, ist])
        c = c + conj(c) - conj(psi[1, ist]) * psi[1, ist]
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\nOrtho check w.r.t state #1:\n")
    for ist = 2:Nstates
        c = dot(psi[:, ist], psi[:, 1])
        c = c + conj(c) - conj(psi[1, ist]) * psi[1, 1]
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\n")
    return
end

function overlap_gamma(psi1::Array{ComplexF64,2}, psi2::Array{ComplexF64,2})
    C = psi2' * psi1
    Nstates = size(psi1, 2)
    v1g = zeros(ComplexF64, Nstates)
    v2g = zeros(ComplexF64, Nstates)
    for ist in 1:Nstates
        v1g[ist] = psi2[1, ist] # psi is the array that is conj transposed in the orig expression
        v2g[ist] = psi1[1, ist]  # v2g is the array that will be conj transposed
    end
    return C + conj(C) - v1g * v2g'
end

function ortho_sqrt_gamma!(psis::BlochWavefuncGamma)
    Nspin = size(psis, 1)
    Nstates = size(psis.data, 2)
    C = zeros(ComplexF64, Nstates, Nstates)
    for i in 1:Nspin
        C = overlap_gamma(psis.data[i], psis.data[i])
        @views psis.data[i][:, :] = psis.data[i][:, :] * inv(sqrt(C))
    end
    return
end

function ortho_sqrt_gamma!(psi::Array{ComplexF64,2})
    C = overlap_gamma(psi, psi)
    @views psi[:, :] = psi[:, :] * inv(sqrt(C))
    return
end


function ortho_gram_schmidt!(psis::BlochWavefuncGamma)
    Nspin = length(psis)
    for ispin in 1:Nspin
        ortho_GS_gamma!(psis.data[ispin])
    end
    return
end

function ortho_GS_gamma(vin::Array{ComplexF64,2})
    Ncol = size(vin)[2]
    v = copy(vin)
    for ii = 1:Ncol
        zz = 2 * dot(v[:, ii], v[:, ii]) - conj(v[1, ii]) * v[1, ii]
        v[:, ii] = v[:, ii] / sqrt(zz)
        for jj = ii+1:Ncol
            puv = _prj_gamma(v[:, ii], v[:, jj])
            v[:, jj] = v[:, jj] - puv * v[:, ii]
        end
    end
    return v
end

function ortho_GS_gamma!(v::Array{ComplexF64,2})
    Ncol = size(v)[2]
    for ii = 1:Ncol
        c = dot(v[:, ii], v[:, ii])
        zz = c + conj(c) - conj(v[1, ii]) * v[1, ii]
        v[:, ii] = v[:, ii] / sqrt(zz)
        for jj = ii+1:Ncol
            puv = _prj_gamma(v[:, ii], v[:, jj])
            v[:, jj] = v[:, jj] - puv * v[:, ii]
        end
    end
    return
end

function _prj_gamma(u, v)
    c = dot(u, v)
    num = c + conj(c) - conj(u[1]) * v[1]
    c = dot(u, u)
    denum = c + conj(c) - conj(u[1]) * u[1]
    return num / denum
end