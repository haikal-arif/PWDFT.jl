import Base: show
using Printf

function show(io::IO, pspotNL::PsPotNL)
    @printf(io, "\n")
    @printf(io, "NbetaNL = %d", pspotNL.NbetaNL)
end
show(pspotNL::PsPotNL) = show(stdout, pspotNL)