simulate(model::HPoisson, region::Area) = SPP(simulate(model.μ, region), region)

estimate(::HPoisson, spp::SPP) = HPoisson(length(spp) / measure(spp.region))

function intensity(model::HPoisson, spp::SPP)
    cif(x) = model.μ
    ∫cif(A) = model.μ * measure(A)
    return DomainFunction(cif, spp.region, integral=∫cif)
end