## SPDE_data.jl : data struct for lognormal diffusion problem

## import statements ##
import Base.getindex

## user data ##
struct SPDE_Data{V}
    fields::V
end

getindex(s::SPDE_Data,index::Index) = s.fields[index]
