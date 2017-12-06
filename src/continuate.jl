## continuate.jl : implements continuation algorithms for the simulate function 
simulate(s::S) where S<:{Simulator} = throw(ArgumentError("simulate undefined for simulator of type $(typeof(s))"))
