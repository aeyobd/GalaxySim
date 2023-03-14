module GalaxySim
export testrun

# Write your package code here.
include("physics.jl")

function testrun()
    Physics.main()
end

end
