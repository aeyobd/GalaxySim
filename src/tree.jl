module Tree

using ..Particles
using ..Constants
using ..Density
using LinearAlgebra
using Printf

export make_tree, dv_G!, find_neighbors!, find_within_r, U_G


const max_depth = 30 # Adjust the maximum depth of the octree according to your problem


mutable struct OctreeNode
    children::Union{Nothing, Vector{OctreeNode}}
    particles::Vector{Particle}
    center::Vector{F}
    center_of_mass::Vector{F}

    mass::F
    size::F
    depth::Int
end


function create_octree_node(center, size, depth=0)
    return OctreeNode(nothing, [], center, zeros(3), 0, size, depth)
end


function create_children(node::OctreeNode)
    child_size = node.size / 2
    offsets = [
        [+1, +1, +1],
        [+1, +1, -1],
        [+1, -1, +1],
        [+1, -1, -1],
        [-1, +1, +1],
        [-1, +1, -1],
        [-1, -1, +1],
        [-1, -1, -1],
    ]
    return [create_octree_node(node.center .+ offset .* (child_size / 2), 
                               child_size,
                              node.depth + 1) 
            for offset in offsets]
end


function insert_particle!(node::OctreeNode, particle::Particle, depth::Int)
    node.center_of_mass = (node.mass*node.center_of_mass  + particle.x*particle.m)/(node.mass + particle.m)
    node.mass += particle.m

    # If either at maximum depth or no children, we can insert the particle
    if depth >= max_depth || node.children === nothing
        push!(node.particles, particle)
        return
    end

    if node.children === nothing
        node.children = create_children(node)

        for p in node.particles
            insert_particle!(node, p, depth + 1)
        end
        node.particles = []
    end

    index = 1
    index += 4 * (particle.x[1] < node.center[1])
    index += 2 * (particle.x[2] < node.center[2])
    index += 1 * (particle.x[3] < node.center[3])

    insert_particle!(node.children[index], particle, depth + 1)
end



"""
finds all particles within R of the target particle

If nothing is close, returns the nearest
"""
function find_within_r(p::Particle, node::OctreeNode, r::F)
    result = Particle[]
    min_distance = Inf
    min_particle = nothing

    function traverse(node::OctreeNode, r::F)
        if node.children === nothing
            for q in node.particles
                if q == p
                    continue
                end

                d = norm(q.x - p.x)
                if d <= r 
                    push!(result, q)
                    min_distance = -1 # don't need this anymore
                elseif min_distance > 0 && d <= min_distance
                    min_distance = d
                    min_particle = q
                end
            end
            return
        end

        for child in node.children
            d = norm(child.center - p.x) - child.size
            if  d <= r || (min_distance > 0 && d <= min_distance) || min_particle===nothing
                traverse(child, r)
            end
        end
    end

    traverse(node, r)

    if min_distance > 0
        if min_particle == nothing
            println(p)
            throw(error("could not find nearest neighbor, are there NaNs?"))
        end
        return [min_particle]
    end

    return result
end



function find_neighbors!(p::Particle, tree, params)
    nearby = find_within_r(p, tree, p.h)

    p.neighbors = NParticle[]
    for q in nearby
        q1 = interpolate(q, p.t, params)
        q1.w = W(p, q1)
        q1.dw = ∇W(p, q1)
        push!(p.neighbors, q1)
    end
    return p.neighbors
end


"""
Gravitational acceleration between two masses
"""
function dv_G!(p::Particle, q::Particle, params)
    eps = min(p.h, q.h) * params.r_plummer 
    r = q.x .- p.x
    p.dv_G .+= G*q.m * normalize(r)/(sum(r.*r) + eps^2)
    return p.dv_G
end


function dv_G!(p::Particle, node::OctreeNode, params)
    if node.children === nothing
        if node.particles != [] 
            for q in node.particles
                if q != p
                    dv_G!(p, q, params)
                end
            end
        end
        return p.dv_G
    end

    distance = norm(node.center - p.x)
    if node.size / distance < params.theta
        # Node is sufficiently far away, use a single force calculation for the node
        # Don't use softening here
        p.dv_G .+= G * node.mass * (node.center_of_mass - p.x) / distance^3
    else
        # Node is not sufficiently far away, recurse to children
        for child in node.children
            dv_G!(p, child, params)
        end
    end

    return p.dv_G
end


"""
Directly calculates gravitational acceleration at particle p
used to test BH implementation above
"""
function a_G_direct(p::Particle, particles::Vector{Particle})
    a = zeros(3)

    for q in particles
        if p != q
            a += a_G(p, q)
        end
    end
    return a
end

""" 
Makes an octotree
"""
function make_tree(particles::Vector{Particle}, params)
    size = 5*params.R_virial

    tree = create_octree_node(zeros(3), size, 1) # Assuming particles are in the unit cube

    for particle in particles
        insert_particle!(tree, particle, 0)
    end 

    return tree
end


"""
Gravitational potential between two masses
"""
function U_G(p::Particle, q::Particle, params)
    eps = min(p.h, q.h) * params.r_plummer 
    r = q.x .- p.x
    u = -G*q.m*p.m /sqrt(norm(r)^2 + eps^2)
    return u
end


function U_G(p::Particle, node::OctreeNode, params)
    u = 0
    if node.children === nothing
        if node.particles != [] 
            for q in node.particles
                if q != p
                    u += U_G(p, q, params)
                end
            end
        end
        return u
    end

    distance = norm(node.center - p.x)
    if node.size / distance < params.theta
        # Node is sufficiently far away, use a single force calculation for the node
        # Don't use softening here
        u += -G * p.m* node.mass / distance
    else
        # Node is not sufficiently far away, recurse to children
        for child in node.children
            u += U_G(p, child, params)
        end
    end

    return u
end





function Base.show(io::IO, tree::OctreeNode)
    if tree.children == nothing 
        for q in  tree.particles
            print(io, "   |"^(tree.depth -1))
            print(io, "--")
            show(io, "text/plain", q)
            println(io)
        end
    else
        print(io, "   |"^(tree.depth-1))
        print(io, "-")
        @printf io " (%0.1f, %0.1f, %0.1f)\n"  (tree.center/pc)...

        for child in tree.children
            if child.particle == nothing && child.children != nothing
                print(io, "   |"^(tree.depth))
                println(io)
            end
            if child.children != nothing || child.particle != nothing
                print(io,  child)
            end
        end

    end

    return io
end



end 

