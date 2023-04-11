module Tree

using ..Particles
using ..Constants
using ..Init
using LinearAlgebra
using Printf

const max_depth = 30 # Adjust the maximum depth of the octree according to your problem

mutable struct OctreeNode
    children::Union{Nothing, Vector{OctreeNode}}
    particle::Union{Nothing, Particle}
    center::Vector{T}
    center_of_mass::Vector{T}

    mass::T
    size::T
    depth::Int
end


function create_octree_node(center, size, depth=0)
    return OctreeNode(nothing, nothing, center, zeros(3), 0, size, depth)
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
    if depth >= max_depth || node.children === nothing && node.particle === nothing
        node.particle = particle
        return
    end

    if node.children === nothing
        node.children = create_children(node)

        insert_particle!(node, node.particle, depth + 1)
        node.particle = nothing
    end

    index = 1
    index += 4 * (particle.x[1] < node.center[1])
    index += 2 * (particle.x[2] < node.center[2])
    index += 1 * (particle.x[3] < node.center[3])

    insert_particle!(node.children[index], particle, depth + 1)
end



function particles_within_range(p::Particle, node::OctreeNode, r::T)
    result = Particle[]

    function traverse(node::OctreeNode, r::T)
        if node.children === nothing
            if node.particle !== nothing && node.particle != target && norm(node.particle.x - target.x) <= r
                push!(result, node.particle)
            end
        end

        for child in node.children
            if norm(child.center - target.x) - child.size <= r
                traverse(child, r)
            end
        end
    end

    traverse(node, r)
    return result
end


"""
Gravitational acceleration between two masses
"""
function a_G(p::Particle, q::Particle)
    r_vec = q.x - p.x
    r_mag = norm(r_vec)
    eps = min(p.h, q.h)/2
    force_mag = G * q.m / sqrt(r_mag^2 + eps^2)
    return force_mag * r_vec / r_mag
end


function a_G(p::Particle, node::OctreeNode, theta::T)
    if node.children === nothing
        if node.particle === nothing || node.particle == p
            return zeros(3)
        else
            return a_G(p, node.particle)
        end
    end

    distance = norm(node.center - p.x)
    if node.size / distance < theta
        # Node is sufficiently far away, use a single force calculation for the node
        return G * node.mass * (node.center_of_mass - p.x) / distance^3
    else
        # Node is not sufficiently far away, recurse to children
        total_accel = zeros(3)
        for child in node.children
            total_accel += a_G(p, child, theta)
        end
        return total_accel
    end
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
function make_tree(particles::Vector{Particle})
    size = 20*Init.Rp

    tree = create_octree_node(zeros(3), size, 1) # Assuming particles are in the unit cube

    for particle in particles
        insert_particle!(tree, particle, 0)
    end 

    return tree
end






function Base.show(io::IO, tree::OctreeNode)
    if tree.children == nothing 
        if tree.particle != nothing
            print(io, "  "^(tree.depth))
            show(io, "text/plain", tree.particle)
        end
    else
        print(io, "  "^(tree.depth))
        print(io, "node $(tree.depth)  ")
        @printf io "(%0.1f, %0.1f, %0.1f)\n"  (tree.center/pc)...
        for child in tree.children
            print(io,  child)
        end
    end

    return io
end



end 

