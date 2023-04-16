module Tree

using ..Particles
using ..Constants
using ..Init
using LinearAlgebra
using Printf

export make_tree, a_G, find_within_r, a_G!

const max_depth = 30 # Adjust the maximum depth of the octree according to your problem

mutable struct OctreeNode
    children::Union{Nothing, Vector{OctreeNode}}
    particle::Union{Nothing, Particle}
    center::Vector{F}
    center_of_mass::Vector{F}

    mass::F
    size::F
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



function find_within_r(p::Particle, node::OctreeNode, r::F)
    result = Particle[]

    function traverse(node::OctreeNode, r::F)
        if node.children === nothing
            if node.particle !== nothing && node.particle != p && norm(node.particle.x - p.x) <= r
                push!(result, node.particle)
            end
            return
        end

        for child in node.children
            if norm(child.center - p.x) - child.size <= r
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

function a_G!(a, p::Particle, q::Particle, params)
    eps = min(p.h, q.h) * params.r_plummer 
    r = q.x .- p.x
    a .+= G*q.m * normalize(r)/(sum(r.*r) + eps^2)
    return a
end


function a_G!(a, p::Particle, node::OctreeNode, params)
    if node.children === nothing
        if node.particle === nothing || node.particle == p
        else
            a_G!(a, p, node.particle, params)
        end
        return a
    end

    distance = norm(node.center - p.x)
    if node.size / distance < params.theta
        # Node is sufficiently far away, use a single force calculation for the node
        # Don't use softening here
        a .+= G * node.mass * (node.center_of_mass - p.x) / distance^3
    else
        # Node is not sufficiently far away, recurse to children
        for child in node.children
            a_G!(a, p, child, params)
        end
    end

    return a
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






function Base.show(io::IO, tree::OctreeNode)
    if tree.children == nothing 
        if tree.particle != nothing
            print(io, "   |"^(tree.depth -1))
            print(io, "--")
            show(io, "text/plain", tree.particle)
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

