#!/usr/bin/env julia

# Copyright 2022 John T. Foster
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
module project2

using MPI
using RecipesBase
using Plots

struct Grid{T}
    dx::T
    dy::T
    xrange::NTuple{2, T}
    yrange::NTuple{2, T}
    u::Matrix{T}
end

function Grid(nx::Integer=10, ny::Integer=10, 
              xrange::NTuple{2, Real}=(0.0, 1.0),
              yrange::NTuple{2, Real}=(0.0, 1.0))
    xmin, xmax = xrange
    ymin, ymax = yrange
    dx = (xmax - xmin) / (nx - 1)
    dy = (ymax - ymin) / (ny - 1)
    t = typeof(dx)
    u = zeros(t, (nx, ny))
    Grid{t}(dx, dy, xrange, yrange, u)
end


function set_boundary_condition!(g::Grid, f::Function=(x,y) -> 0.0; 
        sides::Union{String, Tuple{Vararg{String}}})

    xmin, xmax = g.xrange
    ymin, ymax = g.yrange
    nx = size(g.u)[2]
    ny = size(g.u)[1]

    if sides == "all"
        sides = ("bottom", "top", "left", "right")
    elseif typeof(sides) == String
        sides = (sides,)
    end

    for side in sides

        if side == "left"
            g.u[:, 1] = f.(xmin, LinRange(ymin, ymax, ny))
        elseif side == "right"
            g.u[:, end] = f.(xmax, LinRange(ymin, ymax, ny))
        elseif side == "bottom"
            g.u[1, :] = f.(LinRange(xmin, xmax, nx), ymin)
        elseif side == "top"
            g.u[end, :] = f.(LinRange(xmin, xmax, nx), ymax)
        end
    end

end

function reset!(g::Grid)
    g.u .= zero(typeof(g.dx))
end

function iterate!(g::Grid)

    u = g.u
    nx, ny = size(u)

    dx2, dy2 = g.dx^2, g.dy^2

    err = 0.0

    for i = 2:(nx-1)
        for j = 2:(ny-1)

            tmp = u[i, j]

            u[i,j] = ((u[i-1, j] + u[i+1, j]) * dy2 +
                      (u[i, j-1] + u[i, j+1]) * dx2) / (dx2 + dy2) / 2
                
            diff = u[i,j] - tmp
                
            err += diff * diff
        end
    end

    sqrt(err)
end

function solve!(g::Grid; max_iterations::Integer=10000, 
                         tolerance::AbstractFloat=1.0e-12,
                         quiet::Bool=false)

    for i=1:max_iterations
        error = iterate!(g)
        if error < tolerance
            if !quiet
                println("Solution converged in $(i) iterations")
            end
            break
        end
    end
end


function partition_grid(g::Union{Grid, Nothing}, comm::MPI.Comm)::Grid
   # Return Grid{t}(dx, dy, xrange, yrange, my_u) on each rank
   rank = MPI.Comm_rank(comm)
   size_world = MPI.Comm_size(comm)
   ny = nothing
   dx = nothing
   dy = nothing
   xrange = nothing
   yrange = nothing
   lengths = nothing

    if rank == 0
        nx = size(g.u)[2]
        ny = size(g.u)[1]
        dx = g.dx
        dy = g.dy

        xrange = g.xrange
        yrange = g.yrange

        remainder = nx % size_world
        chunk_size = nx ÷ size_world
        
        lengths = ones(Int, size_world) .* chunk_size .* ny #number of entries
        lengths[1:remainder] .+= ny
        # lengths[2:end] .+= ny
    end

    lengths = MPI.bcast(lengths, 0, comm)
    ny = MPI.bcast(ny, 0, comm)
    dx = MPI.bcast(dx, 0, comm)
    dy = MPI.bcast(dy, 0, comm)
    xrange = MPI.bcast(xrange, 0, comm)
    yrange = MPI.bcast(yrange, 0, comm)

    t = typeof(dx)
    if rank == 0 || rank == size_world -1
        recv_buf = zeros(t, (lengths[rank+1] + ny))
    else
        recv_buf = zeros(t, (lengths[rank+1] + 2 * ny))
    end

    if rank == 0
        send_buf = MPI.VBuffer(g.u, lengths)
        MPI.Scatterv!(send_buf, recv_buf, 0, comm)
    else
        send_buf = nothing
        MPI.Scatterv!(send_buf, view(recv_buf, (ny+1):lastindex(recv_buf)), 0, comm)

    end

    my_u = reshape(recv_buf, ny, :)
    Grid{t}(dx, dy, xrange, yrange,my_u)

    
end

function solve!(my_g::Grid, comm::MPI.Comm; max_iterations::Integer=10000, 
                                            tolerance::AbstractFloat=1.0e-12,
                                            quiet::Bool=false)
   #Returns nothing, should update my_g in place
   rank = MPI.Comm_rank(comm)
   size_world = MPI.Comm_size(comm)
   for i=1:max_iterations
        if rank > 0
            MPI.Send(view(my_g.u, :, 2), comm; dest = rank-1)
            MPI.Recv!(view(my_g.u, :, 1), comm; source=rank-1)
        end
        if rank < size_world -1
            MPI.Send(view(my_g.u, :, size(my_g.u)[2]-1), comm; dest=rank+1)
            MPI.Recv!(view(my_g.u, :, size(my_g.u)[2]), comm; source=rank+1)
        end
        
        my_error = iterate!(my_g) ^ 2
        if i % 50 ==0
            error = MPI.Allreduce(my_error, +, comm)
            
            total_error = sqrt(sum(error))
            
            if total_error < tolerance
                if !quiet && rank == 0
                    println("Solution converged in $(i) iterations")
                end
                break
            end
        end
    end
end

function get_solution!(g::Union{Grid, Nothing}, my_g::Grid, comm::MPI.Comm)
   #Returns nothing, should update g in place
   rank = MPI.Comm_rank(comm)
   size_world = MPI.Comm_size(comm)

    if rank == 0
        send_buf = view(my_g.u, :, 1:size(my_g.u)[2]-1)
    elseif rank == size_world - 1
        send_buf = view(my_g.u, :, 2:size(my_g.u)[2])
    else
        send_buf = view(my_g.u, :, 2:size(my_g.u)[2]-1)
    end
    
    len = MPI.Allgather(length(send_buf), comm)
    

    if rank == 0
        recv_buf = MPI.VBuffer(g.u,len)
    else
        recv_buf = nothing
    end

    MPI.Gatherv!(send_buf, recv_buf, 0, comm)   
   
end

@recipe function f(g::Grid)
    xmin, xmax = g.xrange
    ymin, ymax = g.yrange
    nx = size(g.u)[2]
    ny = size(g.u)[1]
    LinRange(xmin, xmax, nx), LinRange(ymin, ymax, ny), g.u
end

function run_serial()
    g = Grid(101, 101)
    set_boundary_condition!(g, (x, y) -> 10, sides=("top", "bottom"))
    solve!(g)

    ENV["GKSwstype"]="nul"
    p = contourf(g, aspect_ratio=:equal, xlims=g.xrange);
    savefig(p, "laplace_serial.png");
end


function run_parallel()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    ENV["GKSwstype"]="nul"

    if rank == 0
        g = Grid(101, 101)
        set_boundary_condition!(g, (x, y) -> 10, sides=("top", "bottom"))
    else
        g = nothing
    end
    my_g = partition_grid(g, comm)
    solve!(my_g, comm)
    # p = contourf(my_g, xlims=my_g.xrange, aspect_ratio=:equal)
    # savefig(p, "laplace_parallel_$(rank).png");
    get_solution!(g, my_g, comm)
    if rank == 0
        p = contourf(g, aspect_ratio=:equal, xlims=g.xrange);
        savefig(p, "laplace_parallel.png");
    end
end

export run_serial, run_parallel, set_boundary_condition!, partition_grid, Grid, solve!, get_solution!

end
