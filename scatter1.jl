using MPI
MPI.Init()

import Base.+

struct Point{T}
    x::T
    y::T
end

+(A::Point{T}, B::Point{T}) where T = Point{T}(A.x + B.x, A.y + B.y)


comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

p = Point(rank, rank) 
print("Original Point on rank $(rank):\n  $(p)\n")
    

recv_buf = MPI.Reduce(p, +, comm; root=0)

if rank == 0
    print("\nNew Point on rank 0:\n  $(recv_buf)\n")
end
