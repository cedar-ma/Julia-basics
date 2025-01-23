using MPI

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

v = collect(1:size)
print("Original array on rank $(rank):\n  $(v)\n")
    

recv_buf = MPI.Reduce(v, +, comm; root=0)

if rank == 0
    print("New array on rank 0:\n  $(recv_buf)\n")
    total_sum = sum(recv_buf)
    print("Total sum:\n  $(total_sum)")
end
