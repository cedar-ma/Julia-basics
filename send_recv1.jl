using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

send_buf = Array{Float64}(undef,10)
recv_buf = Array{Float64}(undef,10)

fill!(send_buf, Float64(rank))

if rank == 0
    send_status = MPI.Isend(send_buf,comm; dest=1)
elseif rank>0
    recv_status = MPI.Irecv!(recv_buf, comm; source=(rank-1))
    send_status = MPI.Isend(send_buf,comm; dest = (rank+1)%size)
end
if rank ==0
    recv_status = MPI.Irecv!(recv_buf, comm; source = (size-1))
end

status = MPI.Waitall!([send_status, recv_status])

print("MY rank is $(rank)\n I received this: $(recv_buf)\n")
