#! /usr/bin/env -S julia --startup-file=no --color=yes
using KagomeDSL
using Carlo
using Carlo.JobTools
using Dates
using LinearAlgebra: eigvals

tm = TaskMaker()
tm.thermalization = 2000
tm.sweeps = 100000 
tm.binsize = 100
tm.n1 = 4 
tm.n2 = 2
ns = tm.n1 * tm.n2 * 3
tm.PBC = (false, false)

# pre check shell
lat = DoubleKagome(1.0, tm.n1, tm.n2, tm.PBC)
H = KagomeDSL.Hmat(lat)
E = sort(eigvals(H))
num = findlast(x -> isapprox(x, E[1], atol = 1e-10), E)
# the number of N_up and N_down should be at least > num
# if num < ns ÷ 2, i is from ns-num to num

for i = (ns-num):-1:num
    task(tm; N_up = i, N_down = ns - i)
end


dir = @__DIR__
# savepath = dir * "/../data/" * Dates.format(Dates.now(), "mm-ddTHH-MM-SS")
savepath = dir * "/../data/" * "mpi$(tm.n1)x$(tm.n2)"
job = JobInfo(
    savepath,
    KagomeDSL.MC;
    tasks = make_tasks(tm),
    checkpoint_time = "5:00",
    run_time = "24:00:00",
)

Carlo.start(job, ARGS)
