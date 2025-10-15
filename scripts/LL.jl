#! /usr/bin/env -S julia --project --startup-file=no --color=yes
# to run, use:
# ./job.jl run -r
using KagomeDSL
using Carlo
using Carlo.JobTools
using Dates
using LinearAlgebra

tm = TaskMaker()
tm.thermalization = 1000
tm.sweeps = 5_000
tm.binsize = 2
tm.n1 = 8 
tm.n2 = 8 
ns = tm.n1 * tm.n2 * 3
tm.PBC = (true, true)
tm.antiPBC = (false, true)
tm.lattice = DoubleKagome
imbalances = [0, 2, 4, 8, 10, 14]
for imbalance in imbalances
    S = tm.n1 * tm.n2 * 2√3
    tm.B = imbalance * π / S
    tm.imbalance = imbalance
    task(tm; N_up = ns ÷ 2 , N_down = ns ÷ 2 )
end

dir = @__DIR__
# savepath = dir * "/../data/" * Dates.format(Dates.now(), "mm-ddTHH-MM-SS")
savepath = dir * "/../data/" * "LL-$(tm.n1)x$(tm.n2)"
job = JobInfo(
    savepath,
    KagomeDSL.MC;
    tasks = make_tasks(tm),
    checkpoint_time = "30:00",
    run_time = "24:00:00",
)

Carlo.start(job, ARGS)
