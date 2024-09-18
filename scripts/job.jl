#!/usr/bin/env -S julia --color=yes --startup-file=no

using KagomeDSL
using Carlo
using Carlo.JobTools
using Dates

tm = TaskMaker()
tm.thermalization = 0
tm.sweeps = 100
tm.binsize = 1
tm.n1 = 8
tm.n2 = 8
ns = tm.n1 * tm.n2 * 3
tm.PBC = (true, false)
tm.χ = 1.0
for i = ns÷2:-1:2
    task(tm; N_up = i, N_down = ns - i)
end

dir = @__DIR__
# savepath = dir * "/../data/" * Dates.format(Dates.now(), "mm-ddTHH-MM-SS")
savepath = dir * "/../data/" * "$(tm.n1)x$(tm.n1)"
job = JobInfo(
    savepath,
    KagomeDSL.MC;
    tasks = make_tasks(tm),
    checkpoint_time = "30:00",
    run_time = "24:00:00",
)

Carlo.start(job, ARGS)
