#! /usr/bin/env -S julia --startup-file=no --color=yes
using KagomeDSL
using Carlo
using Carlo.JobTools
using Dates
using LinearAlgebra
using ArnoldiMethod

tm = TaskMaker()
tm.thermalization = 5000
tm.sweeps = 10_000
tm.binsize = 10
tm.n1 = 12
tm.n2 = 12
tm.link_in = Dict(
    (1, 2) => 1,
    (1, 3) => 1,
    (2, 3) => 1,
    (2, 4) => 1,
    (4, 6) => 1,
    (4, 5) => 1,
    (5, 6) => 1,
    (2, 1) => 1,
    (3, 1) => 1,
    (3, 2) => 1,
    (4, 2) => 1,
    (6, 4) => 1,
    (5, 4) => 1,
    (6, 5) => 1,
)
tm.link_inter = Dict(
    (3, 5, -1, 1) => 1,
    (3, 1, 0, 1) => 1,
    (6, 2, 0, 1) => 1,
    (6, 4, 0, 1) => 1,
    (5, 1, 1, 0) => 1,
    (1, 5, -1, 0) => 1,
    (1, 3, 0, -1) => 1,
    (2, 6, 0, -1) => 1,
    (4, 6, 0, -1) => 1,
    (5, 3, 1, -1) => 1,
)
ns = tm.n1 * tm.n2 * 3
tm.PBC = (true, true)
tm.antiPBC = (true, false)
task(tm; N_up = ns ÷ 2, N_down = ns ÷ 2)

dir = @__DIR__
# savepath = dir * "/../data/" * Dates.format(Dates.now(), "mm-ddTHH-MM-SS")
savepath = dir * "/../data/" * "zero-$(tm.n1)x$(tm.n2)"
job = JobInfo(
    savepath,
    KagomeDSL.MC;
    tasks = make_tasks(tm),
    checkpoint_time = "30:00",
    run_time = "24:00:00",
)

Carlo.start(job, ARGS)
