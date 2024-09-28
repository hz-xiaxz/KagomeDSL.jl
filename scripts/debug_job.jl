using KagomeDSL
using Carlo
using Carlo.JobTools
using Dates
using LinearAlgebra: eigvals
using Logging

tm = TaskMaker()
tm.thermalization = 100
tm.sweeps = 10
tm.binsize = 1
tm.n1 = 8
tm.n2 = 8
ns = tm.n1 * tm.n2 * 3
tm.PBC = (true, false)

# pre check shell
lat = DoubleKagome(1.0, tm.n1, tm.n2, tm.PBC)
H = KagomeDSL.Hmat(lat)
E = sort(eigvals(H))
shell_pool = []
# iteratively find degenerate spaces
start_shell = 1
while start_shell < length(E)
    global start_shell
    num = findlast(x -> isapprox(x, E[start_shell], atol = 1e-10), E)
    push!(shell_pool, (start_shell, num))
    start_shell = num + 1
end
# the number of N_up and N_down should be at least > num
# if num < ns ÷ 2, i is from ns-num to num
first_num = shell_pool[1][2]
for i = first_num:(ns÷2)
    # find i in the shell scope
    shell = filter(
        x -> (x[1] <= i && x[2] > i) || (x[1] <= (ns - i) && x[2] > (ns - i)),
        shell_pool,
    )
    if isempty(shell)
        task(tm; N_up = i, N_down = ns - i)
    end
end


dir = @__DIR__
# savepath = dir * "/../data/" * Dates.format(Dates.now(), "mm-ddTHH-MM-SS")
savepath = dir * "/../data/" * "debug$(tm.n1)x$(tm.n2)PBC"
job = JobInfo(
    savepath,
    KagomeDSL.MC;
    tasks = make_tasks(tm),
    checkpoint_time = "30:00",
    run_time = "24:00:00",
)

Carlo.cli_delete(job, Dict())
with_logger(Carlo.default_logger()) do
    start(Carlo.SingleScheduler, job)
end
