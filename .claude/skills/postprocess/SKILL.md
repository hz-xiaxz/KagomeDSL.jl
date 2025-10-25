---
name: postprocess
description: Reuse scripts and postprocess results (project)
---

# KagomeDSL.jl Postprocessing Protocol

## Overview

This skill provides **minimal, reusable templates** for running simulations and analyzing results. Templates are **directly editable** - Claude modifies parameters and runs them in place. No script duplication needed.

## Philosophy: Templates as Direct Tools

**Key principle:** Templates in `.claude/skills/postprocess/scripts/` are minimal, complete scripts. When you request a simulation or analysis:

1. **Claude copies template to `scratch/`** (gitignored folder)
2. **Edits parameters** in the scratch copy
3. **Runs from scratch/** - allows multiple simultaneous runs
4. **Results saved** to `/home/hzxiaxz/.julia/dev/KagomeDSL/data/` with full metadata
5. **Scratch files auto-cleaned** or kept for debugging

## Available Templates

### Simulation Templates

All save to `/home/hzxiaxz/.julia/dev/KagomeDSL/data/`:

| Template | Purpose | Key Parameters |
|----------|---------|----------------|
| `single_LL.jl` | Single Landau Level run | `n1`, `n2`, `imbalance` |
| `single_FP.jl` | Single Fermi Point run | `n1`, `n2`, `imbalance` |
| `LL_imbalance_sweep.jl` | LL with imbalance sweep | `n1`, `n2`, `imbalances` array |
| `LL_size_sweep.jl` | LL with size sweep | `sizes` array, `imbalance` |
| `FP_imbalance_sweep.jl` | FP with imbalance sweep | `n1`, `n2`, `imbalances` array |
| `zero_flux.jl` | Zero flux configuration | `n1`, `n2` |

**Key differences:**
- **LL (Landau Level)**: `B ≠ 0`, `N_up = N_down = ns÷2` (balanced particles)
- **FP (Fermi Point)**: `B = 0`, `N_up ≠ N_down` (imbalanced particles)

### Analysis Template

| Template | Purpose | Key Parameters |
|----------|---------|----------------|
| `analyze_results.jl` | Load and inspect results | `data_file` (e.g., "LL-8x8.results.json") |

## Claude Workflow

### When User Requests Simulation

**Example:** "Run LL simulation for 4x4 lattice with imbalances [0, 2, 4]"

**Claude actions:**
```bash
# 1. Copy template to scratch/:
cp .claude/skills/postprocess/scripts/LL_imbalance_sweep.jl scratch/LL_4x4_imb024.jl

# 2. Edit parameters in scratch copy:
#    - n1 = 4
#    - n2 = 4
#    - imbalances = [0, 2, 4]

# 3. Run from scratch/:
julia --project scratch/LL_4x4_imb024.jl run -r

# 4. Results automatically saved to data/LL-4x4.results.json
```

**Workflow:**
1. Identify appropriate template
2. Copy to `scratch/` with descriptive name using Bash
3. Edit only the "MODIFY THESE" section in scratch copy
4. Run from scratch/
5. Metadata is automatically saved in the results file

### When User Requests Analysis

**Example:** "Analyze LL-8x8 results and plot energy vs imbalance"

**Claude actions:**
```bash
# 1. Copy template to scratch/:
cp .claude/skills/postprocess/scripts/analyze_results.jl scratch/analyze_LL_8x8.jl

# 2. Edit in scratch copy:
#    - data_file = "LL-8x8.results.json"
#    - Add custom plotting code for energy vs imbalance

# 3. Run:
julia --project scratch/analyze_LL_8x8.jl

# 4. Figures saved to figures/
```

**Workflow:**
1. Copy analyze_results.jl to `scratch/` with descriptive name
2. Edit `data_file` parameter in scratch copy
3. Add custom analysis/plotting code if requested
4. Run and check output in `figures/`

## Standard Analysis Patterns

### Load Data

```julia
using Carlo, Carlo.ResultTools, DataFrames, Measurements

data_path = "/home/hzxiaxz/.julia/dev/KagomeDSL/data/LL-8x8.results.json"
df = DataFrame(Carlo.ResultTools.dataframe(data_path))

# Common columns: Splus_vector_real, Splus_vector_imag, energy, n1, n2, imbalance, B
```

### Extract Vector Observables

```julia
# Get S+ measurements with error bars
Sp_real = df[!, :Splus_vector_real][1]
Sp_imag = df[!, :Splus_vector_imag][1]

# Calculate phases (error propagation automatic)
phases = [atan(r, i)/π for (r, i) in zip(Sp_real, Sp_imag)]
```

### Parameter Sweep Analysis

```julia
# For parameter sweeps (multiple rows)
for row in eachrow(df)
    imbalance = row.imbalance
    energy = row.energy
    # Process each point
end
```

### Visualization with CairoMakie

```julia
using CairoMakie

fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], xlabel="Site", ylabel="Phase (π)")

x = 1:length(phases)
y = Measurements.value.(phases)
y_err = Measurements.uncertainty.(phases)

errorbars!(ax, x, y, y_err, whiskerwidth=10)
scatter!(ax, x, y, markersize=10)

save("/home/hzxiaxz/.julia/dev/KagomeDSL/figures/phases.png", fig)
```

** Don't read the picture **, just tell where to find.

## File Organization

```
KagomeDSL/
├── .claude/skills/postprocess/
│   └── scripts/              # Clean templates (never modified)
│       ├── single_LL.jl
│       ├── LL_imbalance_sweep.jl
│       ├── analyze_results.jl
│       └── ...
├── scratch/                  # Working copies (gitignored)
│   ├── LL_4x4_imb024.jl
│   ├── analyze_LL_8x8.jl
│   └── ...
├── data/                     # Simulation results (with metadata)
│   ├── LL-8x8.results.json
│   └── ...
├── figures/                  # Generated plots
│   └── ...
└── scripts/                  # Project-specific scripts (if needed)
```

**Why this works:**
- Templates stay clean - never modified
- Scratch copies are gitignored - edit freely
- Descriptive filenames in scratch/ for clarity
- Results include all metadata - no need to keep simulation scripts
- Absolute paths - run from anywhere
- Minimal token usage - only edit parameters, not whole scripts

## Common Pitfalls

1. **System size**: Must use `n1 = n2 = 4*n` (e.g., 4, 8, 12, 16)
2. **Imbalance**: Must be even numbers, ≤ `ns÷2` where `ns = n1 * n2 * 3`
3. **Error bars**: Always use Measurements.jl operations to preserve uncertainties
4. **File existence**: analyze_results.jl checks if data file exists before loading

## Quick Reference

**Run LL simulation:**
```bash
# Copy template to scratch/
cp .claude/skills/postprocess/scripts/single_LL.jl scratch/my_sim.jl
# Edit scratch/my_sim.jl (n1, n2, imbalance)
julia --project scratch/my_sim.jl run -r
```

**Analyze results:**
```bash
# Copy template to scratch/
cp .claude/skills/postprocess/scripts/analyze_results.jl scratch/my_analysis.jl
# Edit scratch/my_analysis.jl (data_file)
julia --project scratch/my_analysis.jl
```

**Check what's in a results file:**
```julia
using Carlo.ResultTools, DataFrames
df = DataFrame(ResultTools.dataframe("data/LL-8x8.results.json"))
names(df)  # Show all available measurements
```
