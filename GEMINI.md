# Agent Development Guide for KagomeDSL.jl

This document provides guidelines for an AI agent developing this project. The goal is to ensure code quality, reproducibility, and adherence to Julia and computational physics best practices for simulating quantum spin systems on the Kagome lattice.


## Project Organization

The project follows the standard Julia package structure:

-   `src/`: Main source code for the project.
    -   `KagomeDSL.jl`: The main module file.
    -   `Lattice.jl`: Defines the Kagome lattice structure and connectivity.
    -   `Hamiltonian.jl`: Code related to defining the physical model (e.g., Heisenberg Hamiltonian).
    -   `MonteCarlo.jl`: Implementations of the Stochastic Series Expansion (SSE) Monte Carlo algorithm.
    -   `Measurements.jl`: Functions for calculating physical observables (e.g., spin correlations, energy).
-   `test/`: Project tests.
    -   `runtests.jl`: The main test runner.
    -   Tests are organized per file to mirror the `src/` directory structure (e.g., `test-Hamiltonian.jl`).
-   `docs/`: Documentation sources for the `Documenter.jl`-based website.
-   `scripts/`: Standalone scripts for running simulations and post-processing results. This is the primary way to generate the data found in `data/`.
-   `data/`: For storing simulation results. Raw data is not committed to Git.
-   `Project.toml` & `Manifest.toml`: Defines project dependencies and their exact versions. (!NEVER EDIT `Project.toml`!)

## Development Workflow

Follow these steps for any change, from a small bugfix to a new feature.

### 1. Activate the Project Environment

Before starting, ensure the project's dependencies are loaded correctly. In a Julia REPL from the project root:

```julia
julia --project
```

### 2. Code, Test

The core development loop is to:
1.  Write/modify code in `src/`.
2.  Write or update corresponding tests in `test/`.

## Writing Code: Best Practices

### Julia Best Practices

-   **Type Stability**: This is critical for performance. Write functions where the output type can be inferred from the input types. Use `@code_warntype` to check for type instabilities in performance-critical code. Refer to the [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/).
-   **Dispatch, Not Conditionals**: Use Julia's multiple dispatch to handle different models or algorithms. Instead of `if model_type == :heisenberg`, define separate methods: `run_mc(model::Heisenberg, ...)` and `run_mc(model::Ising, ...)`.
-   **In-place Operations**: For large arrays (like spin configurations), use in-place functions (e.g., `mul!`, `copyto!`, functions ending in `!`) to reduce memory allocations and improve speed.
-   **Modularity**: Keep different concepts in separate files. The `Hamiltonian` definition is separate from the `MonteCarlo` solver, allowing the same solver to be used for different models.

### Physics & Algorithm-Specific Guidance

-   **Model Definition**: Physical models are defined by the `Hamiltonian` and the `Lattice`. These should contain all necessary parameters (e.g., `J`, lattice size), making them easy to pass to solver functions.
-   **Algorithm Abstraction**: The core logic of an algorithm (e.g., an SSE update step) is separated from the specific model it's applied to. This is a key design principle of the project.
-   **Reproducibility**:
    -   All simulations must be runnable from a script in `scripts/` that saves the parameters used.
    -   **Always seed the random number generator** (using `Random; Random.seed!(1234)`) for Monte Carlo simulations to ensure results are reproducible.
    -   Save results with the parameters that generated them. Use formats like JLD2.jl or BSON.jl to save both the data and the parameter struct in the same file.

## Testing

Comprehensive testing is non-negotiable for ensuring the physical correctness of the simulations.

### Running Tests

To run the full test suite, execute 

```
julia --project -e 'using Pkg; Pkg.test()'
```

### What to Test

-   **Analytical Limits**: Test your code against known analytical solutions. For example, the non-interacting limit or limits for small system sizes (e.g., a single triangle).
-   **Symmetries**: If your model has a conservation law (e.g., total `Sz`), write a test to ensure your algorithm respects it.
-   **Benchmarks**: Compare results against established results from papers or other well-tested codes where possible.
-   **Type Stability**: You can add `@inferred` from the `Test` standard library to your tests to ensure key functions are type-stable.

## Commit Message Formatting

Follow these guidelines to maintain a clean and informative git history.

-   **Title**: Use the format `scope: Brief summary`. The scope should be a component of the project.
    -   **Examples**:
        -   `feat(mc): Add parallel tempering support to SSE`
        -   `fix(hamiltonian): Correct sign in transverse field term`
        -   `docs(readme): Update installation instructions`
        -   `test(lattice): Add test for periodic boundary conditions`
        -   `refactor(measure): Improve performance of correlation function`
-   **Body**: In the commit body, explain the "what" and "why" of the change. Do not just repeat the "how". If the change fixes a GitHub issue, add `Fixes #123` to the end of the body.
