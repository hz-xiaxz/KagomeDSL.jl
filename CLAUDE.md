# Claude AI Development Guide for KagomeDSL.jl

This document provides specific guidelines for Claude AI when developing this project. The goal is to ensure code quality, reproducibility, and adherence to Julia and computational physics best practices for simulating quantum spin systems on the Kagome lattice.

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

If any Package is missing, install it using:

```julia
julia --project -e 'import Pkg; Pkg.add("PackageName")'
```

Never edit `Project.toml` directly; use the Julia REPL to manage dependencies.

If running file in subdirectories, create new ``Project.toml`` files in those directories to avoid conflicts with the main project environment. 
```julia
julia --project=./subdirectory -e 'import Pkg; Pkg.activate("."); Pkg.instantiate()'
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

To run the full test suite, execute:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

### What to Test

-   **Analytical Limits**: Test your code against known analytical solutions. For example, the non-interacting limit or limits for small system sizes (e.g., a single triangle).
-   **Symmetries**: If your model has a conservation law (e.g., total `Sz`), write a test to ensure your algorithm respects it.
-   **Benchmarks**: Compare results against established results from papers or other well-tested codes where possible.
-   **Type Stability**: You can add `@inferred` from the `Test` standard library to your tests to ensure key functions are type-stable.

## Claude-Specific Guidelines

### File Operations
- Always use the Read tool to examine existing files before making changes
- Use Edit or MultiEdit tools to modify existing files instead of creating new ones
- Follow the project's existing code style and conventions

### Task Management
- Use the TodoWrite tool to plan and track complex tasks
- Mark tasks as in_progress when starting work and completed when finished
- Break down large features into smaller, manageable steps

### Code Quality
- Run tests after making changes using: `julia --project -e 'using Pkg; Pkg.test()'`
- Check for type stability in performance-critical code
- Follow Julia naming conventions and existing project patterns

### Documentation
- Only create documentation files if explicitly requested
- Focus on code implementation rather than extensive documentation
- Use docstrings for public API functions

### Version Control
- Only commit changes when explicitly asked by the user
- Use meaningful commit messages following the project's format
- Check git status before making commits

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