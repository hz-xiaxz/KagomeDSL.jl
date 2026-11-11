# Claude AI Development Guide for KagomeDSL.jl

This document provides specific guidelines for Claude AI when developing this project. The goal is to ensure code quality, reproducibility, and adherence to Julia and computational physics best practices for simulating quantum spin systems on the Kagome lattice.

## Project Organization

The project follows the standard Julia package structure:

-   `src/`: Main source code for the project.
    -   `KagomeDSL.jl`: The main module file.
    -   `Lattice.jl`: Defines the Kagome lattice structure and connectivity.
    -   `Hamiltonian.jl`: Code related to defining the physical model (e.g., Heisenberg Hamiltonian).
    -   `MonteCarlo.jl`: Implementations of the Variational Monte Carlo (VMC) algorithm.
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

### Formatting Guidelines

- When writing Julia code, use `4 whitespaces` for indentation and try to keep
  the maximum line length under `92` characters.
- When writing Markdown text, use `2 whitespaces` for indentation and try to
  keep the maximum line length under `80` characters.
- When writing commit messages, follow the format `component: brief summary` for
  the title. In the body of the commit message, provide a brief prose summary of
  the purpose of the changes made.
  Also, ensure that the maximum line length never exceeds 72 characters.
  When referencing external GitHub PRs or issues, use proper GitHub interlinking
  format (e.g., `owner/repo#123` for PRs/issues).

### Julia Best Practices

-   **Type Stability**: This is critical for performance. Write functions where the output type can be inferred from the input types. Use `@code_warntype` to check for type instabilities in performance-critical code. Refer to the [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/).
-   **Dispatch, Not Conditionals**: Use Julia's multiple dispatch to handle different models or algorithms. Instead of `if model_type == :heisenberg`, define separate methods: `run_mc(model::Heisenberg, ...)` and `run_mc(model::Ising, ...)`.
-   **In-place Operations**: For large arrays (like spin configurations), use in-place functions (e.g., `mul!`, `copyto!`, functions ending in `!`) to reduce memory allocations and improve speed.
-   **Modularity**: Keep different concepts in separate files. The `Hamiltonian` definition is separate from the `MonteCarlo` solver, allowing the same solver to be used for different models.

### Coding Rules

- When writing functions, use the most restrictive signature type possible.
  This allows JET to easily catch unintended errors.
  Of course, when prototyping, it's perfectly fine to start with loose type
  declarations, but for the functions you ultimately commit, it's desirable to
  use type declarations as much as possible.
  Especially when AI agents suggest code, please make sure to clearly
  specify the argument types that functions expect.
  In situations where there's no particular need to make a function generic, or
  if you're unsure what to do, submit the function with the most restrictive
  signature type you can think of.

- For function calls with keyword arguments, use an explicit `;` for clarity.
  For example, code like this:
  ```julia
  ...
  Position(; line=i-1, character=m.match.offset-1)
  ...
  ```
  is preferred over:
  ```julia
  ...
  Position(line=i-1, character=m.match.offset-1)
  ...
  ```

- **ONLY INCLUDE COMMENTS WHERE TRULY NECESSARY**.
  When the function name or implementation clearly indicates its purpose or
  behavior, redundant comments are unnecessary.

- On the other hand, for general utilities that expected to be used in multiple
  places, it's fine to use docstrings to clarify their behavior. However, even
  in these cases, if the function name and behavior are self-explanatory, no
  special docstring is needed.

### Physics & Algorithm-Specific Guidance

-   **Model Definition**: Physical models are defined by the `Hamiltonian` and the `Lattice`. These should contain all necessary parameters (e.g., `J`, lattice size), making them easy to pass to solver functions.
-   **Algorithm Abstraction**: The core logic of an algorithm (e.g., a VMC update step) is separated from the specific model it's applied to. This is a key design principle of the project.
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

#### Running Specific Tests

If explicit test file or code is provided, prioritize running that.
Otherwise, you can run the entire test suite for the project by executing
`using Pkg; Pkg.test()` from the root directory of this repository.

For example, if you receive a prompt like this:
> Improve the error message of diagnostics.
> Use test/test_diagnostics for the test cases.

The command you should run is:
```bash
julia --startup-file=no -e 'using Test; @testset "test_diagnostics" include("test/test_diagnostics")'
```
Note that the usage of the `--startup-file=no` flag, which avoids loading
unnecessary startup utilities.

### What to Test

-   **Analytical Limits**: Test your code against known analytical solutions. For example, the non-interacting limit or limits for small system sizes (e.g., a single triangle).
-   **Symmetries**: If your model has a conservation law (e.g., total `Sz`), write a test to ensure your algorithm respects it.
-   **Benchmarks**: Compare results against established results from papers or other well-tested codes where possible.
-   **Type Stability**: You can add `@inferred` from the `Test` standard library to your tests to ensure key functions are type-stable.

### Test Code Organization

Testing functionality is challenging. To fully test such functionality, you need
to set up proper test environments and send requests that mimic realistic usage.

However, writing comprehensive tests can be tricky. Therefore, unless explicitly
requested by the core developers, you don't need to write test code to fully
test newly implemented features. It's generally sufficient to test important
subroutines that are easy to test in the implementation.

Test code should be written in files that define independent module spaces with
a `test_` prefix. Then include these files from [`test/runtests.jl`](./test/runtests.jl).
This ensures that these files can be run independently from the REPL.

In each test file, you are encouraged to use `@testset "testset name"` to
organize our tests cleanly. For code clarity, unless specifically necessary,
avoid using `using`, `import`, and `struct` definitions  inside `@testset`
blocks, and instead place them at the top level.

Also, you are encouraged to use `let`-blocks to ensure that names aren't
unintentionally reused between multiple test cases.

For example, here is what good test code looks like:
```julia
module test_physics

using Test

function testcase_util(param)
    # Test utility function - calculate something physics-related
    return param^2  # Example: energy calculation
end

@testset "physics_calculations" begin
    let param = 1.0
        result = testcase_util(param)
        @test result ≈ 1.0
    end
    let param = 2.0
        result = testcase_util(param)
        @test result ≈ 4.0
    end
end

end # module test_physics
```

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

## Environment-Related Issues

For AI agents: **NEVER MODIFY [Project.toml](./Project.toml) OR [test/Project.toml](./test/Project.toml) BY YOURSELF**.
If you encounter errors that seem to be environment-related when running tests,
in most cases this is due to working directory issues, so first `cd` to the root directory of this project
and re-run the tests. Never attempt to fix environment-related issues yourself.
If you cannot resolve the problem, inform the human engineer and ask for instructions.

## About Modifications to Code You've Written

If you, as an AI agent, add or modify code, and the user appears to have made
further manual changes to that code after your response, please respect those
modifications as much as possible.
For example, if the user has deleted a function you wrote, do not reintroduce
that function in subsequent code generation.
If you believe that changes made by the user are potentially problematic,
please clearly explain your concerns and ask the user for clarification.

## Commit Message Formatting

Follow these guidelines to maintain a clean and informative git history.

-   **Title**: Use the format `scope: Brief summary`. The scope should be a component of the project.
    -   **Examples**:
        -   `feat(mc): Add parallel tempering support to VMC`
        -   `fix(hamiltonian): Correct sign in transverse field term`
        -   `docs(readme): Update installation instructions`
        -   `test(lattice): Add test for periodic boundary conditions`
        -   `refactor(measure): Improve performance of correlation function`
-   **Body**: In the commit body, explain the "what" and "why" of the change. Do not just repeat the "how". If the change fixes a GitHub issue, add `Fixes #123` to the end of the body.
# Tips

## Performance Tips
1. Use `eachindex(iterator)` instead of `1:length(iterator)` to avoid unnecessary allocations
2. Use `@views` to avoid unnecessary allocations when slicing arrays
3. Use `axes(array, dim)` instead of `1:size(array, dim)` for better performance and idiomatic Julia code.

## Visualization Tips  
3. Always use `Makie` ecosystem for visualization:
   - 2D plots: Use `CairoMakie`
   - 3D/interactive plots: Use `GLMakie` or `WGLMakie`
4. Use `Figure(size=...)` instead of deprecated `Figure(resolution=...)`
5. Use `arrows2d`/`arrows3d` instead of deprecated `arrows` function

## Code Style Tips
6. When writing scripts, use `!` instead of invalid `\!` notation

## Documentation Tips
7. Use double backticks `` `` to document LaTeX-formatted mathematical equations
   - Example: ``⟨x|cᵢ⁺cⱼ|ψ⟩/⟨x|ψ⟩`` instead of ⟨x|cᵢ⁺cⱼ|ψ⟩/⟨x|ψ⟩