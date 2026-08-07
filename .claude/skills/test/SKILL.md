name: Kagome-test
description: Test based developer protocol for KagomeDSL.jl, especially when adding new features or debugging.
---

# KagomeDSL.jl Developer Protocol

## Overview

Your mission is to develop `KagomeDSL.jl`, a high-performance Julia package for simulating quantum spin systems on the Kagome lattice. You will act as an expert in both Julia and computational physics.

**Core principles**: Every line of code you write must be:
1.  **Physically Correct**: Adheres to the model and simulation constraints.
2.  **Reproducible**: Can be run by others to get the exact same result.
3.  **Maintainable**: Clear, well-tested, and idiomatic.

Violating the letter of the mandatory rules is violating the spirit of the rules.

## The Iron Laws

1.  **NO MANUAL EDITS TO `Project.toml` OR `Manifest.toml`**. Use the Julia Pkg REPL for all dependency management. Period.
2.  **NO UN-REPRODUCIBLE SIMULATIONS**. All Monte Carlo scripts in `scripts/` MUST set a random seed.
3.  **NO PRODUCTION CODE IN `src/` WITHOUT A CORRESPONDING TEST IN `test/`**.

Wrote code that violates these laws? Delete it. Start over.

## The Development Cycle

The Three-Step Cycle (Red-Green-Refactor)
Follow these three steps in order for every small feature you add.
Step 1: RED - Write a Failing Test
Your Action: Write a new, small test that describes the feature you are about to build.
The Goal: This test must fail because you haven't written the code for the feature yet. This is the "Red" phase.
Verification: Run all tests and confirm that only the new test fails, and it fails for the expected reason (e.g., a function is not defined, or the result is incorrect).
Step 2: GREEN - Make the Test Pass
Your Action: Write the simplest possible code to make the failing test pass.
The Goal: Do not add any extra logic or features. Just do the minimum work required to turn the test from "Red" to "Green".
Verification: Run all tests again and confirm that they all pass now.
Step 3: REFACTOR - Improve the Code
Your Action: Now that you have a passing test as a safety net, you can clean up and improve the code you just wrote. You can also improve any other related code.
The Goal: Make the code more readable, efficient, and remove any duplication. You do this without changing what the code does.
Verification: Run the tests frequently as you refactor. They must stay "Green" at all times. If a test fails, you know you broke something and should fix it immediately.

## Running Tests and Interpreting Output

To run the test suite, execute the following command from the project root:
```bash
julia --project -e "using Pkg; Pkg.test()"
```

### Understanding Test Output

When you run the tests, you'll see output indicating the status of each test set.

*   **Passing Tests**: A passing test set will be marked with `Test Summary: | Pass`.
*   **Failing Tests**: A failing test will be clearly marked with `Test Summary: | Fail`. The output will show the expression that failed, the expected outcome, and the actual result. Look for the `@test` macro that failed.

Example of a failing test output:
```
Test Failed at /path/to/KagomeDSL/test/test-Lattice.jl:25
  Expression: get_neighbor(lattice, 1, 1) == 2
   Evaluated: 3 == 2
```
This tells you:
1.  The failure occurred in `test/test-Lattice.jl` on line 25.
2.  The test `get_neighbor(lattice, 1, 1) == 2` failed.
3.  The actual result was `3`, but the test expected `2`.

### Handling Failing Tests

When a test fails, your goal is to make it pass. However, sometimes you might encounter failures in tests that are unrelated to your current changes.

1.  **Focus on the relevant failure**: Identify the test that is failing due to your recent changes.
2.  **Isolate the test**: If other tests are failing and creating noise, you can temporarily comment them out. This helps you focus on the specific problem you are trying to solve.
    *   In the test file (e.g., `test/test-Lattice.jl`), find the `@testset` or `@test` that is failing.
    *   You can comment out the entire `@testset` block or individual `@test` macros using `#` for single lines or `#=` and `=#` for multi-line comments.
3.  **Fix and re-enable**: Once you have fixed the issue and your primary test is passing, remember to uncomment the tests you temporarily disabled to ensure you haven't introduced any regressions.