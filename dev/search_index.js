var documenterSearchIndex = {"docs":
[{"location":"91-developer/#dev_docs","page":"Developer documentation","title":"Developer documentation","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"note: Contributing guidelines\nIf you haven't, please read the Contributing guidelines first.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If you want to make contributions to this package that involves code, then this guide is for you.","category":"page"},{"location":"91-developer/#First-time-clone","page":"Developer documentation","title":"First time clone","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: If you have writing rights\nIf you have writing rights, you don't have to fork. Instead, simply clone and skip ahead. Whenever upstream is mentioned, use origin instead.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If this is the first time you work with this repository, follow the instructions below to clone the repository.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fork this repo\nClone your repo (this will create a git remote called origin)\nAdd this repo as a remote:\ngit remote add upstream https://github.com/hz-xiaxz/KagomeDSL.jl","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"This will ensure that you have two remotes in your git: origin and upstream. You will create branches and push to origin, and you will fetch and update your local main branch from upstream.","category":"page"},{"location":"91-developer/#Linting-and-formatting","page":"Developer documentation","title":"Linting and formatting","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Install a plugin on your editor to use EditorConfig. This will ensure that your editor is configured with important formatting settings.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We use https://pre-commit.com to run the linters and formatters. In particular, the Julia code is formatted using JuliaFormatter.jl, so please install it globally first:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # Press ]\npkg> activate\npkg> add JuliaFormatter","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To install pre-commit, we recommend using pipx as follows:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"# Install pipx following the link\npipx install pre-commit","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"With pre-commit installed, activate it as a pre-commit hook:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit install","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To run the linting and formatting manually, enter the command below:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit run -a","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Now, you can only commit if all the pre-commit tests pass.","category":"page"},{"location":"91-developer/#Testing","page":"Developer documentation","title":"Testing","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"As with most Julia packages, you can just open Julia in the repository folder, activate the environment, and run test:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # press ]\npkg> activate .\npkg> test","category":"page"},{"location":"91-developer/#Working-on-a-new-issue","page":"Developer documentation","title":"Working on a new issue","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We try to keep a linear history in this repo, so it is important to keep your branches up-to-date.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fetch from the remote and fast-forward your local main\ngit fetch upstream\ngit switch main\ngit merge --ff-only upstream/main\nBranch from main to address the issue (see below for naming)\ngit switch -c 42-add-answer-universe\nPush the new local branch to your personal remote repository\ngit push -u origin 42-add-answer-universe\nCreate a pull request to merge your remote branch into the org main.","category":"page"},{"location":"91-developer/#Branch-naming","page":"Developer documentation","title":"Branch naming","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If there is an associated issue, add the issue number.\nIf there is no associated issue, and the changes are small, add a prefix such as \"typo\", \"hotfix\", \"small-refactor\", according to the type of update.\nIf the changes are not small and there is no associated issue, then create the issue first, so we can properly discuss the changes.\nUse dash separated imperative wording related to the issue (e.g., 14-add-tests, 15-fix-model, 16-remove-obsolete-files).","category":"page"},{"location":"91-developer/#Commit-message","page":"Developer documentation","title":"Commit message","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Use imperative or present tense, for instance: Add feature or Fix bug.\nHave informative titles.\nWhen necessary, add a body with details.\nIf there are breaking changes, add the information to the commit message.","category":"page"},{"location":"91-developer/#Before-creating-a-pull-request","page":"Developer documentation","title":"Before creating a pull request","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: Atomic git commits\nTry to create \"atomic git commits\" (recommended reading: The Utopic Git History).","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Make sure the tests pass.\nMake sure the pre-commit tests pass.\nFetch any main updates from upstream and rebase your branch, if necessary:\ngit fetch upstream\ngit rebase upstream/main BRANCH_NAME\nThen you can open a pull request and work with the reviewer to address any issues.","category":"page"},{"location":"91-developer/#Building-and-viewing-the-documentation-locally","page":"Developer documentation","title":"Building and viewing the documentation locally","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Following the latest suggestions, we recommend using LiveServer to build the documentation. Here is how you do it:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Run julia --project=docs to open Julia in the environment of the docs.\nIf this is the first time building the docs\nPress ] to enter pkg mode\nRun pkg> dev . to use the development version of your package\nPress backspace to leave pkg mode\nRun julia> using LiveServer\nRun julia> servedocs()","category":"page"},{"location":"91-developer/#Making-a-new-release","page":"Developer documentation","title":"Making a new release","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To create a new release, you can follow these simple steps:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Create a branch release-x.y.z\nUpdate version in Project.toml\nUpdate the CHANGELOG.md:\nRename the section \"Unreleased\" to \"[x.y.z] - yyyy-mm-dd\" (i.e., version under brackets, dash, and date in ISO format)\nAdd a new section on top of it named \"Unreleased\"\nAdd a new link in the bottom for version \"x.y.z\"\nChange the \"[unreleased]\" link to use the latest version - end of line, vx.y.z ... HEAD.\nCreate a commit \"Release vx.y.z\", push, create a PR, wait for it to pass, merge the PR.\nGo back to main screen and click on the latest commit (link: https://github.com/hz-xiaxz/KagomeDSL.jl/commit/main)\nAt the bottom, write @JuliaRegistrator register","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"After that, you only need to wait and verify:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Wait for the bot to comment (should take < 1m) with a link to a RP to the registry\nFollow the link and wait for a comment on the auto-merge\nThe comment should said all is well and auto-merge should occur shortly\nAfter the merge happens, TagBot will trigger and create a new GitHub tag. Check on https://github.com/hz-xiaxz/KagomeDSL.jl/releases\nAfter the release is create, a \"docs\" GitHub action will start for the tag.\nAfter it passes, a deploy action will run.\nAfter that runs, the stable docs should be updated. Check them and look for the version number.","category":"page"},{"location":"95-reference/#reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"95-reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Modules = [KagomeDSL]","category":"page"},{"location":"95-reference/#KagomeDSL.DoubleKagome-Tuple{Float64, Int64, Int64, Tuple{Bool, Bool}}","page":"Reference","title":"KagomeDSL.DoubleKagome","text":"double triangle unit cell Kagome lattice\n\nnote n1 here is still the number of repititions of triangle in the a1 direction, so n1 is asserted to be even. The total number is n1 * n2 * 3 sites.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.MC-Tuple{AbstractDict}","page":"Reference","title":"KagomeDSL.MC","text":"MC(params::AbstractDict)\n\n\n\nCreate a Monte Carlo object\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#Carlo.init!-Tuple{MC, MCContext, AbstractDict}","page":"Reference","title":"Carlo.init!","text":"Carlo.init!(mc::MC, ctx::MCContext, params::AbstractDict)\n\n\n\nInitialize the Monte Carlo object params\n\nn1 : Int number of cells in x direction\nn2 : Int number of cells in y direction\nPBC : Tuple{Bool,2} boundary condition, e.g. (false, false)\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#Carlo.sweep!-Tuple{MC, MCContext}","page":"Reference","title":"Carlo.sweep!","text":"Carlo.sweep!(mc::MC, ctx::MCContext) -> Nothing\n\nPerform one Monte Carlo sweep for the Mott state simulation. This implements a two-spin swap update where electrons can exchange positions between different spin states at occupied sites.\n\nNote: The W matrices are re-evaluated periodically (every n_occupied sweeps) to maintain numerical stability. If the re-evaluation fails due to singular matrices, a warning is issued and the simulation continues.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.Hmat-Tuple{DoubleKagome}","page":"Reference","title":"KagomeDSL.Hmat","text":"Hmat(lat::DoubleKagome) -> Matrix{Float64}\n\nReturn the Spinor Hamiltonian matrix for a DoubleKagome lattice.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.Sz-Tuple{Int64, Vector{Int64}, Vector{Int64}}","page":"Reference","title":"KagomeDSL.Sz","text":"Sz(i::Int, kappa_up::Vector{Int}, kappa_down::Vector{Int}) -> Float64\n\nSz = 12 (f^_ f_ - f^_ f_) Calculate the z-component of spin at site i given up and down spin configurations.\n\nReturns:     +1/2 for up spin     -1/2 for down spin\n\nThrows:     ArgumentError: if site is doubly occupied or empty     BoundsError: if i is outside the valid range\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.getOL-Tuple{Carlo.AbstractMC, Vector{Int64}, Vector{Int64}}","page":"Reference","title":"KagomeDSL.getOL","text":"The observable O_L = fracxHpsi_Gxpsi_G The Hamiltonian should be the real one!\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.getxprime-Tuple{Hamiltonian, Vector{Int64}, Vector{Int64}}","page":"Reference","title":"KagomeDSL.getxprime","text":"return x = Hx  where H is the Heisenberg Hamiltonian Note x here should also be a Mott state.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.is_occupied-Tuple{Vector{Int64}, Int64}","page":"Reference","title":"KagomeDSL.is_occupied","text":"is_occupied(kappa::Vector{Int}, l::Int) -> Bool\n\nCheck if site l is occupied in the kappa configuration vector.\n\nThrows:     BoundsError: if l is outside the valid range of kappa\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.spinInteraction!-Tuple{Dict, Vector{Int64}, Vector{Int64}, Int64, Int64}","page":"Reference","title":"KagomeDSL.spinInteraction!","text":"spinInteraction!(xprime::Dict, kappa_up::Vector{Int}, kappa_down::Vector{Int}, i::Int, j::Int)\n\nCompute the spin flip term 1/2(S+i S-j + S-i S+j) contribution to xprime.\n\nThe function handles two cases:\n\nS+i S-j: when j has up spin and i has down spin\nS-i S+j: when i has up spin and j has down spin\n\nEach case contributes with coefficient 1/2.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.tilde_U-Tuple{AbstractMatrix, Vector{Int64}}","page":"Reference","title":"KagomeDSL.tilde_U","text":"tilde_U(U::AbstractMatrix, kappa::Vector{Int})\n\n\n\nCreates a tilde matrix by rearranging rows of U according to kappa indices.\n\nParameters:\n\nU: Source matrix of size (n × m)\nkappa: Vector of indices where each non-zero value l indicates that row Rl of U         should be placed at row l of the output\n\nReturns:\n\nA matrix of size (m × m) with same element type as U\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.update_W!-Tuple{AbstractMatrix}","page":"Reference","title":"KagomeDSL.update_W!","text":"update_W!(W::AbstractMatrix; l::Int, K::Int)\n\n\n\nUpdate the W matrix W_Ij = W_Ij - W_Il  W_Kl * (W_Kj - delta_lj)\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#KagomeDSL.update_W_matrices-Tuple{MC}","page":"Reference","title":"KagomeDSL.update_W_matrices","text":"update_W_matrices(mc::MC; K_up::Int, K_down::Int, l_up::Int, l_down::Int)\n\n\n\nUpdate the W matrices\n\n\n\n\n\n","category":"method"},{"location":"90-contributing/#contributing","page":"Contributing guidelines","title":"Contributing guidelines","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"First of all, thanks for the interest!","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"We welcome all kinds of contribution, including, but not limited to code, documentation, examples, configuration, issue creating, etc.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"Be polite and respectful, and follow the code of conduct.","category":"page"},{"location":"90-contributing/#Bug-reports-and-discussions","page":"Contributing guidelines","title":"Bug reports and discussions","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you think you found a bug, feel free to open an issue. Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.","category":"page"},{"location":"90-contributing/#Working-on-an-issue","page":"Contributing guidelines","title":"Working on an issue","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you found an issue that interests you, comment on that issue what your plans are. If the solution to the issue is clear, you can immediately create a pull request (see below). Otherwise, say what your proposed solution is and wait for a discussion around it.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"tip: Tip\nFeel free to ping us after a few days if there are no responses.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If your solution involves code (or something that requires running the package locally), check the developer documentation. Otherwise, you can use the GitHub interface directly to create your pull request.","category":"page"},{"location":"","page":"KagomeDSL","title":"KagomeDSL","text":"CurrentModule = KagomeDSL","category":"page"},{"location":"#KagomeDSL","page":"KagomeDSL","title":"KagomeDSL","text":"","category":"section"},{"location":"","page":"KagomeDSL","title":"KagomeDSL","text":"Documentation for KagomeDSL.","category":"page"},{"location":"#Contributors","page":"KagomeDSL","title":"Contributors","text":"","category":"section"},{"location":"","page":"KagomeDSL","title":"KagomeDSL","text":"<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->\n<!-- prettier-ignore-start -->\n<!-- markdownlint-disable -->\n\n<!-- markdownlint-restore -->\n<!-- prettier-ignore-end -->\n\n<!-- ALL-CONTRIBUTORS-LIST:END -->","category":"page"}]
}
