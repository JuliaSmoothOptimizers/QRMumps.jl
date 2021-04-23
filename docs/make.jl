using Documenter, QRMumps

makedocs(
  modules = [QRMumps],
  checkdocs = :exports,
  doctest = true,
  linkcheck = true,
  strict = false,
  format = Documenter.HTML(assets = ["assets/style.css"], prettyurls = get(ENV, "CI", nothing) == "true"),
  sitename = "QRMumps.jl",
  pages = ["Introduction" => "index.md",
           "Features" => "features.md",
           "Optional features" => "optional_features.md",
           "API" => "api.md",
           "Control parameters" => "control_parameters.md",
           "Information parameters" => "information_parameters.md",
           "Performance tuning" => "performance.md",
           "Tutorials" => ["Symmetric and positive definite linear systems" => "tutorials/spd.md",
                           "Least-squares problems" => "tutorials/ls.md",
                           "Least-norm problems" => "tutorials/ln.md"],
           "Reference" => "reference.md",
          ]
)

deploydocs(repo = "github.com/JuliaSmoothOptimizers/QRMumps.jl.git")
