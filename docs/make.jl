using Documenter, StructuresKit

makedocs(
    sitename="StructuresKit.jl",
    authors ="Cris Moen and contributors",
    pages = [
    "Home" => "index.md",
    "Primitives" => Any[
        "Beam.md",
        "BeamColumn.md",
        "Connections.md"],
    "Assembly" => "Assembly.md",
    "Analysis" => Any[
        "Solve.md",
        "InternalForces.md",
         "CrossSection.md"],
    "Codes and Standards" => Any[
        "AISIS10016.md",
        "AISIS10024.md"],
    "Visualization" => "Visualize.md",
    "Applications" => Any["PurlinDesigner.md"]
    ]
)

deploydocs(
    repo = "github.com/runtosolve/StructuresKit.jl.git",
)



