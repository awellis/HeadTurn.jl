using Documenter
using HeadTurn

makedocs(
    sitename = "HeadTurn",
    format = Documenter.HTML(),
    modules = [HeadTurn]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
