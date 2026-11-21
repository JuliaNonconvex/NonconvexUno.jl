using DocumenterFragments: build_fragment

# The documentation is built as a DocumenterFragments fragment: it assembles
# into the org-wide Nonconvex.jl documentation site. Links to sibling packages
# (e.g. `../core/...`) only resolve in that assembled site, so tolerate broken
# local cross-references in this standalone preview build (which is also what
# the repository's GitHub Pages deploy serves).
build_fragment(@__DIR__; checkdocs=:all, doctest=true, warnonly=[:cross_references])
