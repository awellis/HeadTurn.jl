publication_theme() = Theme(
    fontsize=24, font="sans",
    Axis=(xlabelsize=20, xgridstyle=:dash, ygridstyle=:dash,
        xtickalign=1, ytickalign=1, yticksize=10, xticksize=10,
        xlabelpadding=-5, ylim=(0, 5)),
    Legend=(framecolor=(:black, 0.5), bgcolor=(:white, 0.5)),
    Colorbar=(ticksize=16, tickalign=1, spinewidth=0.5))

savefig(fig, name::String) = begin
    save("figures/$name.pdf", fig)
    save("figures/$name.png", fig)
end

