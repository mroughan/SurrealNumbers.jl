# parameters used in common in plots
#   https://www.matecdev.com/posts/julia-plotting-font-size.html
#   https://docs.juliaplots.org/latest/generated/attributes_plot/
#   https://docs.juliaplots.org/latest/generated/attributes_series/
#   https://docs.juliaplots.org/latest/generated/colorschemes/#scientific
# 
using Colors
using ColorSchemes
# colors = ColorSchemes.darktest.colors
# colors = palette(:Dark2_8)
# colors = palette(:Set1_9)
# color_palette = palette(:seaborn_dark6)
color_palette = palette(:seaborn_colorblind)
default(color_palette=color_palette)

using LaTeXStrings 
latex_example = L"\sin(2\pi x)"

lw = 1.5 # default linewidth
default(linewidth=lw)

ms = 8   # default markersize
sms = 4  # size for small markers
ma = 0.05 # alpha for faint markers
default(markersize=ms)

marker = [:circle, :cross, :square]
mc = palette(:seaborn_muted)

# fontsizes = (
#     titlefontsize = 200,# plot title (should not use)
#     tickfontsize=16,    # axis tick labels
#     guidefontsize=18,   # axis labels
#     legendfontsize=32,  # legend    
#     legendtitlefontsize=200,
# )
# colorbar_tickfontsize
Plots.scalefontsizes()
Plots.scalefontsizes(1.8) # default font is quite a bit too small for default figsize

# default legend location
default(legend=:topright)
# :outertopright


# line=(3,:green,:dash,:sticks)  :solid, :dash, :dot, :dashdot
# marker=(:circle,8,:green,:green) 
#    markeralpha (interior)
#    markershape [:none, :auto, :circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x].
#    
# seriestype =  [:none, :line, :path, :steppre, :stepmid, :steppost, :sticks, :scatter, :heatmap, :hexbin, :barbins, :barhist, :histogram, :scatterbins, :scatterhist, :stepbins, :stephist, :bins2d, :histogram2d, :histogram3d, :density, :bar, :hline, :vline, :contour, :pie, :shape, :image, :path3d, :scatter3d, :surface, :wireframe, :contour3d, :volume, :mesh3d] 
# 
# 

# default(titlefont = (20, "times"), legendfontsize = 18, guidefont = (18, :darkgreen), tickfont = (12, :orange), guide = "x", framestyle = :zerolines, yminorgrid = true)
# default(guidefont = (18, :darkgreen), tickfont = (12, :orange), guide = "x", color_palette=color_palette, linewidth=lw)
