using PGFPlotsX
using ColorSchemes
using Colors
using Statistics
using SurrogateModelOptim

function branin_2D(x)
    out = 1*(x[2] - 5.1/(4*π^2)*x[1]^2 + 5/π*x[1] - 6)^2 + 10*(1-1/8/π)*cos(x[1]) + 10
end
search_range = [(-5.0, 10.0), (0.0, 15.0)]

# Julia color scheme
j_blue = Colors.RGB([64,99,215]./255...);
j_green = Colors.RGB([56,152,38]./255...);
j_purple = Colors.RGB([149,88,178]./255...);
j_red = Colors.RGB([203,60,51]./255...);
julia_scheme = ColorScheme([j_blue,j_purple,j_green,j_red],
               "custom", "julia colors");
j_scheme = ColorScheme([get(julia_scheme, i) for i in 0.0:0.05:1.0]);


# result = smoptimize(branin_2D, search_range;
#                     options=SurrogateModelOptim.Options(
#                     iterations=15,
#                     num_interpolants=20, #Preferably even number of added processes
#                     num_start_samples=10,
#                     infill_funcs=[:median,:std],
#                     sampling_plan_opt_gens=10_000,
#                         ));
# push!(PGFPlotsX.CUSTOM_PREAMBLE,"\\pgfplotsset{
#     % define the layers you need.
#     % (Don't forget to add `main' somewhere in that list!!)
#     layers/my layer set/.define layer set={
#         background,
#         pre main,
#         main,
#         foreground
#     }{
#         % you could state styles here which should be moved to
#         % corresponding layers, but that is not necessary here.
#         % That is why wo don't state anything here
#     },
#     % activate the newly created layer set
#     set layers=my layer set,
# }")

function funcplotX(fun,j_scheme,search_range)
    
    N = 31    
    x = range(search_range[1][1], stop = search_range[1][2], length = N)
    y = range(search_range[2][1], stop = search_range[2][2], length = N)
    grid(x,y) = [x,y]
    z = @. fun(grid(x,y'))

    points_all = [result.lhc_plan result.infill_plan]
    points_vals = [fun(x) for x in eachcol(points_all)]
    
    ax = @pgf Axis(
        {
            view = (25, 17),
            "point meta min" = (minimum(z)+(maximum(z)-minimum(z))*0.05),
            "point meta max" = (maximum(z)-(maximum(z)-minimum(z))*0.39),
            "hide axis",
            "axis background/.style={fill=none}",
        },

        Plot3(
            {
                "only marks",
                "mark=*",
                "mark size" = 0.8,
                color = "black!70",
                "on layer=background",
                "mark options" = {fill="black!70"},
            },
            Coordinates(points_all[1,:],points_all[2,:],[-95 for i in 1:length(points_vals)])
        ),

        Plot3(
            {
                surf,
                mesh,
                "on layer=foreground",
                "line width=0.8pt",
            },
            Coordinates(x,y,z)
        ),

        
    )
    @pgf push!(ax, ("julia_scheme", j_scheme[:]))    

    for (x,y,z) in zip(points_all[1,:],points_all[2,:],points_vals)
        p = @pgf Plot3(
                {
                    "no marks",
                    "dashed",
                    thin,
                    color = "black!70",
                    "on layer=pre main",
                },
                Coordinates([x,x],[y,y],[-95,z])
            )
        push!(ax,p)
    end

    return ax
end

ax = funcplotX(x->mean(result.sm_interpolant(x)),j_scheme,search_range)
pgfsave("src/assets/logo.pdf", ax;dpi=600)