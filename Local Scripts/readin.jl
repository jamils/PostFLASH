using Plots
using HDF5
using Unitful
using UnitfulRecipes
using ProgressBars
using DataFrames
using GLM

# ===============
# Reading in Data
# ===============

shotval = 16;
item = "ni";

shotstring = string(shotval);
path = ".\\shot" * shotstring * "\\";
file = item*".h5";
group = item*"group/"*item;
tdat = item*"group/time";
pdat = item*"group/params";

data = h5read(path*file, group);
timesteps = h5read(path*file, tdat);
params = h5read(path*file, pdat);
xmax, ymax = params;
tmax = timesteps[end];

# ===============
# Setting units
# ===============

nex, ney, net = size(data);
timestamps = 1:net;
xlen = Base.range(0, xmax, length=nex);
ylen = Base.range(0, ymax, length=ney);
tlen = Base.range(0, tmax, length=net);

xlen = @. xlen*u"cm";
ylen = @. ylen*u"cm";
tlen = @. tlen*u"s";

xlen = @. u"mm"(xlen);
ylen = @. u"mm"(ylen);

# ===============
# Computing Gradient
# ===============

function gradient(ℵ)
    len, time = size(ℵ);
    ℶ = zeros(len-1, time);

    for j ∈ 1:len-1
        ℶ[j,:] = @. ℵ[j+1,:] - ℵ[j,:];
    end

    return ℶ
end

graddat = gradient(data[begin,:,:]);

shock_locations = zeros(net);

for i ∈ 1:net
    dmp, shock_locations[i] = findmin(graddat[:,i]);
end

shock_locations = @. Int64(shock_locations);

# xrange = tlen[88:146];
# yrange = ylen[shock_locations[88:146]];


# xrange = tlen[80:140];
# yrange = ylen[shock_locations[80:140]];

xrange = tlen[71:200];
yrange = ylen[shock_locations[71:200]];

xrange = @. u"ns"(xrange); xrangeu = ustrip(xrange);
yrange = @. u"μm"(yrange); yrangeu = ustrip(yrange);

# Compute slope
data = DataFrame(y = yrangeu, x = xrangeu);
slope = lm(@formula(y ~ x), data);

b, m = coef(slope);

f(x) = @. (m*x + b) * u"μm";

slp = m * u"μm/ns";