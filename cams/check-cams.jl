import Pkg; Pkg.activate(".");

using NCDatasets
using CairoMakie
using TimeSeries
using Statistics

fldr_data = "data/cams/";

# Load cams dataset
ds = NCDataset(fldr_data*"cams_fuentes_2004-01-02-2022-11-22.nc");

x = ds["time"][:];
y = ds["BNI"][:][1,1,1,:];

ta = TimeArray(x,y);

vals = values(ta);
tempo = string.(timestamp(ta));
lentime = length(tempo);
slice_dates = range(1, lentime, step=lentime ÷ 8)

fig = Figure(resolution=(600, 400), font="CMU Serif");
ax = Axis(fig[1, 1], xlabel="", ylabel=L"W/m^2");
line1 = lines!(ax, 1:lentime, vals; color=:black, linewidth=0.85);
ax.xticks = (slice_dates, tempo[slice_dates]);
ax.xticklabelrotation = π / 4;
ax.xticklabelalign = (:right, :center);
fig


kk = findall(y .!= 0);
hist(y[kk])