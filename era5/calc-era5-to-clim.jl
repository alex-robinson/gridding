import Pkg; Pkg.activate("$(homedir())/.JuliaEnvironments/myanalysis");

using NCDatasets
using DataStructures
using Glob
using JSON 

function calc_multifile_average(files,var_name)

    # Calculate the average of the given variable
    nc = NCDataset(files[1]);
    var = nc[var_name][:];
    close(nc);

    for k = 2:length(files)
        nc = NCDataset(files[k]);
        now = nc[var_name][:];
        var .= var .+ now;
        close(nc);
    end

    var .= var ./ length(files);

    return var
end

function calc_clim_dataset(var_file_name;clim_range = (1961,1990),pres=nothing)

    # Generate correct filename for output
    clim_range_str = string.(clim_range);

    
    if isnothing(pres)
        var_file_name_now = var_file_name;
    else
        var_file_name_now = var_file_name*"_"*string(pres);
    end

    file_out = fldr_out*"era5_monthly-single-levels_"*var_file_name_now*"_"*clim_range_str[1]*"-"*clim_range_str[2]*".nc";


    #var_name  = "t2m";
    #long_name = "2m_temperature";

    # Get list of relevant files for this variable
    files = glob("era5_monthly-single-levels_"*var_file_name_now*"*.nc",fldr_data);

    # Define the variable name of interest 
    nc = NCDataset(files[1]);
    var_names_all = keys(nc);
    k = findall( .! in(["longitude","latitude","time"]).(var_names_all) )
    var_name = var_names_all[k][1];
    println("var_name = $var_name")

    # Also get dimension information 
    lon  = nc["longitude"][:];
    lat  = nc["latitude"][:];
    time = nc["time"][:];

    close(nc); 

    # Calculate the average over multiple files
    var = calc_multifile_average(files,var_name);

    # Generate code to make a new NetCDF Dataset, and include it here
    ncgen(files[1],"tmp.jl",newfname=file_out);
    tmp  = read("tmp.jl",String);
    
    # Modify any defintions of Int16 to Int32 to avoid errors
    tmp = replace(tmp,"Int16" => "Int32");

    # Modify definitions or variables with periods in the name
    tmp = replace(tmp,"ncp54.162" => "ncp54162");
    tmp = replace(tmp,"ncp63.162" => "ncp63162");
    tmp = replace(tmp,"ncp53.162" => "ncp53162");

    f = open("tmp.jl", "w");
    write(f, tmp);
    close(f);

    # Load the script that will generate the new NetCDF file
    include("tmp.jl");

    # Load the new NetCDF Dataset file
    ds = NCDataset(file_out,"a");

    # Define actual variable values in new file
    ds["longitude"][:]  = lon;
    ds["latitude"][:]   = lat;
    ds["time"][:]       = time;
    ds[var_name][:]     = var;

    close(ds);

    println("Produced file: ", file_out);

end

# Load information about all variables of interest 
info      = JSON.parsefile("era5_config_monthly.json");
fldr_data = "data/era5/monthly-single-levels/";
fldr_out  = "data/era5/monthly-single-levels/clim/";

begin
if true
    # Are we working with a specific pressure level (else `nothing`)
    pres = nothing;
    #vars_all = info["vars"];
    vars_all = info["vars"][37:39];

    #pres = 750;
    #vars_all = info["vars_pres"];

    # Define climatology range of interest
    clim_range = (1961,1990);
    #clim_rante = (1981,2010);

    # Loop over all filename variables
    for var_file_name in vars_all
        calc_clim_dataset(var_file_name;clim_range=clim_range,pres=pres)
    end
end
end

# Take a look at the final file...
#ds = NCDataset(file_out);
