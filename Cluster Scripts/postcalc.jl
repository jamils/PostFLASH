# @everywhere using ProgressBars
@everywhere using SharedArrays

@everywhere include("postfunc.jl")

# shotval = 2;

# path = "/work/cascades/jamils/FLASH4.6.2/object/$shotval/";
@everywhere savepath = "/work/cascades/jamils/FLASH4.6.2/postdir/shot$shotval/"
@everywhere dir = "/work/cascades/jamils/FLASH4.6.2/object/";
@everywhere dir = dir * "shot$shotval/";

# i = 001;
# maxt = 179;

# Refinement factor
@everywhere reffactor = 16;
@everywhere xmat, ymat = matsize(dir);
@everywhere ne = SharedArray{Float64}(xmat, ymat, maxt);
@everywhere ni = SharedArray{Float64}(xmat, ymat, maxt);
@everywhere dens = SharedArray{Float64}(xmat, ymat, maxt);
@everywhere YE = SharedArray{Float64}(xmat, ymat, maxt);
@everywhere SUMY = SharedArray{Float64}(xmat, ymat, maxt);
@everywhere velx = SharedArray{Float64}(xmat, ymat, maxt);
@everywhere vely = SharedArray{Float64}(xmat, ymat, maxt);
@everywhere nemax = SharedArray{Float64}(maxt);
@everywhere timesteps = SharedArray{Float64}(maxt);
@everywhere dtsteps = SharedArray{Float64}(maxt);


@sync @distributed for i ∈ 1:maxt
    #println(i)
    path = pathcalc(dir, i);
    prep_vals = prep(path, reffactor);
    ye_mat = readin(path, "ye", reffactor, prep_vals);
    ye_temp = interp(ye_mat, prep_vals);

    if size(ye_temp') ≠ size(YE[:,:,i])
        error(i)
        YE[:,:,i] = ye_temp;
    else
        YE[:,:,i] = ye_temp'; # '
    end

    sumy_mat = readin(path, "sumy", reffactor, prep_vals);
    sumy_temp = interp(sumy_mat, prep_vals);

    velx_mat = readin(path, "velx", reffactor, prep_vals);
    velx_temp = interp(velx_mat, prep_vals);
    velx[:,:,i] = velx_temp';

    vely_mat = readin(path, "vely", reffactor, prep_vals);
    vely_temp = interp(vely_mat, prep_vals);
    vely[:,:,i] = vely_temp';

    dens_mat = readin(path, "dens", reffactor, prep_vals);
    dens_temp = interp(dens_mat, prep_vals);
    if size(dens_temp') ≠ size(dens[:,:,i])
        dens[:,:,i] = dens_temp;
    else
        dens[:,:,i] = dens_temp'; # '
    end


    # To pull simulation time and timestep
    time, dt = readtime(path);
    timesteps[i] = time;
    dtsteps[i] = dt;

    Nₐ = 6.02214076e23; # Avagadro constant

    # n_ion = @. Nₐ * SUMY * dens; # ion number density
    # nₑ = @. Nₐ * YE * dens; # electron number density
    #!!!!!!!!

    if size(ye_mat) != size(dens_mat)
        error("ye_mat is size ", ye_mat, " while dens is size ", size(dens_mat))
    end

    nₑ_mat = @. ye_mat * dens_mat;
    nₑ = interp(nₑ_mat, prep_vals);
    nₑ = @. Nₐ * nₑ;

    nᵢ_mat = @. sumy_mat * dens_mat;
    nᵢ = interp(nᵢ_mat, prep_vals);
    nᵢ = @. Nₐ * nᵢ;

    ###
    #=
    if size(nₑ') ≠ size(ne[:,:,i])
        error("This shit is breaking at time step ", i, ". size(nₑ') = ", size(nₑ'), "; size(ne[:,:,i]) = ", size(ne[:,:,1]));
    end
    =#

    if size(nₑ') ≠ size(ne[:,:,i])
        ne[:,:,i] = nₑ;
    else
        ne[:,:,i] = nₑ'; # '
    end

    ni[:,:,i] =  nᵢ';
    
    #ne[:,:,i] = nₑ'; # '
    nemax[i] = maximum(ne[4,:,i]);
end

#####
loc = zeros(maxt);

for i ∈ 1:maxt
    for j ∈ 1:ymat
        if ne[4,j,i] == nemax[i];
            loc[i] = j;
        end
    end
end

# Putting x, y, and t dimensions into array to store in h5
readpath = pathcalc(dir, 1);
xmax, ymax = readdomain(readpath);
params = [xmax ymax];

# Saving ne to hdf5 file
fid = h5open(savepath*"ne.h5", "cw");
create_group(fid, "negroup");
g = fid["negroup"];
g["ne"] = ne;
g["time"] = timesteps;
g["params"] = params;
close(fid);
#

fid = h5open(savepath*"ni.h5", "cw");
create_group(fid, "nigroup");
g = fid["nigroup"];
g["ni"] = ni;
g["time"] = timesteps;
g["params"] = params;
close(fid);

fid = h5open(savepath*"dens.h5", "cw");
create_group(fid, "densgroup");
g = fid["densgroup"];
g["dens"] = dens;
g["time"] = timesteps;
g["params"] = params;
close(fid);

fid = h5open(savepath*"velx.h5", "cw");
create_group(fid, "velxgroup");
g = fid["velxgroup"];
g["velx"] = velx;
g["time"] = timesteps;
g["params"] = params;
close(fid);

fid = h5open(savepath*"vely.h5", "cw");
create_group(fid, "velygroup");
g = fid["velygroup"];
g["vely"] = vely;
g["time"] = timesteps;
g["params"] = params;
close(fid);

#=
fid = h5open(savepath*"YE.h5", "cw");
create_group(fid, "YEgroup");
g = fid["YEgroup"];
g["YE"] = YE;
g["time"] = timesteps;
g["params"] = params;
close(fid);
=#
