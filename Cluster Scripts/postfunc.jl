@everywhere using MATLAB
@everywhere using HDF5
using DataFrames
using GLM
using DelimitedFiles
using CSV

@everywhere function prep(path, reffactor)
    blksize = h5read(path, "block size");
    coord = h5read(path, "coordinates");

    bl = size(blksize,2);

    # Get center coordinates of each cell and store them to x and y 4 cells in one block

    x = zeros(reffactor*bl + reffactor);
    y = zeros(reffactor*bl + reffactor);

    for j ∈ 1:size(blksize, 2)
        for i ∈ 1:reffactor
            x[reffactor*(j-1)+i] = coord[1,j] - blksize[1,j]/2 + i*blksize[1,j]/reffactor - blksize[1,j]/(2*reffactor);
            y[reffactor*(j-1)+i] = coord[2,j] - blksize[2,j]/2 + i*blksize[2,j]/reffactor - blksize[2,j]/(2*reffactor);
        end
    end

    # Convert vector x and y to matrixes X and Y (to satisfy requirement of griddata)

    X = zeros(reffactor, reffactor*bl);
    Y = zeros(reffactor, reffactor*bl);

    for j ∈ 1:size(blksize,2)
        for i ∈ 1:reffactor
            X[i, reffactor*(j-1)+1:reffactor*j] .= x[reffactor*(j-1)+i];
            Y[:,reffactor*(j-1)+i] .= y[reffactor*(j-1)+i];
        end
    end

    # Generate interpolation grid and generate xminimum, xmax, yminimum, ymax for plot

    xres = minimum(minimum(@. abs(-X[1:end-1,:] + X[2:end,:]))); # finest mash size in simulation
    yres = minimum(minimum(@. abs(-Y[:,1:end-1] + Y[:,2:end]))); # finest mash size in simulation

    xmin = minimum(@. coord[1,:] - blksize[1,:]/2); # Boundary position, not center
    xmax = maximum(@. coord[1,:] + blksize[1,:]/2);
    ymin = minimum(@. coord[2,:] - blksize[2,:]/2);
    ymax = maximum(@. coord[2,:] + blksize[2,:]/2);

    xq = (xmin + xres/2:xres:xmax-xres/2);
    yq = (ymin + yres/2:yres:ymax-yres/2);
    yq = yq';

    return blksize, bl, x, y, X, Y, xq, yq
end

@everywhere function readin(path, var::String, reffactor, prep_vals)
    blksize, bl, x, y, X, Y, xq, yq = prep_vals;

    var_mat = zeros(reffactor, reffactor*bl);
    var_read = h5read(path, "/$var");

    for j ∈ 1:bl
        for i ∈ 1:reffactor
            var_mat[:,reffactor*(j-1) + 1:reffactor*j] = var_read[:,:,1,j];
        end
    end

    return var_mat
end

@everywhere function readtime(path)
    time = h5read(path, "/real scalars")[1][:value];
    dt = h5read(path, "/real scalars")[2][:value];
    
    return time, dt
end

@everywhere function readdomain(path)
    xmax = h5read(path, "/real runtime parameters")[593][:value];
    ymax = h5read(path, "/real runtime parameters")[596][:value];

    return xmax, ymax
end

@everywhere function pathcalc(dir,i)
    precurse = "lasslab_hdf5_plt_cnt_";
    dir = dir * precurse;
    if i<10
        path = dir * "000";
    elseif i>9 && i<100
        path = dir * "00";
    elseif i>99 && i<1000
        path = dir * "0";
    elseif i>999
        path = dir;
    end

    item = string(i);
    path = path * item;

    return path
end

@everywhere function interp(var, prep_vals)
    blksize, bl, x, y, X, Y, xq, yq = prep_vals;

    X = @. real(X);
    Y = @. real(Y);
    var = @. real(var);
    xq = @. real(xq);
    yq = @. real(yq);

    @mput X Y var xq yq

    eval_string("vr = griddata(X,Y,var,xq,yq,'nearest');");

    @mget vr

    return vr
end

@everywhere function matsize(dir)
    path = pathcalc(dir,0)
    prep_vals = prep(path, reffactor);
    ye_mat = readin(path, "ye", reffactor, prep_vals);
    YE = interp(ye_mat, prep_vals);

    xmat = length(YE[1,:]);
    ymat = length(YE[:,1]);

    return xmat, ymat
end
