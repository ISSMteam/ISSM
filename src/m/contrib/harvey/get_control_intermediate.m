function [X_all, G_all] = get_control_intermediate(dir, iter, num_controls)

    X_all = cell(1, num_controls);
    G_all = cell(1, num_controls);

    for c = 0:(num_controls-1)
        filename = sprintf([dir, 'control%d_iter_%d.bin'], c, iter);

        fid = fopen(filename, 'rb', 'ieee-le');
        if fid == -1
            error('no file: %s', filename);
        end

        s = fread(fid, 1, 'uint64');
        t = fread(fid, 1, 'uint64');

        n = s * t;

        X_flat = fread(fid, n, 'double');
        G_flat = fread(fid, n, 'double');

        fclose(fid);

        X_all{c+1} = reshape(X_flat, [s, t]);
        G_all{c+1} = reshape(G_flat, [s, t]);
    end
end
