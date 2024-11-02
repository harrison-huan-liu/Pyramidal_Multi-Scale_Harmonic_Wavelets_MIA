function [draw_wavelet, ortho]=wavelet_analysis(NodeNum, nhub, scale, multiPhi, Phi_ave, p)
    %% Visualization of Wavelets
    draw_wavelet = zeros(NodeNum, scale);
    for i = 1:NodeNum
        for j = 1:scale
            if abs(multiPhi(i,6*(j-1)+1)) > 0.1 % 0.02 % 0.01
                draw_wavelet(i,j) = multiPhi(i,6*(j-1)+1);
            end
        end
    end

    %% orthogonality
    ortho = cell(nhub, 3*scale);
    for i = 1:nhub
        indices = [0, cumsum(p(i, :))];
        for j = 1:scale
            ortho{i, j} = multiPhi(:, indices(j)+1:indices(j+1)).' * multiPhi(:, indices(j+1)+1:indices(j+2));
            ortho{i, j+scale} = multiPhi(:, indices(j)+1:indices(j+1)).' * multiPhi(:, indices(j)+1:indices(j+1));
            ortho{i, j+scale*2} = multiPhi(:, indices(j)+1:indices(j+1)).' * Phi_ave;
        end
    end
end