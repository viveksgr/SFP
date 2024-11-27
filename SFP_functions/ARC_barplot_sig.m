function ARC_barplot_sig(rsa_P1, rsa_Pp, sw)
    if nargin < 3
        sw = true;
    end

    S_mat = squeeze(nanmean(rsa_P1));
    S_err = squeeze(nanstd(rsa_P1))./sqrt(size(rsa_P1, 1));  % Corrected to consider dynamic number of subjects

    if sw
        figure('Position', [100, 100, 640, 480])  % Adjusted for better visibility
        hold on
    end

    ngroups = size(S_mat, 1);
    nbars = size(S_mat, 2);
    b = bar(S_mat);  % Capture handle to bar for setting colors
    colororder("earth")
    
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    x_m = [];
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
        errorbar(x, S_mat(:,i), S_err(:,i), 'k.');
        x_m = [x_m; x];
    end

    % Add significance markers based on p-values
    offsets = linspace(-0.1, 0.1, nbars); % Offset each set of stars slightly for visibility
    for i = 1:nbars
        for j = 1:ngroups
            % sig_x = x_m(i, j) + offsets(i);  % Adjust x position for each group
             sig_x = x_m(i, j); 
            sig_y = S_mat(j, i) + S_err(j, i) * 1.5;  % Place stars a bit above the error bars
            if rsa_Pp(j, i) < 0.001
                text(sig_x, sig_y, '***', 'HorizontalAlignment', 'center', 'Color', 'r')
            elseif rsa_Pp(j, i) < 0.01
                text(sig_x, sig_y, '**', 'HorizontalAlignment', 'center', 'Color', 'r')
            elseif rsa_Pp(j, i) < 0.05
                text(sig_x, sig_y, '*', 'HorizontalAlignment', 'center', 'Color', 'r')
            end
        end
    end


    c_s = {'r','g','b'}; % Data dots for subjects
    for ii = 1:size(rsa_P1,2) % For bars for perceptual, chemical and combinations
        for jj = 1:3
            plot(x_m(:,ii),squeeze(rsa_P1(jj,ii,:)),c_s{jj},'handle','off')
        end
    end


    if sw
        hold off
    end

end
