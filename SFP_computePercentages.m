function output = SFP_computePercentages(rsa_P1_, thr)
    [rows, cols, ~] = size(rsa_P1_);  % Get the size of the cell array
    output = zeros(rows, cols, 3);    % Initialize the output array

    for i = 1:rows
        for j = 1:cols
            % Extract the arrays for the current cell
            array1 = rsa_P1_{i, j, 1};
            array2 = rsa_P1_{i, j, 2};

            % Calculate the percentages
            % % (1) Elements in array1 that exceed thr but not in array2
            % count1 = sum(array1 > thr & array2 <= thr);
            % output(i, j, 1) = (count1 / numel(array1)) * 100;
            % 
            % % (2) Elements in array2 that exceed thr but not in array1
            % count2 = sum(array2 > thr & array1 <= thr);
            % output(i, j, 2) = (count2 / numel(array2)) * 100;
            % 
            % % (3) Elements that are higher than thr in both arrays
            % count3 = sum(array1 > thr & array2 > thr);
            % output(i, j, 3) = (count3 / numel(array1)) * 100; % Assuming array1 and array2 have the same size


            count1 = sum(array1 - array2 > thr);
            output(i, j, 1) = (count1 / numel(array1)) * 100;

            % (2) Elements in array2 that exceed thr but not in array1
            count2 = sum(array2 - array1 > thr);
            output(i, j, 2) = (count2 / numel(array2)) * 100;

            % (3) Elements that are higher than thr in both arrays
            count3 = sum(abs(array1-array2)< thr);
            output(i, j, 3) = (count3 / numel(array1)) * 100; % Assuming array1 and array2 have the same size


        end
    end
end
