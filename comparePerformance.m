function performance = comparePerformance()
    % Number of runs
    num_runs = 10;

    % Size of the matrices for multiplication and other operations
    matrix_size = 1000;

    % Size of the vector for Hilbert transform and other operations
    vector_size = 10000000;

    % Initialize arrays to store elapsed times
    times = zeros(1, num_runs);

    for i = 1:num_runs
        % Generate random matrices
        A = rand(matrix_size);
        B = rand(matrix_size);

        % Generate random vector
        V = rand(vector_size, 1);

        % Start timing
        tic;

        % Perform matrix multiplication
        C = A * B;

        % Perform Hilbert transform
        H = hilbert(V);

        % Perform FFT (Fast Fourier Transform)
        F = fft(V);

        % Perform IFFT (Inverse Fast Fourier Transform)
        IF = ifft(F);

        % Stop timing and store elapsed time
        times(i) = toc;
    end

    % Calculate mean and standard deviation
    mean_time = mean(times);
    std_time = std(times);

    performance = struct('MeanTime', mean_time, 'StdTime', std_time);
end
