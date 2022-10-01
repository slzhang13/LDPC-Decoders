function H_sparse = IEEE80211n(N, R)

    switch N
        case 1944

            switch R
                case 1/2
                    H = [57 -1 -1 -1 50 -1 11 -1 50 -1 79 -1 1 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                        3 -1 28 -1 0 -1 -1 -1 55 7 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                        30 -1 -1 -1 24 37 -1 -1 56 14 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1;
                        62 53 -1 -1 53 -1 -1 3 35 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1;
                        40 -1 -1 20 66 -1 -1 22 28 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1;
                        0 -1 -1 -1 8 -1 42 -1 50 -1 -1 8 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1;
                        69 79 79 -1 -1 -1 56 -1 52 -1 -1 -1 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1;
                        65 -1 -1 -1 38 57 -1 -1 72 -1 27 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1;
                        64 -1 -1 -1 14 52 -1 -1 30 -1 -1 32 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1;
                        -1 45 -1 70 0 -1 -1 -1 77 9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
                        2 56 -1 57 35 -1 -1 -1 -1 -1 12 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0;
                        24 -1 61 -1 60 -1 -1 27 51 -1 -1 16 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0];
                case 2/3
                    H = [61 75 4 63 56 -1 -1 -1 -1 -1 -1 8 -1 2 17 25 1 0 -1 -1 -1 -1 -1 -1;
                        56 74 77 20 -1 -1 -1 64 24 4 67 -1 7 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1;
                        28 21 68 10 7 14 65 -1 -1 -1 23 -1 -1 -1 75 -1 -1 -1 0 0 -1 -1 -1 -1;
                        48 38 43 78 76 -1 -1 -1 -1 5 36 -1 15 72 -1 -1 -1 -1 -1 0 0 -1 -1 -1;
                        40 2 53 25 -1 52 62 -1 20 -1 -1 44 -1 -1 -1 -1 0 -1 -1 -1 0 0 -1 -1;
                        69 23 64 10 22 -1 21 -1 -1 -1 -1 -1 68 23 29 -1 -1 -1 -1 -1 -1 0 0 -1;
                        12 0 68 20 55 61 -1 40 -1 -1 -1 52 -1 -1 -1 44 -1 -1 -1 -1 -1 -1 0 0;
                        58 8 34 64 78 -1 -1 11 78 24 -1 -1 -1 -1 -1 58 1 -1 -1 -1 -1 -1 -1 0];
                case 3/4
                    H = [48 29 28 39 9 61 -1 -1 -1 63 45 80 -1 -1 -1 37 32 22 1 0 -1 -1 -1 -1;
                        4 49 42 48 11 30 -1 -1 -1 49 17 41 37 15 -1 54 -1 -1 -1 0 0 -1 -1 -1;
                        35 76 78 51 37 35 21 -1 17 64 -1 -1 -1 59 7 -1 -1 32 -1 -1 0 0 -1 -1;
                        9 65 44 9 54 56 73 34 42 -1 -1 -1 35 -1 -1 -1 46 39 0 -1 -1 0 0 -1;
                        3 62 7 80 68 26 -1 80 55 -1 36 -1 26 -1 9 -1 72 -1 -1 -1 -1 -1 0 0;
                        26 75 33 21 69 59 3 38 -1 -1 -1 35 -1 62 36 26 -1 -1 1 -1 -1 -1 -1 0];
                case 5/6
                    H = [13 48 80 66 4 74 7 30 76 52 37 60 -1 49 73 31 74 73 23 -1 1 0 -1 -1;
                        69 63 74 56 64 77 57 65 6 16 51 -1 64 -1 68 9 48 62 54 27 -1 0 0 -1;
                        51 15 0 80 24 25 42 54 44 71 71 9 67 35 -1 58 -1 29 -1 53 0 -1 0 0;
                        16 29 36 41 44 56 59 37 50 24 -1 65 4 65 52 -1 4 -1 73 52 1 -1 -1 0];
                otherwise
                    disp('No such rate when N = 1944.');
            end

        case 1296

            switch R
                case 1/2
                    H = [40 -1 -1 -1 22 -1 49 23 43 -1 -1 -1 1 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                        50 1 -1 -1 48 35 -1 -1 13 -1 30 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                        39 50 -1 -1 4 -1 2 -1 -1 -1 -1 49 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1;
                        33 -1 -1 38 37 -1 -1 4 1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1;
                        45 -1 -1 -1 0 22 -1 -1 20 42 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1;
                        51 -1 -1 48 35 -1 -1 -1 44 -1 18 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1;
                        47 11 -1 -1 -1 17 -1 -1 51 -1 -1 -1 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1;
                        5 -1 25 -1 6 -1 45 -1 13 40 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1;
                        33 -1 -1 34 24 -1 -1 -1 23 -1 -1 46 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1;
                        1 -1 27 -1 1 -1 -1 -1 38 -1 44 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
                        -1 18 -1 -1 23 -1 -1 8 0 35 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0;
                        49 -1 17 -1 30 -1 -1 -1 34 -1 -1 19 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0];
                case 2/3
                    H = [39 31 22 43 -1 40 4 -1 11 -1 -1 50 -1 -1 -1 6 1 0 -1 -1 -1 -1 -1 -1;
                        25 52 41 2 6 -1 14 -1 34 -1 -1 -1 24 -1 37 -1 -1 0 0 -1 -1 -1 -1 -1;
                        43 31 29 0 21 -1 28 -1 -1 2 -1 -1 7 -1 17 -1 -1 -1 0 0 -1 -1 -1 -1;
                        20 33 48 -1 4 13 -1 26 -1 -1 22 -1 -1 46 42 -1 -1 -1 -1 0 0 -1 -1 -1;
                        45 7 18 51 12 25 -1 -1 -1 50 -1 -1 5 -1 -1 -1 0 -1 -1 -1 0 0 -1 -1;
                        35 40 32 16 5 -1 -1 18 -1 -1 43 51 -1 32 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
                        9 24 13 22 28 -1 -1 37 -1 -1 25 -1 -1 52 -1 13 -1 -1 -1 -1 -1 -1 0 0;
                        32 22 4 21 16 -1 -1 -1 27 28 -1 38 -1 -1 -1 8 1 -1 -1 -1 -1 -1 -1 0];
                case 3/4
                    H = [39 40 51 41 3 29 8 36 -1 14 -1 6 -1 33 -1 11 -1 4 1 0 -1 -1 -1 -1;
                        48 21 47 9 48 35 51 -1 38 -1 28 -1 34 -1 50 -1 50 -1 -1 0 0 -1 -1 -1;
                        30 39 28 42 50 39 5 17 -1 6 -1 18 -1 20 -1 15 -1 40 -1 -1 0 0 -1 -1;
                        29 0 1 43 36 30 47 -1 49 -1 47 -1 3 -1 35 -1 34 -1 0 -1 -1 0 0 -1;
                        1 32 11 23 10 44 12 7 -1 48 -1 4 -1 9 -1 17 -1 16 -1 -1 -1 -1 0 0;
                        13 7 15 47 23 16 47 -1 43 -1 29 -1 52 -1 2 -1 53 -1 1 -1 -1 -1 -1 0];
                case 5/6
                    H = [48 29 37 52 2 16 6 14 53 31 34 5 18 42 53 31 45 -1 46 52 1 0 -1 -1;
                        17 4 30 7 43 11 24 6 14 21 6 39 17 40 47 7 15 41 19 -1 -1 0 0 -1;
                        7 2 51 31 46 23 16 11 53 40 10 7 46 53 33 35 -1 25 35 38 0 -1 0 0;
                        19 48 41 1 10 7 36 47 5 29 52 52 31 10 26 6 3 2 -1 51 1 -1 -1 0];
                otherwise
                    disp('No such rate when N = 1296');
            end

        case 648

            switch R
                case 1/2
                    H = [0 -1 -1 -1 0 0 -1 -1 0 -1 -1 0 1 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                        22 0 -1 -1 17 -1 0 0 12 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1;
                        6 -1 0 -1 10 -1 -1 -1 24 -1 0 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1;
                        2 -1 -1 0 20 -1 -1 -1 25 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1;
                        23 -1 -1 -1 3 -1 -1 -1 0 -1 9 11 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1;
                        24 -1 23 1 17 -1 3 -1 10 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1;
                        25 -1 -1 -1 8 -1 -1 -1 7 18 -1 -1 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1;
                        13 24 -1 -1 0 -1 8 -1 6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1;
                        7 20 -1 16 22 10 -1 -1 23 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1;
                        11 -1 -1 -1 19 -1 -1 -1 13 -1 3 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
                        25 -1 8 -1 23 18 -1 14 9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0;
                        3 -1 -1 -1 16 -1 -1 2 25 5 -1 -1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0];
                case 2/3
                    H = [25 26 14 -1 20 -1 2 -1 4 -1 -1 8 -1 16 -1 18 1 0 -1 -1 -1 -1 -1 -1;
                        10 9 15 11 -1 0 -1 1 -1 -1 18 -1 8 -1 10 -1 -1 0 0 -1 -1 -1 -1 -1;
                        16 2 20 26 21 -1 6 -1 1 26 -1 7 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1;
                        10 13 5 0 -1 3 -1 7 -1 -1 26 -1 -1 13 -1 16 -1 -1 -1 0 0 -1 -1 -1;
                        23 14 24 -1 12 -1 19 -1 17 -1 -1 -1 20 -1 21 -1 0 -1 -1 -1 0 0 -1 -1;
                        6 22 9 20 -1 25 -1 17 -1 8 -1 14 -1 18 -1 -1 -1 -1 -1 -1 -1 0 0 -1;
                        14 23 21 11 20 -1 24 -1 18 -1 19 -1 -1 -1 -1 22 -1 -1 -1 -1 -1 -1 0 0;
                        17 11 11 20 -1 21 -1 26 -1 3 -1 -1 18 -1 26 -1 1 -1 -1 -1 -1 -1 -1 0];
                case 3/4
                    H = [16 17 22 24 9 3 14 -1 4 2 7 -1 26 -1 2 -1 21 -1 1 0 -1 -1 -1 -1;
                        25 12 12 3 3 26 6 21 -1 15 22 -1 15 -1 4 -1 -1 16 -1 0 0 -1 -1 -1;
                        25 18 26 16 22 23 9 -1 0 -1 4 -1 4 -1 8 23 11 -1 -1 -1 0 0 -1 -1;
                        9 7 0 1 17 -1 -1 7 3 -1 3 23 -1 16 -1 -1 21 -1 0 -1 -1 0 0 -1;
                        24 5 26 7 1 -1 -1 15 24 15 -1 8 -1 13 -1 13 -1 11 -1 -1 -1 -1 0 0;
                        2 2 19 14 24 1 15 19 -1 21 -1 2 -1 24 -1 3 -1 2 1 -1 -1 -1 -1 0];
                case 5/6
                    H = [17 13 8 21 9 3 18 12 10 0 4 15 19 2 5 10 26 19 13 13 1 0 -1 -1;
                        3 12 11 14 11 25 5 18 0 9 2 26 26 10 24 7 14 20 4 2 -1 0 0 -1;
                        22 16 4 3 10 21 12 5 21 14 19 5 -1 8 5 18 11 5 5 15 0 -1 0 0;
                        7 7 14 14 4 16 16 24 24 10 1 7 15 6 10 26 8 18 21 14 1 -1 -1 0];
                otherwise
                    disp('No such rate when N = 648');
            end

        otherwise
            disp('No such block-length');
    end

    block_size = N / size(H, 2);
    I_matrix = eye(block_size, 'uint8');
    M = block_size * size(H, 1);
    H_sparse = zeros(M, N, 'uint8');

    for i = 1:size(H, 1)

        for j = 1:size(H, 2)

            if H(i, j) ~= -1
                H_sparse((i - 1) * block_size + 1:i * block_size, (j - 1) * block_size + 1:j * block_size) = circshift(I_matrix, H(i, j), 2);
            end

        end

    end

    H_sparse = sparse(logical(H_sparse));
end
