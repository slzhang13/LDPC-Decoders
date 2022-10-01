function H_dec = H_preprocessing(H)

    [M, N] = size(H);
    dc_list = full(sum(H, 2));
    dc_max = max(dc_list);
    cn_neighbor_idx = zeros(M, dc_max);

    for cn_idx = 1:M
        cn_neighbor_idx(cn_idx, 1:dc_list(cn_idx)) = find(H(cn_idx, :));
    end

    H_dec.M = int32(M);
    H_dec.N = int32(N);
    H_dec.dc_max = int32(dc_max);
    H_dec.dc_list = int32(dc_list);
    H_dec.cn_neighbor_idx = int32(cn_neighbor_idx);

end
