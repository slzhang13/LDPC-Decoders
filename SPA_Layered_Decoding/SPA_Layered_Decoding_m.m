function [vn_llr_app, cn_llr_ext, iter_termi] = SPA_Layered_Decoding_m(H_dec, vn_llr_app, cn_llr_ext, iter_max, termi_method)

    % Reference: MATLAB document of the built-in function ldpcDecode (R2021b+)

    M = H_dec.M;
    dc_list = H_dec.dc_list;
    cn_neighbor_idx = H_dec.cn_neighbor_idx;

    iter_termi = 0;

    for iter_cnt = 1:iter_max

        for m = 1:M

            A = 0;
            S = 1;

            for n = 1:dc_list(m)
                vn_llr_app(cn_neighbor_idx(m, n)) = vn_llr_app(cn_neighbor_idx(m, n)) - cn_llr_ext(m, n);
                A = A + psi_func(vn_llr_app(cn_neighbor_idx(m, n)));
                S = S * sign(vn_llr_app(cn_neighbor_idx(m, n)));
            end

            for n = 1:dc_list(m)
                Amn = A - psi_func(vn_llr_app(cn_neighbor_idx(m, n)));
                Smn = S * sign(vn_llr_app(cn_neighbor_idx(m, n)));
                cn_llr_ext(m, n) = -Smn * psi_func(Amn);
                vn_llr_app(cn_neighbor_idx(m, n)) = vn_llr_app(cn_neighbor_idx(m, n)) + cn_llr_ext(m, n);
            end

        end

        codeword = vn_llr_app < 0;
        parity_checks = zeros(M, 1);

        for m = 1:M

            parity = 0;

            for n = 1:dc_list(m)
                parity = parity + codeword(cn_neighbor_idx(m, n));
            end

            parity_checks(m) = mod(parity, 2);

            %             parity_checks(m) = mod(sum(codeword(cn_neighbor_idx(m, :))), 2);
        end

        termi_flag = 1;

        for m = 1:M

            if parity_checks(m) == 1
                termi_flag = 0;
                break;
            end

        end

        %         termi_flag = sum(parity_checks) == 0;

        if termi_flag && termi_method == "early"
            iter_termi = iter_cnt;
            break;
        end

    end

end

function y = psi_func(x)

    t = abs(tanh(x / 2));

    if t == 0
        y = -38.14;
    else
        y = log(t);
    end

end
