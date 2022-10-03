function [vn_llr_app, cn_llr_ext, iter_termi] = SPA_m(H_dec, vn_llr_app, cn_llr_ext, iter_max, termi_method)

    % Reference: Channel Codes Classical and Modern (Sec. 5.4)

    M = H_dec.M;
    dc_list = H_dec.dc_list;
    dc_max = H_dec.dc_max;
    cn_neighbor_idx = H_dec.cn_neighbor_idx;

    iter_termi = 0;

    for iter_cnt = 1:iter_max

        vn_llr_app_bak = vn_llr_app;

        vn_llr_ext = zeros(M, dc_max);
        S_list = ones(M, 1);
        A_list = zeros(M, 1);

        for m = 1:M

            for n = 1:dc_list(m)
                vn_idx = cn_neighbor_idx(m, n);
                vn_llr_ext(m, n) = vn_llr_app_bak(vn_idx) - cn_llr_ext(m, n);
                S_list(m) = S_list(m) * sign(vn_llr_ext(m, n));
                A_list(m) = A_list(m) + phi_func(abs(vn_llr_ext(m, n)));
                vn_llr_app(vn_idx) = vn_llr_app(vn_idx) - cn_llr_ext(m, n);
            end

        end

        for m = 1:M

            for n = 1:dc_list(m)
                vn_idx = cn_neighbor_idx(m, n);
                tmp = A_list(m) - phi_func(abs(vn_llr_ext(m, n)));
                cn_llr_ext(m, n) = S_list(m) * sign(vn_llr_ext(m, n)) * phi_func(tmp);
                vn_llr_app(vn_idx) = vn_llr_app(vn_idx) + cn_llr_ext(m, n);
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

function y = phi_func(x)

    t = tanh(x / 2);

    if t == 0
        y = 38.14;
    else
        y = -log(t);
    end

end
