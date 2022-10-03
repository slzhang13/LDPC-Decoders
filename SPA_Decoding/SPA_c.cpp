#include <math.h>
#include <string.h>
#include "mex.h"

// Input Arguments
#define H_DEC prhs[0]
#define VN_LLR_APP prhs[1]
#define CN_LLR_EXT prhs[2]
#define ITER_MAX prhs[3]
#define TERMI_METHOD prhs[4]

// Output Arguments
#define VN_LLR_APP_OUT plhs[0]
#define CN_LLR_EXT_OUT plhs[1]
#define ITER_TERMI_OUT plhs[2]

const double eps = 1e-18;

#define SIGN(x) (2 * ((x) > 0) - 1)
#define PHI(x) (tanh((x) / 2) > eps ? -log(tanh((x) / 2)) : 38.14)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    VN_LLR_APP_OUT = mxDuplicateArray(VN_LLR_APP);
    CN_LLR_EXT_OUT = mxDuplicateArray(CN_LLR_EXT);

    int M, N, dc_max, iter_max, iter_termi = 0;
    int *dc_list, *cn_neighbor_idx;
    double *vn_llr_app, *cn_llr_ext, *vn_llr_ext;
    int *codeword, *parity_checks;

    M = (int)mxGetScalar(mxGetField(H_DEC, 0, "M"));
    N = (int)mxGetScalar(mxGetField(H_DEC, 0, "N"));
    dc_max = (int)mxGetScalar(mxGetField(H_DEC, 0, "dc_max"));
    dc_list = mxGetInt32s(mxGetField(H_DEC, 0, "dc_list"));
    cn_neighbor_idx = mxGetInt32s(mxGetField(H_DEC, 0, "cn_neighbor_idx"));

    vn_llr_app = mxGetDoubles(VN_LLR_APP_OUT);
    cn_llr_ext = mxGetDoubles(CN_LLR_EXT_OUT);
    iter_max = (int)mxGetScalar(ITER_MAX);

    int buff_len = mxGetNumberOfElements(TERMI_METHOD) + 1;
    char *termi_method = (char *)mxCalloc(buff_len, sizeof(char));
    mxGetString(TERMI_METHOD, termi_method, buff_len); // "early" or "max"

    bool termi_early = strcmp(termi_method, "early") == 0;

    codeword = (int *)mxCalloc(N, sizeof(int));
    parity_checks = (int *)mxCalloc(M, sizeof(int));

    for (int iter_cnt = 0; iter_cnt < iter_max; ++iter_cnt)
    {
        double *A_list = (double *)mxCalloc(M, sizeof(double));
        int *S_list = (int *)mxCalloc(M, sizeof(int));
        for (int i = 0; i < M; ++i)
        {
            S_list[i] = 1;
        }
        double *vn_llr_ext = (double *)mxCalloc(M * dc_max, sizeof(double));

        for (int m = 0; m < M; ++m)
        {
            for (int n = 0; n < dc_list[m]; ++n)
            {
                int vn_idx = cn_neighbor_idx[m + n * M] - 1;
                vn_llr_ext[m + n * M] = vn_llr_app[vn_idx] - cn_llr_ext[m + n * M];
                S_list[m] *= SIGN(vn_llr_ext[m + n * M]);
                A_list[m] += PHI(abs(vn_llr_ext[m + n * M]));
            }
        }
        for (int m = 0; m < M; ++m)
        {
            for (int n = 0; n < dc_list[m]; ++n)
            {
                int vn_idx = cn_neighbor_idx[m + n * M] - 1;
                vn_llr_app[vn_idx] -= cn_llr_ext[m + n * M];
            }
        }

        for (int m = 0; m < M; ++m)
        {
            for (int n = 0; n < dc_list[m]; ++n)
            {
                int vn_idx = cn_neighbor_idx[m + n * M] - 1;
                double tmp = A_list[m] - PHI(abs(vn_llr_ext[m + n * M]));
                cn_llr_ext[m + n * M] = S_list[m] * SIGN(vn_llr_ext[m + n * M]) * PHI(tmp);
                vn_llr_app[vn_idx] += cn_llr_ext[m + n * M];
            }
        }

        mxFree(A_list);
        mxFree(S_list);
        mxFree(vn_llr_ext);

        for (int n = 0; n < N; ++n)
        {
            codeword[n] = vn_llr_app[n] < 0;
        }

        for (int m = 0; m < M; ++m)
        {
            int parity = 0;
            for (int n = 0; n < dc_list[m]; ++n)
            {
                int vn_idx = cn_neighbor_idx[m + n * M] - 1;
                parity += codeword[vn_idx];
            }
            parity_checks[m] = parity & 0x1;
        }

        bool valid_codeword = true;
        for (int m = 0; m < M; ++m)
        {

            if (parity_checks[m] > 0)
            {
                valid_codeword = false;
                break;
            }
        }

        if (valid_codeword && termi_early)
        {
            iter_termi = iter_cnt + 1;
            break;
        }
    }

    mxFree(codeword);
    mxFree(parity_checks);
    mxFree(termi_method);

    ITER_TERMI_OUT = mxCreateDoubleScalar(iter_termi);
}