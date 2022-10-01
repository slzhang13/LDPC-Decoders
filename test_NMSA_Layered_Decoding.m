clear;

addpath("NMSA_Layered_Decoding/")

% code parameters
N = 64800;
R = 2/3;
K = N * R;
M = N - K;

H = dvbs2ldpc(R); % 64800
% H = IEEE80211n(N, R);

% built-in codec configuration
enc_cfg = ldpcEncoderConfig(H);
dec_cfg = ldpcDecoderConfig(H, 'norm-min-sum');

termi_method = 'early';
% termi_method = 'max';

iter_max = 10;
iter_max = int32(iter_max); % int32

rng(12345); % reproducibility

% info bits and coded bits
b = randi([0, 1], K, 1);
c = ldpcEncode(b, enc_cfg);

% BPSK
s = 1 - 2 * c;

% AWGN
snr_db = 4.0;
snr = 10^(0.1 * snr_db);
sigma2 = 1 / snr;
noise = randn(size(s)) * sqrt(sigma2);

y = s + noise;

% soft infomation from channel can be easily calculated for BPSK/QPSK
Lch = 2 * y / sigma2;

alpha = 0.7;

% preprocess PCM
H_dec = H_preprocessing(H);

% initialization
% vn_llr_app_m = Lch;
% vn_llr_app_c = Lch;
cn_llr_ext_m = zeros(M, H_dec.dc_max);
% cn_llr_ext_c = zeros(M, H_dec.dc_max);

tic
[vn_llr_app_m, cn_llr_ext_m, iter_termi_m] = NMSA_Layered_Decoding_m(H_dec, Lch, cn_llr_ext_m, iter_max, termi_method, alpha);
toc

tic
vn_llr_app_matlab = ldpcDecode(Lch, dec_cfg, iter_max, 'OutputFormat', 'whole', 'DecisionType', 'soft', 'Termination', termi_method, 'MinSumScalingFactor', alpha);
toc

norm(vn_llr_app_matlab - vn_llr_app_m)