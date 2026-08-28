import numpy as np

def load(fname):
    return np.genfromtxt(fname, skip_header=1)

# map delay tags (as printed in filenames) to the matching reference trSPEC file (by tm index)
delay_map = [
    ('-5.654867', 'trSPEC_XY_0010.dat'),
    ('-3.769911', 'trSPEC_XY_0310.dat'),
    ('-1.884956', 'trSPEC_XY_0510.dat'),
    ('0.000000',  'trSPEC_XY_0910.dat'),
    ('1.884956',  'trSPEC_XY_1310.dat'),
    ('3.769911',  'trSPEC_XY_1510.dat'),
    ('5.654867',  'trSPEC_XY_1810.dat'),
]

for tag, ref_file in delay_map:
    ref = load(f'/tmp/fortran_toy_reference/{ref_file}')
    ch0 = load(f'/tmp/modeb_toy_test2/trRIXS2_chan0_t{tag}.txt')
    ch1 = load(f'/tmp/modeb_toy_test2/trRIXS2_chan1_t{tag}.txt')
    ch2 = load(f'/tmp/modeb_toy_test2/trRIXS2_chan2_t{tag}.txt')
    total = ch0[:,2] + ch1[:,2] + ch2[:,2]
    ref_peak_idx = ref[:,2].argmax()
    mine_peak_idx = total.argmax()
    print(f"=== delay={tag} ===")
    print("  ref : peak I=%.6e at ab=%.3f loss=%.2f" % (ref[ref_peak_idx,2], ref[ref_peak_idx,0], ref[ref_peak_idx,1]))
    print("  mine: peak I=%.6e at ab=%.3f loss=%.2f" % (total[mine_peak_idx], ch0[mine_peak_idx,0], ch0[mine_peak_idx,1]))
    # also print full loss profile at the reference's resonant ab, for both
    ref_ab = ref[ref_peak_idx,0]
    ref_sub = ref[np.isclose(ref[:,0], ref_ab, atol=0.01)]
    ref_sub = ref_sub[np.argsort(ref_sub[:,1])]
    mine_mask = np.isclose(ch0[:,0], ref_ab, atol=0.01)
    mine_sub_loss = ch0[mine_mask,1]
    mine_sub_I = total[mine_mask]
    order = np.argsort(mine_sub_loss)
    print("  loss-axis @ ab=%.3f:" % ref_ab)
    for i in range(len(ref_sub)):
        print("    loss=%5.2f  ref=%.4e  mine=%.4e" % (ref_sub[i,1], ref_sub[i,2], mine_sub_I[order][i]))
    print()
