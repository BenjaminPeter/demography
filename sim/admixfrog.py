def admixfrog_input(snps, coverage, contamination, prefix='admixfrog'):
    snps['pos'] = snps.index.astype(int)
    snps2 = snps.drop_duplicates(subset=['pos'])
    afr_cols = [col for col in snps2.columns if col.startswith("afr")] 
    asn_cols = [col for col in snps2.columns if col.startswith("asn")] 
    eur_cols = [col for col in snps2.columns if col.startswith("eur")] 
    n_afr, n_asn, n_eur = len(afr_cols), len(asn_cols), len(eur_cols)
    D = dict()
    D['chrom'] = '1'
    D['pos']=snps2.index.astype(int)
    D['map']=snps2.index / 1e6
    D['ref']='A'
    D['alt']='G'
    D["AFR_alt"] = np.sum(snps2[afr_cols], 1)
    D["AFR_ref"] = n_afr - D['AFR_alt']
    D["ALT_alt"] = snps2.neand0 + snps2.neand1
    D["CHA_alt"] = snps2.neand2 + snps2.neand3
    D["VIN_alt"] = snps2.neand4 + snps2.neand5
    D["ALT_ref"] = 2 - D['ALT_alt']
    D["CHA_ref"] = 2 - D['CHA_alt']
    D["VIN_ref"] = 2 - D['VIN_alt']
    D["PAN_alt"] = snps2.chimp0 * 2
    D["PAN_ref"] = 2 - D['PAN_alt']

    if n_asn:
        D["ASN_alt"] = np.sum(snps2[asn_cols], 1)
        D["ASN_ref"] = n_asn - D['ASN_alt']
    ref = pd.DataFrame.from_dict(D)
    ref.to_csv(f"{prefix}.panel.xz", float_format="%.5f", index=False, compression="xz")

    n_libs = len(coverage)
    assert len(contamination) == n_libs
    n_samples = int(n_eur / 2)
    libs = [f"lib{i}" for i in range(n_libs)]

    for i in range(n_samples):
        ids = eur_cols[slice(2*i, 2*(i+1))]
        admixfrog_sample(ids, ref, snps2, coverage, contamination, libs, i, prefix)
    admixfrog_sample(['neand0', 'neand1'], ref, snps2, coverage, contamination, libs,
                     'ALT', prefix)


def admixfrog_sample(ids, ref, snps2, coverage, contamination, libs, name, prefix):
    print(f"samples {i}, {ids}")
    S = []
    for cov, cont, lib in zip(coverage, contamination, libs):
        print(f'Sample{i}\tLib:{lib}\tCov:{cov}\tcont:{cont}', end="\t")
        data = ref[['chrom', 'pos', 'map']].copy()
        print(f'{list(snps2[ids].columns)}')
        data['true_alt'] = np.sum(snps2[ids],1)
        data['true_ref'] = 2 - data['true_alt']
        data['lib'] = lib

        
        cov_real = poisson.rvs(cov * (1. - cont), size=data.shape[0])
        cov_cont = poisson.rvs(cov * cont, size=data.shape[0])
        p = data['true_alt'] / (data['true_ref'] + data['true_alt'])
        p_cont = ref.ASN_alt / (ref.ASN_ref + ref.ASN_alt)
        data['ralt'] = binom.rvs(cov_real, p)
        data['rref'] = cov_real - data['ralt']

        data['calt'] = binom.rvs(cov_cont, p_cont)
        data['cref'] = cov_cont - data['calt']

        data['talt'] = data.ralt + data.calt
        data['tref'] = data.rref + data.cref

        print(f"alt:\t{lib}\t{np.mean(data['ralt']):.3f}\t{np.mean(data['calt']):.3f}\t{np.mean(data['talt']):.3f}")
        print(f"ref:\t{lib}\t{np.mean(data['rref']):.3f}\t{np.mean(data['cref']):.3f}\t{np.mean(data['tref']):.3f}")

        data = data[data.tref+data.talt>0]
        S.append(data)
    data = pd.concat(S).sort_values(['chrom', 'pos', 'map', 'lib'])
    data.to_csv(f"{prefix}.{name}.sample.xz", float_format="%.5f", index=False, compression="xz")

debug_str = ""\
"--neand 55000,25,0,0.03 "\
"--eur-ages 45000 45000 40000 40000 30000 30000 "\
"--seq-len 10_000_000 "\
"--stats true_neand "\
"--introgression neand-eur "\
"--nafrA 40 "\
"--nafrB 20 "\
"--nasn 40 "\
"--snps "\
"--admixfrog "\
"--coverage 1 1 1 1 1 1 1 1 1 1 "\
"--contamination 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 "\
"--output-prefix tmp/cont".split()
#"--vcf --vcf-ploidy 1 "\
