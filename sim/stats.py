def calc_freq(snps, sample_names=None):
    """Calculate the frequencies of SNPs."""
    if not sample_names:
        sample_names = list(snps.columns)

    if not isinstance(sample_names, list):
        sample_names = [sample_names]

    allele_counts = snps.loc[:, sample_names].sum(axis=1)
    return allele_counts / len(sample_names)


def f4(snps, w, x, y, z):
    w_freq = calc_freq(snps, w)
    x_freq = calc_freq(snps, x)
    y_freq = calc_freq(snps, y)
    z_freq = calc_freq(snps, z)

    return ((w_freq - x_freq) * (y_freq - z_freq)).sum()


def f4ratio(snps, x, a, b, c, o):
    return f4(snps, a, o, x, c) / f4(snps, a, o, b, c)


def d(snps, w, x, y, z):
    w_freq = calc_freq(snps, w)
    x_freq = calc_freq(snps, x)
    y_freq = calc_freq(snps, y)
    z_freq = calc_freq(snps, z)

    nom = ((w_freq - x_freq) * (y_freq - z_freq)).sum()
    denom = ((w_freq + x_freq - 2 * w_freq * x_freq) *
             (y_freq + z_freq - 2 * y_freq * z_freq)).sum()

    return nom / denom
