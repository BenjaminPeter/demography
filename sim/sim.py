import random
from collections import defaultdict

import pandas as pd
import numpy as np
import msprime as msp

from . import stats

def years_to_gen(years, gen_time=25):
    return int(years / gen_time)


def get_all_snps(ts, samples):
    """Extract all simulated SNPs from the tree sequence replicates."""
    if isinstance(ts, msp.TreeSequence):
        ts = [ts]

    return pd.concat(
        pd.DataFrame(
            t.genotype_matrix(),
            columns=samples,
            index=(v.position for v in t.variants())
        ) for t in ts
    )


def get_asc_snps(snps, neand, afr):
    """Filter for fixed differences between Africans and Neandertals."""
    afr_freq = stats.calc_freq(snps, afr)
    neand_freq = stats.calc_freq(snps, neand)

    return snps.loc[((afr_freq == 0.0) & (neand_freq == 1.0)) | ((afr_freq == 1.0) & (neand_freq == 0.0))]


def get_true_snps(ts, all_snps, pop_params):
    """Examine coalescent trees and find true Neanderthal-derived SNPs."""
    neand_names = pop_inds("neand", pop_params)

    mut_pops = assign_pops(ts)
    neand_snps = all_snps[mut_pops == pop_params["neand"]["id"]]
    fixed_neand_snps = neand_snps[neand_snps[neand_names].mean(axis=1) == 1.0]

    return fixed_neand_snps


def assign_times(ts):
    """Randomly assign time of origin to each mutation.
    Inspired by an example from msprime's documentation."""
    rng = random.Random()
    mut_times = np.zeros(ts.num_mutations)
    for tree in ts.trees():
        for mutation in tree.mutations():
            a = tree.time(mutation.node)
            b = tree.time(tree.parent(mutation.node))
            mut_times[mutation.id] = rng.uniform(a, b)
    return mut_times


def gather_migrations(ts):
    """Gather all migrations in each node in a tree."""
    node_migrations = defaultdict(list)
    for migration in ts.migrations():
        node_migrations[migration.node].append(migration)
    return node_migrations


def assign_pops(ts):
    """Assign population of origin to each mutation.
    Inspired by an example from msprime's documentation."""
    mut_times = assign_times(ts)
    node_migrations = gather_migrations(ts)

    pop_assign = np.repeat(-1, ts.num_mutations)
    for tree in ts.trees():
        for site in tree.sites():
            for mut in site.mutations:
                pop_assign[mut.id] = tree.population(mut.node)
                for mig in node_migrations[mut.node]:
                    if mig.left <= site.position < mig.right:
                        if mig.time < mut_times[mut.id]:
                            assert pop_assign[mut.id] == mig.source
                            pop_assign[mut.id] = mig.dest
    return pop_assign


def pop_inds(pop, pop_params):
    """Generate list of string names of population samples
    (i.e. ["eur0", "eur1",..., "eurN"]).
    """
    names = all_inds(pop_params)
    return [name for name in names if name.startswith(pop)]


def all_inds(pop_params):
    """Generate list of all simulated sample names."""
    sample_names = []
    for p in pop_params:
        n_pop = len(pop_params[p]["t_sample"])
        sample_names.extend([f"{p}{i}" for i in range(n_pop)])
    return sample_names


def define_samples(pop_params):
    """Generate list of sample definitions for msprime."""
    sample_names = []
    for i, pop in enumerate(pop_params):
        times = [years_to_gen(t) for t in pop_params[pop]["t_sample"]]
        sample_names.extend([msp.Sample(population=i, time=t) for t in times])
    return sample_names


def define_popconfig(pop_params):
    """Generate list of population configurations for msprime."""
    return [msp.PopulationConfiguration(
        initial_size=pop_params[p]["Ne"]) for p in pop_params]


class GeneFlow:
    def __init__(self, start, duration, prop1, prop2, pop1, pop2):
        """Start - "real" start forward in time."""
        self._start = years_to_gen(start)
        self._duration = years_to_gen(duration)
        self._end = self._start - self._duration
        self._populations = pop1, pop2
        self._rates = {
            (pop1, pop2): prop2 / self._duration if self._duration else 0,
            (pop2, pop1): prop1 / self._duration if self._duration else 0
        }

    def __str__(self):
        return f"""Geneflow properties
    start: {self._end} generations BP
    duration: {self._duration} generations
    populations: {self._populations}
    ancestry rates: {self._rates.values()}"""

    def _geneflow(self, pop1, pop2, start, rate):
        return msp.MigrationRateChange(
            time=start,
            rate=rate,
            matrix_index=(pop1, pop2)
        )

    def end(self, pop1, pop2):
        if len(set([pop1, pop2]) & set(self._populations)) != 2:
            raise TypeError("Invalid geneflow population pair")
        return self._geneflow(pop1, pop2, self._end, self._rates[(pop1, pop2)])

    def start(self, pop1, pop2):
        if len(set([pop1, pop2]) & set(self._populations)) != 2:
            raise TypeError("Invalid geneflow population pair")
        return self._geneflow(pop1, pop2, self._start, 0.0)
