#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict

import msprime as msp
import pandas as pd

from sim import stats, sim


def sample_ages(ages):
    """Generate DataFrame with simulated sample names and ages."""
    df_list = []
    for p in ages:
        df = pd.DataFrame(list(enumerate(ages[p])), columns=["name", "age"])
        df.name = p + df.name.astype(str)
        df_list.append(df)
    return pd.concat(df_list, ignore_index=True)



def simulate(neand, eurA, eurB, afr, pop_params, seq_len, debug):
    id_chimp, id_neand, id_afrA, id_afrB, id_eur, id_asn = \
        [pop_params[p]["id"] for p in pop_params]

    t_chimp, t_neand, t_afrB, t_eur, t_asn = [
        sim.years_to_gen(pop_params[p]["t_split"])
        for p in ["chimp", "neand", "afrB", "eur", "asn"]
    ]

    eur_neand = sim.GeneFlow(*neand, id_eur, id_neand)
    asn_neand = sim.GeneFlow(*neand, id_asn, id_neand)
    eurA = sim.GeneFlow(*eurA, id_eur, id_afrA)
    eurB = sim.GeneFlow(*eurB, id_eur, id_afrB)
    afr = sim.GeneFlow(*afr, id_afrA, id_afrB)

    demography = [
        # geneflow "start" (backwards in time)
        eurA.end(id_eur, id_afrA), eurA.end(id_afrA, id_eur),
        eurB.end(id_eur, id_afrB), eurB.end(id_afrB, id_eur),
        afr.end(id_afrA, id_afrB), afr.end(id_afrB, id_afrA),
        eur_neand.end(id_eur, id_neand), asn_neand.end(id_asn, id_neand),

        # geneflow "end" (backwards in time)
        eurA.start(id_eur, id_afrA), eurA.start(id_afrA, id_eur),
        eurB.start(id_eur, id_afrB), eurB.start(id_afrB, id_eur),
        afr.start(id_afrA, id_afrB), afr.start(id_afrB, id_afrA),
        eur_neand.start(id_eur, id_neand), asn_neand.start(id_asn, id_neand),

        # EUR-ASN split
        msp.MassMigration(t_asn, id_asn, id_eur, 1.0),

        # population size during the bottleneck
        msp.PopulationParametersChange(
            time=t_eur - sim.years_to_gen(pop_params["eur"]["t_bottle"]),
            initial_size=pop_params["eur"]["Ne_bottle"],
            population_id=id_eur),

        # out of Africa migration
        msp.MassMigration(t_eur, id_eur, id_afrB, 1.0),

        # split of the second African lineage
        msp.MassMigration(t_afrB, id_afrB, id_afrA, 1.0),

        # split of the Neanderthal lineage from modern humans
        msp.MassMigration(t_neand, id_neand, id_afrA, 1.0),

        # split of the chimpanzee lineage from humans
        msp.MassMigration(t_chimp, id_chimp, id_afrA, 1.0),
    ]
    demography.sort(key=lambda event: event.time)

    samples = sim.define_samples(pop_params)
    pop_config = sim.define_popconfig(pop_params)

    # effective population sizes
    Ne0 = 10000

    if debug:
        print("Population identifiers:")
        for p, i in pop_params.items():
            print(f"    {p}: {i['id']}")
        msp.DemographyDebugger(
            Ne=10000,
            population_configurations=pop_config,
            demographic_events=demography
        ).print_history()
        sys.exit()
    else:
        return msp.simulate(
            length=seq_len,
            Ne=Ne0,
            mutation_rate=1e-8,
            recombination_rate=1e-8,
            samples=samples,
            population_configurations=pop_config,
            demographic_events=demography,
            record_migrations=True
        )


def geneflow_type(arg):
    args = arg.split(",")
    if len(args) != 4:
        raise argparse.ArgumentTypeError("Invalid geneflow format")
    return int(args[0]), int(args[1]), float(args[2]), float(args[3])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--eurA", help="EUR <-> AFR_A geneflow",
                        metavar="start,duration,rate1,rate2",
                        type=geneflow_type, default=(0, 0, 0, 0))
    parser.add_argument("--eurB", help="EUR <-> AFR_B geneflow",
                        metavar="start,duration,rate1,rate2",
                        type=geneflow_type, default=(0, 0, 0, 0))
    parser.add_argument("--afr", help="AFR_A <-> AFR_B geneflow",
                        metavar="start,duration,rate1,rate2",
                        type=geneflow_type, default=(0, 0, 0, 0))
    parser.add_argument("--neand", help="nonAfr <-> Neanderthal geneflow",
                        metavar="start,duration,rate1,rate2",
                        type=geneflow_type, default=(0, 0, 0, 0))

    parser.add_argument("--Ne-nonafr", type=int, default=10000)
    parser.add_argument("--Ne-afr", type=int, default=10000)

    parser.add_argument("--seq-len", help="Sequence length", type=int, required=True)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--eur-ages", nargs="+", type=int,
                        help="Ages of European samples [years BP]")
    group.add_argument("--neur", type=int,
                        help="Number of European chromosomes to simulate")
    parser.add_argument("--nasn", type=int, default=0,
                        help="Number of Asian chromosomes to simulate")
    parser.add_argument("--nafrA", type=int, default=2,
                        help="Number of AfrB chromosomes to simulate")
    parser.add_argument("--nafrB", type=int, default=2,
                        help="Number of AfrA chromosomes to simulate")

    parser.add_argument("--snps", action="store_true", help="Save all SNPs to a file")
    parser.add_argument("--stats", nargs="+",
                        choices=["true_neand", "asc_neand", "indirect", "direct", "afr_f4"],
                        help="Which statistics to calculate?")

    parser.add_argument("--output-prefix", metavar="FILE", help="Prefix of output files")

    parser.add_argument("--debug", action="store_true", help="Debug info")

    args = parser.parse_args()

    neand_ages = 2 * [125000] + 2 * [90000] + 2 * [55000]
    eur_ages = args.eur_ages if args.eur_ages else [0] * args.neur
    asn_ages = [0] * args.nasn

    pop_params = {
        "chimp": {"id": 0, "Ne": 10000, "t_sample": 1 * [0], "t_split": 6_000_000},
        "neand":  {"id": 1, "Ne": 1000,  "t_sample": neand_ages, "t_split": 600_000},
        "afrA": {"id": 2, "Ne": args.Ne_afr, "t_sample": args.nafrA * [0]},
        "afrB": {"id": 3, "Ne": args.Ne_afr, "t_sample": args.nafrB * [0], "t_split": 150_000},
        "eur": {"id": 4, "Ne": args.Ne_nonafr, "t_sample": eur_ages, "t_split": 70_000, "t_bottle": 10000, "Ne_bottle": 2000},
        "asn": {"id": 5, "Ne": args.Ne_nonafr, "t_sample": asn_ages, "t_split": 45_000},
    }

    # prepare all simulation parameters
    samples = sample_ages({"eur": eur_ages, "asn": asn_ages})

    # simulate the data
    ts = simulate(args.neand, args.eurA, args.eurB, args.afr,
                  pop_params, args.seq_len, args.debug)

    # process the simulations into different tables of SNPs
    all_snps = sim.get_all_snps(ts, sim.all_inds(pop_params))

    if args.snps:
        all_snps.to_csv(f"{args.output_prefix}_snps.tsv", sep="\t", index=False)

    if args.stats:
        # calculate admixture statistics and bind them into a DataFrame
        afrA = [f"afrA{i}" for i in range(args.nafrA)]
        afrB = [f"afrB{i}" for i in range(args.nafrB)]
        altai, chag, vindija = ["neand0", "neand1"], ["neand2", "neand3"], ["neand4", "neand5"]

        admix_snps = sim.get_asc_snps(all_snps, altai, afrA)
        true_snps = sim.get_true_snps(ts, all_snps, pop_params)

        df = defaultdict(list)
        for s in samples.name:
            if "true_neand" in args.stats: df["true_neand"].append(true_snps[s].mean())
            if "asc_neand"  in args.stats: df["asc_neand"].append((admix_snps[s] == admix_snps[altai[0]]).mean())
            if "direct"     in args.stats: df["direct"].append(stats.f4ratio(all_snps, s, a=altai, b=vindija, c=afrA, o="chimp0"))
            if "indirect"   in args.stats: df["indirect"].append(1 - stats.f4ratio(all_snps, s, a=afrA, b=afrB, c=altai, o="chimp0"))
            if "afr_f4"     in args.stats: df["afr_f4"].append((1 / len(all_snps)) * stats.f4(all_snps, w="eur0", x=s, y=afrA, z="chimp0") if s != "eur0" else "NA")
        df = pd.DataFrame(df)
        df["name"] = samples.name

        final_df = pd.merge(samples, df, on="name")

        final_df.to_csv(f"{args.output_prefix}_stats.tsv", sep="\t", index=False)
