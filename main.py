import os
import re

from typing import Mapping, Any, Optional
from urllib.error import URLError
from urllib.request import urlopen
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd

DATA_DIR = 'data/CD_70'
UNIPROT_URL = 'https://www.uniprot.org/uniprot/'
RHODO_TIER_1 = ['rhodopsin', 'rho', 'visual_purple', 'rhodpsn', 'rhodop', 'a0a167psn1']
RHODO_TIER_2 = ['csnbad1', 'opn2', 'opsin2', 'opsin 2']
NUM_SEQ = 10
N_THREADS = 10


def get_families():
    families = []

    for file in os.listdir(DATA_DIR):
        if file.endswith('fasta'):
            families.append(file)

    return families


def get_seqs(family):
    with open(DATA_DIR + '/' + family, 'r') as f:
        seqs = re.findall('>(.*?)/', f.read())
    return seqs


def get_page(url):
    try:
        url = urlopen(url)
    except URLError:
        pass
    yield url.read().decode('utf-8')


def write_to_file(df: pd.DataFrame) -> None:
    if df is None:
        return None

    with open('family_check_results.csv', 'w') as f:
        df.to_csv(f, header=False)


class BioVerifier:

    @staticmethod
    def check_single_family(family: str) -> Optional[Mapping[str, Any]]:

        seq_tier_1, seq_tier_2, seq_if_reviewed = ([] for _ in range(3))

        seqs = get_seqs(family)
        if seqs is None:
            return None

        print(f'Scratching around {family}, found {len(seqs)} sequences.')
        for seq in tqdm(seqs):
            fam_url = UNIPROT_URL + seq + '.txt'
            for line in get_page(fam_url):
                res_tier_1 = [x for x in RHODO_TIER_1 if (x in line.lower())]
                res_tier_2 = [x for x in RHODO_TIER_2 if (x in line.lower())]
                if 'reviewed' in line:
                    if_reviewed = 1
                elif 'unreviewed':
                    if_reviewed = 0

            seq_tier_1.append(len(res_tier_1))
            seq_tier_2.append(len(res_tier_2))
            seq_if_reviewed.append(if_reviewed)

            if len(seq_tier_1) >= NUM_SEQ:
                break

        family_results = {'family_name': family,
                          'tier_1_num_avg': np.mean(seq_tier_1),
                          'tier_1_num_med': np.median(seq_tier_1),
                          'tier_2_num_avg': np.mean(seq_tier_2),
                          'tier_2_num_med': np.median(seq_tier_2),
                          'verified_%': sum(seq_if_reviewed) / len(seqs)}


        return family_results


    @staticmethod
    def check_all_families_pipeline() -> Mapping[str, pd.DataFrame]:
        results = [BioVerifier.check_single_family(family) for family in get_families()]
        if results is not None:
            return pd.DataFrame(results)
        return None


if __name__ == '__main__':

    families = get_families()
    check = partial(BioVerifier.check_single_family)
    with Pool(N_THREADS) as p:
        results = p.map(check, families)
    #results = BioVerifier.check_all_families_pipeline()

    write_to_file(results)


