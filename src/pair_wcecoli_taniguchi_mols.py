'''Pair genes from Taniguchi study with wcEcoli bulk molecules'''

import math

import pandas as pd


TANIGUCHI_PATH = "assets/TableS6.xls"
GENE_IDS = "../wcEcoli/validation/ecoli/flat/geneIDs.tsv"
WCECOLI_BULK_MOLS = "wcecoli_bulk_mols.csv"


def main() -> None:
    '''Create CSV of gene-molecule pairs

    First, we use the TSV at ``GENE_ID``S to get the EcoCyc IDs for each
    gene in the Taniguchi study. Then, we find candidate bulk molecules
    in wcEcoli for each of those EcoCyc ID. A wcEcoli molecule matches
    if it conains both the EcoCyc ID and the substring ``MONOMER``.
    The resulting pairs of gene names (as listed in the
    ``TANIGUCHI_PATH`` file) and wcEcoli bulk molecule keys are written
    to an output file ``chosen.csv``.
    '''
    taniguchi = pd.read_excel(TANIGUCHI_PATH)  # type: ignore
    gene_ids = pd.read_csv(GENE_IDS, sep="\t")
    wcecoli_bulk_mols = pd.read_csv(WCECOLI_BULK_MOLS)

    gene_synonym_ecocyc_map = dict()
    for ecocyc_id, synonyms_str in zip(
            gene_ids['FrameID'], gene_ids['Names']):
        if (not isinstance(synonyms_str, str)
                and math.isnan(synonyms_str)):
            continue
        synonyms_str = synonyms_str.replace('"', '')
        synonyms_str = synonyms_str.replace('(', '')
        synonyms_str = synonyms_str.replace(')', '')
        for synonym in synonyms_str.split():
            gene_synonym_ecocyc_map[synonym] = ecocyc_id

    # Tuple of (taniguchi gene name, wcEcoli protein key)
    chosen = []

    for name in taniguchi['Gene Name']:
        ecocyc_id = gene_synonym_ecocyc_map.get(name)
        if not ecocyc_id:
            continue
        matches = []
        for wcecoli_mol in wcecoli_bulk_mols['wcEcoli_bulkMolecules']:
            if ecocyc_id in wcecoli_mol:
                matches.append(wcecoli_mol)
        if not matches:
            continue
        if len(matches) == 1 and '_RNA' in matches[0]:
            # For many genes we don't model their translation
            continue
        filtered = [
            match for match in matches
            if 'MONOMER' in match
        ]
        if not filtered:
            print('Matches {} filtered to empty.'.format(matches))
        elif len(filtered) != 1:
            print('Multiple filtered matches:', filtered)
        else:
            chosen.append((name, filtered[0]))

    chosen_taniguchi, chosen_wcecoli = zip(*chosen)
    df = pd.DataFrame({
        'taniguchi_gene_name': chosen_taniguchi,
        'wcecoli_protein_key': chosen_wcecoli,
    })
    df.to_csv('chosen.csv', index=False)


if __name__ == '__main__':
    main()
