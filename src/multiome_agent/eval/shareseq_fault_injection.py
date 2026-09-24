"""The cell-line-label fault the tenx-cell-ranger dataset couldn't exercise
(single sample, no cell-line structure) -- now buildable against the
shareseq-multi-cell-lines data, plus that dataset's own analog of
`fault_injection.py`'s RNA-ATAC pairing shuffle. Same `FaultMetadata`
convention as `fault_injection.py`. Per CLAUDE.md's data-handling rules for
this dataset, ground truth is stored with cell-line identities anonymized to
letters (A, B, C...), never the real strings -- even in-memory metadata
objects returned by these functions stay safe to print/log by construction,
not just "don't print them."
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from mudata import MuData

SHARESEQ_MODALITIES = ("rna", "rna_hvg", "atac")


@dataclass
class FaultMetadata:
    fault_type: str
    severity: float
    description: str
    ground_truth: dict = field(default_factory=dict)


def _anonymize(values: list[str]) -> dict[str, str]:
    return {v: chr(65 + i) for i, v in enumerate(sorted(set(values)))}


def swap_cell_line_labels(mdata: MuData, fraction: float, seed: int = 0) -> tuple[MuData, FaultMetadata]:
    """Reassign `fraction` of cells' `cell_line_name` to a different (real)
    cell line, simulating a demultiplexing/label-assignment error -- the
    underlying RNA/ATAC data is untouched, only the identity label is wrong,
    same in spirit as `fault_injection.py`'s `shuffle_rna_atac_pairing` but
    for cell identity rather than cross-modal pairing.
    """
    rng = np.random.default_rng(seed)
    true_labels = mdata.mod["rna"].obs["cell_line_name"]
    lines = sorted(true_labels.unique().tolist())
    anon = _anonymize(lines)
    n = len(true_labels)
    n_swap = int(round(fraction * n))
    swap_idx = rng.choice(n, size=n_swap, replace=False)

    new_labels = true_labels.copy()
    for i in swap_idx:
        true = true_labels.iloc[i]
        new_labels.iloc[i] = rng.choice([l for l in lines if l != true])

    faulted = mdata.copy()
    for mod in SHARESEQ_MODALITIES:
        faulted.mod[mod].obs["cell_line_name"] = new_labels.reindex(faulted.mod[mod].obs_names).values

    changes = {int(i): {"true": anon[true_labels.iloc[i]], "fake": anon[new_labels.iloc[i]]} for i in swap_idx}
    return faulted, FaultMetadata(
        fault_type="cell_line_label_swap", severity=fraction,
        description=f"{n_swap}/{n} cells' cell-line label reassigned to a different real cell line.",
        ground_truth={"swapped_cell_positions": [int(i) for i in swap_idx], "label_changes_anonymized": changes},
    )


def shuffle_rna_atac_pairing(mdata: MuData, fraction: float, seed: int = 0) -> tuple[MuData, FaultMetadata]:
    """Shareseq-multi-cell-lines-data analog of `fault_injection.shuffle_rna_atac_pairing`:
    permutes which ATAC profile is paired with which RNA barcode (`rna` and
    `rna_hvg` move together -- they're two views of the same RNA cells --
    only `atac`'s row order changes relative to them). Same reasoning as the
    tenx-cell-ranger version: per-cell/per-modality computations (QC, gene activity,
    motif deviations) don't depend on cross-modal pairing, so only
    cross-modal quantities (cluster agreement, cell-line recovery via ATAC,
    TF-motif correlation) are actually affected by this fault -- a clean
    fixed-core result can be reused directly rather than rerun from scratch.
    """
    rng = np.random.default_rng(seed)
    atac = mdata.mod["atac"].copy()
    n = atac.n_obs
    obs_names = list(atac.obs_names)

    n_shuffle = int(round(fraction * n))
    shuffle_idx = rng.choice(n, size=n_shuffle, replace=False)
    perm = np.arange(n)
    if n_shuffle > 1:
        sub = shuffle_idx.copy()
        offset = rng.integers(1, len(sub))
        perm[shuffle_idx] = np.roll(sub, offset)
    elif n_shuffle == 1:
        n_shuffle = 0  # no valid derangement partner for a single index

    atac_shuffled = atac[perm].copy()
    atac_shuffled.obs_names = obs_names
    n_actually_shuffled = int((perm != np.arange(n)).sum())

    faulted = MuData({"rna": mdata.mod["rna"].copy(), "rna_hvg": mdata.mod["rna_hvg"].copy(), "atac": atac_shuffled})
    faulted.uns.update(mdata.uns)

    return faulted, FaultMetadata(
        fault_type="shuffled_rna_atac_pairing", severity=fraction,
        description=f"{n_actually_shuffled}/{n} cells have their ATAC profile swapped with another cell's.",
        ground_truth={"n_shuffled": n_actually_shuffled},
    )
