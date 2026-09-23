"""The two fault types the tenx-cell-ranger dataset couldn't exercise (single
sample, no cell-line/library structure) -- now buildable against the
shareseq-multi-cell-lines data. Same `FaultMetadata` convention as
`fault_injection.py`. Per CLAUDE.md's data-handling rules for this dataset,
ground truth is stored with cell-line/library identities anonymized to
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


def mix_samples(mdata: MuData, seed: int = 0) -> tuple[MuData, FaultMetadata]:
    """Pick two distinct real sequencing libraries and relabel one's cells
    to claim the other's library identity -- simulating two independent
    libraries being merged/mislabeled as a single sample. A genuinely
    different failure mode from `swap_cell_line_labels`: per-cell
    genotype-based cell-line identity is untouched here; what's wrong is
    which cells count as belonging to the same technical batch/sample.
    """
    rng = np.random.default_rng(seed)
    library_col = mdata.mod["rna"].obs["library"]
    libraries = sorted(library_col.unique().tolist())
    lib_a, lib_b = rng.choice(libraries, size=2, replace=False)
    anon = _anonymize([lib_a, lib_b])

    is_b = (library_col == lib_b).to_numpy()
    faulted = mdata.copy()
    for mod in SHARESEQ_MODALITIES:
        obs = faulted.mod[mod].obs
        mask = (obs["library"] == lib_b).to_numpy()
        obs.loc[mask, "library"] = lib_a
        obs.loc[mask, "sample"] = lib_a

    n_mixed = int(is_b.sum())
    return faulted, FaultMetadata(
        fault_type="mixed_samples", severity=n_mixed / mdata.n_obs,
        description=(
            f"{n_mixed} cells from a second sequencing library (anonymized '{anon[lib_b]}') "
            f"relabeled to appear part of a different sample ('{anon[lib_a]}')."
        ),
        ground_truth={
            "mixed_cell_positions": [int(i) for i in np.where(is_b)[0]],
            "true_library_anonymized": {"kept_identity": anon[lib_a], "mislabeled_into_it": anon[lib_b]},
        },
    )
