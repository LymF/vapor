"""Pure grouping resolver for the co-assembly track. Snakemake-free (pytest-able)."""
import csv
import os


def parse_groups(metadata_path: str, samples: list, mode: str) -> dict:
    mode = (mode or "metadata").strip().lower()
    sample_set = set(samples)

    if mode == "all":
        return {"all": list(samples)}

    # mode == "metadata"
    if not metadata_path or not os.path.exists(metadata_path):
        raise ValueError(
            f"coassembly.grouping=metadata requer um metadata TSV existente "
            f"(coassembly.metadata). Caminho inválido: {metadata_path!r}"
        )
    groups: dict = {}
    with open(metadata_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        if "sample" not in cols or "group" not in cols:
            raise ValueError(
                "metadata TSV precisa das colunas 'sample' e 'group'; "
                f"encontrado: {reader.fieldnames}"
            )
        for row in reader:
            s = (row[cols["sample"]] or "").strip()
            g = (row[cols["group"]] or "").strip()
            if not s or not g or s not in sample_set:
                continue
            groups.setdefault(g, []).append(s)
    if "multisplit" in groups:
        raise ValueError(
            "group name 'multisplit' is reserved (used by cobinning_multisplit "
            "output paths); rename it in the metadata TSV"
        )
    return groups
