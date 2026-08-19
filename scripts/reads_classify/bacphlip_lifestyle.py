"""
Extract detected viral genome sequences and run BACPHLIP for lifestyle prediction.

Usage:
    python bacphlip_lifestyle.py <sylph_results.tsv> <genomes.fasta> <out_dir> <virulence_threshold>

sylph_results.tsv : saida crua do sylph. Os IDs saem de Contig_name E de
                    Genome_file -- ver _detected_ids() para por que os dois.
genomes.fasta     : concatenated FASTA of all reference genomes in the sylph DB
out_dir           : directory to write extracted FASTA + BACPHLIP results
virulence_threshold: float 0–1; scores above this = virulent (default 0.5)
"""
import os
import sys
import subprocess
import tempfile
import pandas as pd


def _detected_ids(df):
    """IDs dos genomas que o sylph detectou.

    Usava so `basename(Genome_file).split(".")[0]`, o que funciona para o
    GTDB -- um arquivo por genoma
    ("GCA_044225315.1_genomic.fna.gz" -> "GCA_044225315") -- mas NAO para o
    IMG/VR, onde todo virus vem do MESMO arquivo. Confirmado nos dados da
    Amazonia em 2026-08-19: 111 linhas de sylph, e a coluna Genome_file diz
    "imgvr_reps.fna" em todas as virais. O conjunto detectado colapsava para
    um unico id, "imgvr_reps", que nao e header de nada -- ou seja, com um
    banco IMG/VR o BACPHLIP extraia zero sequencias e o script abortava.
    Quem distingue os virus ali e Contig_name:
    "IMGVR_UViG_3300047463_000035|3300047463|Ga0466545_004007".

    Emite candidatos das DUAS colunas em vez de escolher: o header da FASTA
    de referencia e que decide qual casa, e nao ha como saber qual banco o
    usuario apontou sem le-lo.
    """
    ids = set()
    if "Contig_name" in df.columns:
        for name in df["Contig_name"].dropna():
            name = str(name).strip()
            if not name:
                continue
            ids.add(name)                       # header completo
            ids.add(name.split("|")[0])         # UViG do IMG/VR
    if "Genome_file" in df.columns:
        for path in df["Genome_file"].dropna():
            base = os.path.basename(str(path))
            ids.add(base)
            ids.add(base.split(".")[0])
    return {i for i in ids if i}


def main():
    # Biopython so existe no env da regra (env_reads_classify). Importar no
    # topo tornava _detected_ids inimportavel fora dele, e portanto nao
    # testavel -- que e como o bug do IMG/VR sobreviveu.
    from Bio import SeqIO

    if len(sys.argv) < 4:
        sys.exit("Usage: bacphlip_lifestyle.py sylph_results.tsv genomes.fasta out_dir [virulence_threshold]")

    tsv_path, fasta_path, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    vir_thresh = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5

    os.makedirs(out_dir, exist_ok=True)

    # Parse detected genome IDs from sylph output
    df = pd.read_csv(tsv_path, sep="\t")
    if "Genome_file" not in df.columns and "Contig_name" not in df.columns:
        sys.exit("[bacphlip_lifestyle] nem Genome_file nem Contig_name no TSV do sylph")

    detected_ids = _detected_ids(df)
    print(f"[bacphlip_lifestyle] {len(detected_ids)} genomes detected by sylph")

    # Extract detected sequences from the reference FASTA
    extracted_fasta = os.path.join(out_dir, "detected_genomes.fasta")
    found = 0
    with open(extracted_fasta, "w") as fh:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            rec_id = rec.id.split(".")[0]
            if rec_id in detected_ids or rec.id in detected_ids:
                SeqIO.write(rec, fh, "fasta")
                found += 1

    print(f"[bacphlip_lifestyle] extracted {found}/{len(detected_ids)} sequences")
    if found == 0:
        sys.exit(
            "[bacphlip_lifestyle] nenhuma sequencia extraida. Os %d ids "
            "detectados nao casaram com nenhum header de %s. Exemplos de id: "
            "%s" % (len(detected_ids), fasta_path,
                    ", ".join(sorted(detected_ids)[:5]) or "(nenhum)")
        )

    # Run BACPHLIP
    bacphlip_out = os.path.join(out_dir, "bacphlip_results.tsv")
    subprocess.run(
        ["bacphlip", "-i", extracted_fasta, "--multi_fasta"],
        check=True
    )
    # BACPHLIP writes to <input>.bacphlip by default; rename
    default_out = extracted_fasta + ".bacphlip"
    if os.path.exists(default_out):
        os.rename(default_out, bacphlip_out)

    # Compute lifestyle summary
    # BACPHLIP writes 3 columns: sequence name (as the DataFrame index) plus
    # both probabilities -- its header line is literally "\tVirulent\tTemperate"
    # (bacphlip.py:271). Naming only 2 columns here raised a Length mismatch.
    bp = pd.read_csv(bacphlip_out, sep="\t")
    bp.columns = ["sequence_name", "virulent_score", "temperate_score"]
    bp["lifestyle"] = bp["virulent_score"].apply(
        lambda s: "Virulent" if s > vir_thresh else "Temperate"
    )

    summary_path = os.path.join(out_dir, "lifestyle_summary.tsv")
    bp.to_csv(summary_path, sep="\t", index=False)

    n_virulent  = (bp["lifestyle"] == "Virulent").sum()
    n_temperate = (bp["lifestyle"] == "Temperate").sum()
    ratio = round(n_virulent / len(bp), 3) if len(bp) > 0 else 0.0
    print(
        f"[bacphlip_lifestyle] Virulent: {n_virulent}  Temperate: {n_temperate}  "
        f"Ratio: {ratio:.1%}"
    )
    print(f"[bacphlip_lifestyle] results written to {summary_path}")


if __name__ == "__main__":
    main()
