"""O BACPHLIP sobre genomas detectados pelo sylph (2026-08-19)."""
import pandas as pd
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "scripts", "reads_classify"))
from bacphlip_lifestyle import _detected_ids


def test_imgvr_nao_colapsa_num_unico_id():
    """Todo virus do IMG/VR compartilha Genome_file "imgvr_reps.fna" -- so
    Contig_name distingue. Medido nos dados da Amazonia."""
    df = pd.DataFrame({
        "Genome_file": ["imgvr_reps.fna", "imgvr_reps.fna"],
        "Contig_name": ["IMGVR_UViG_3300047463_000035|3300047463|Ga0466545_004007",
                        "IMGVR_UViG_3300047401_000254|3300047401|Ga0466547_000631"],
    })
    ids = _detected_ids(df)
    assert "IMGVR_UViG_3300047463_000035" in ids
    assert "IMGVR_UViG_3300047401_000254" in ids
    # o comportamento antigo produzia exatamente isto e mais nada
    assert ids != {"imgvr_reps"}


def test_gtdb_continua_funcionando():
    df = pd.DataFrame({
        "Genome_file": ["gtdb/GCA/044/GCA_044225315.1_genomic.fna.gz"],
        "Contig_name": ["contig_1"],
    })
    ids = _detected_ids(df)
    assert "GCA_044225315" in ids


def test_tolera_ausencia_de_contig_name():
    df = pd.DataFrame({"Genome_file": ["GCF_000334915.1_genomic.fna.gz"]})
    assert "GCF_000334915" in _detected_ids(df)
