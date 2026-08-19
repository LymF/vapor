"""Regressoes das duas juncoes de predicao de hospedeiro (2026-08-19)."""
from report.data_loaders import _strip_fasta_ext, build_host_collapse


def test_strip_ext_cobre_fna_dos_bins_do_vamb():
    """O track co-assembly usa bin_ext=".fna"; sem isto o Host virava "1.fna"
    e nunca casava com o user_genome "1" do GTDB-Tk."""
    assert _strip_fasta_ext("1.fna") == "1"
    assert _strip_fasta_ext("136.fna") == "136"


def test_strip_ext_nao_mutila_fasta():
    """`.replace('.fa','')` antes de `.fasta` transformava "x.fasta" em "xsta"."""
    assert _strip_fasta_ext("contig_k141_9.fasta") == "contig_k141_9"
    assert _strip_fasta_ext("binette_bin4.fa") == "binette_bin4"


def test_strip_ext_preserva_nome_sem_extensao():
    assert _strip_fasta_ext("binette_bin4") == "binette_bin4"
    # ponto interno nao e extensao
    assert _strip_fasta_ext("P05_AMD_08_956.6") == "P05_AMD_08_956.6"


def _phist(sample, virus, host):
    return {"sample": sample, "Virus": virus, "Host": host}


def test_genero_vem_do_gtdb_nao_do_nome_do_arquivo():
    """Todo bin do Binette comeca com "binette_", entao o split('_')[0] antigo
    colapsava TODOS os hospedeiros num unico "genero" chamado "binette"."""
    phist = [_phist("S1", "k141_1", "binette_bin4"),
             _phist("S1", "k141_2", "binette_bin12")]
    links = [{"sample": "S1", "Host": "binette_bin4",  "Host_genus": "Neisseria"},
             {"sample": "S1", "Host": "binette_bin12", "Host_genus": "Prevotella"}]
    abund = {"S1": [{"representative": "k141_1", "rpkm": "10"},
                    {"representative": "k141_2", "rpkm": "5"}]}

    out = build_host_collapse(phist, abund, ["S1"], host_links=links)
    generos = {e["genus"]: e for e in out["S1"]}
    assert set(generos) == {"Neisseria", "Prevotella"}
    assert "binette" not in generos
    assert generos["Neisseria"]["total_rpkm"] == 10.0


def test_hospedeiro_sem_taxonomia_vira_unknown_nao_nome_de_arquivo():
    phist = [_phist("S1", "k141_1", "1")]
    out = build_host_collapse(phist, {"S1": []}, ["S1"], host_links=[])
    assert [e["genus"] for e in out["S1"]] == ["Unknown"]


def test_enriquecimento_do_checkv_resolve_provirus():
    """Desde que o aparo do CheckV voltou a funcionar, o Genome da taxonomia
    pode ser "k141_9_1" enquanto o quality_summary so conhece "k141_9"."""
    from report.data_loaders import enrich_taxonomy_with_checkv

    checkv = {"S1": [{"contig_id": "k141_9", "completeness": "88.5",
                      "checkv_quality": "High-quality", "contig_length": "31000"}]}
    tax = [{"sample": "S1", "Genome": "k141_9_1"},
           {"sample": "S1", "Genome": "k141_9"},
           {"sample": "S1", "Genome": "k141_desconhecido"}]

    out = enrich_taxonomy_with_checkv(tax, checkv)
    assert out[0]["CheckV_quality"] == "High-quality"   # provirus resolvido
    assert out[1]["CheckV_quality"] == "High-quality"   # contig inteiro
    assert out[2]["CheckV_quality"] == ""               # desconhecido nao inventa
