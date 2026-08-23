"""Inventário da matriz de status.

Até 2026-08-23 só três regras per-sample eram rastreadas — e as três eram
justamente VISTAS sobre o catálogo global (`bins/gtdbtk`, `bins/amrfinderplus`,
`bins/rgi` são escritas por `mag_views_sample`, não pela ferramenta). O
resultado era duplamente ruim: as regras que realmente rodam por amostra —
fastp, montagem, mapeamento, binning, detecção viral — não apareciam em lugar
nenhum quando falhavam, e as três que apareciam sugeriam que a ferramenta
tinha rodado naquela amostra, quando ela rodou uma vez sobre o representante
do cluster.

Este teste fixa as duas correções: o inventário cobre as regras per-sample de
verdade, e toda entrada declara se é regra ou vista.
"""
import os

import pytest

from report.data_loaders import (
    STATUS_TRACKED_TOOLS, STATUS_VIEW_TOOLS, STATUS_TRACKED_GLOBAL_TOOLS,
    load_tool_status, GLOBAL_STATUS_LABEL,
)


def _escreve(caminho, texto):
    os.makedirs(os.path.dirname(caminho), exist_ok=True)
    with open(caminho, 'w', encoding='utf-8') as fh:
        fh.write(texto)


def test_o_inventario_cobre_as_etapas_per_sample_de_verdade():
    # O MEGAHIT nao entra: nao escreve done.txt, e uma falha dele derruba o
    # DAG inteiro em vez de emitir tabela vazia -- e a diferenca entre regra
    # soft-fail (que precisa de sentinela) e regra que quebra alto.
    esperadas = {"fastp", "mapping", "metabat2", "semibin2",
                 "binette", "vrhyme", "genomad", "final"}
    faltando = esperadas - set(STATUS_TRACKED_TOOLS)
    assert not faltando, (
        "estas regras rodam por amostra e uma falha nelas nao apareceria em "
        f"lugar nenhum do report: {sorted(faltando)}")


def test_as_vistas_nao_se_disfarcam_de_regra_per_sample():
    # bins/gtdbtk e bins/amrfinderplus sao ESCRITAS pelo mag_views_sample. Um
    # 'ok' ali nao diz que a ferramenta rodou nesta amostra -- diz que a
    # distribuicao do resultado global chegou.
    assert "gtdbtk" in STATUS_VIEW_TOOLS
    assert "amrfinderplus" in STATUS_VIEW_TOOLS
    assert not (set(STATUS_VIEW_TOOLS) & set(STATUS_TRACKED_TOOLS))


def test_cada_entrada_declara_se_e_regra_ou_vista(tmp_path):
    o = str(tmp_path)
    _escreve(os.path.join(o, "S1", STATUS_TRACKED_TOOLS["fastp"]), "ok\n")
    _escreve(os.path.join(o, "S1", STATUS_VIEW_TOOLS["gtdbtk"]), "ok\n")

    status = load_tool_status(o, ["S1"])
    assert status["S1"]["fastp"]["kind"] == "rule"
    assert status["S1"]["gtdbtk"]["kind"] == "view"
    assert status[GLOBAL_STATUS_LABEL]["mag_catalog_derep"]["kind"] == "global"


def test_falha_continua_distinta_de_ausencia(tmp_path):
    o = str(tmp_path)
    _escreve(os.path.join(o, "S1", STATUS_TRACKED_TOOLS["genomad"]),
             "failed: disco cheio\n")
    status = load_tool_status(o, ["S1"])
    assert status["S1"]["genomad"]["state"] == "failed"
    assert status["S1"]["genomad"]["reason"] == "disco cheio"
    # Regra que nunca escreveu done.txt e 'unknown', nunca 'ok'.
    assert status["S1"]["binette"]["state"] == "unknown"


def test_nenhum_caminho_do_inventario_e_absoluto_ou_escapa_do_diretorio():
    # Um caminho com '..' ou com barra inicial faria o status de uma amostra
    # ser lido de fora do diretorio dela.
    for nome, rel in {**STATUS_TRACKED_TOOLS, **STATUS_VIEW_TOOLS,
                      **STATUS_TRACKED_GLOBAL_TOOLS}.items():
        assert not os.path.isabs(rel), nome
        assert ".." not in rel.split("/"), nome
