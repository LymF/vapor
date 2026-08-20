#!/usr/bin/env python3
"""Matriz gene x membro dos clusters elegiveis, em TRES estados.

Ausencia em MAG nao e ausencia no organismo. Um "." pode significar "o
organismo nao tem o gene" ou "o MAG esta 74% completo e a regiao nao
montou"; tratar os dois como o mesmo zero produz conclusao errada com
aparencia de dado -- a mesma familia do done.txt vazio lido como zero
biologico (docs/ROADMAP_SIMPLIFICACAO.md).

Por isso: 'x' presente, '.' ausente, '?' NAO AVALIAVEL. E o '?' fica fora
do denominador da frequencia.

Existe um quarto simbolo, fora do trio acima: '-' e escrito pelo chamador
(rules/pangenome.smk) para um membro que NAO pertence ao cluster da linha
(a matriz e larga, colunas = uniao de todos os membros elegiveis). Este
modulo nao produz '-' -- so 'x'/'./'?'.

'?' cobre TRES causas distintas, e quem chama este modulo precisa saber
separa-las no log (nao aqui, que so decide o estado):
  1. completude abaixo do piso (`min_completeness`);
  2. falha de anotacao por membro -- prodigal falhou naquele genoma, ou o
     genoma nem chegou ao pool. Um membro assim nao aparece no manifesto
     do pangenoma, e o chamador passa isso aqui via `annotated`: quem nao
     esta em `annotated` e sempre '?' em TODAS as linhas (defesa e amr),
     mesmo que a completude esteja ok -- falha de ferramenta nao e
     ausencia biologica;
  3. falha do DefenseFinder por membro (o `|| echo WARNING` em
     `rules/defense_amr.smk`, que deliberadamente nao derruba a regra):
     prodigal funcionou (o membro esta em `annotated`), mas o
     DefenseFinder nao produziu saida para aquele genoma especifico. Isto
     e passado via `defense_failed`, e so afeta as linhas `tipo='defesa'`
     desse membro -- AMR vem de AMRFinderPlus/RGI/DeepARG, que rodam
     sobre proteinas concatenadas (nao em laco por genoma), entao uma
     falha ali e GLOBAL e ja fica coberta pelo `failed:` da propria
     regra, sem equivalente de `defense_failed` necessario aqui.

Cada linha tambem carrega um `tipo` ('defesa' ou 'amr'): sistema de defesa
e ARG entram na mesma tabela de genes, e sem essa distincao uma colisao de
NOME entre os dois fundiria duas linhas em uma so. A chave interna de uma
linha e o par (tipo, nome), nunca so o nome -- ver `_typed()`.

CORE, com menos de 10 membros avaliaveis: `CORE_FRACTION = 0.90` so
discrimina de verdade a partir de ~10 avaliaveis (9/10 e core sob 0.90 e
variavel sob 0.99). Abaixo disso, `0.90 * n_eval` arredonda para
`n_eval` (3 -> 2.7 -> 3; 6 -> 5.4 -> 6): "core" e OPERACIONALMENTE
"presente em TODOS os avaliaveis", sem nenhuma tolerancia de fato. Isto e
esperado no regime de clusters pequenos que a fase 1 produz -- nao ler
`n_genes_core` como se houvesse folga ate ter clusters de ~10+ membros.
"""

MIN_COMPLETENESS = 70.0   # mesmo piso do mag_bakta
CORE_FRACTION = 0.90      # 99% zera o core com MAG (metaFun)


def load_defensefinder_summary(path):
    """Le defensefinder_summary.tsv (bin TAB status) e devolve os FALHOS.

    So 'ok' conta como sucesso -- qualquer outro valor de status
    (`'failed'`, ou um valor desconhecido/futuro que este parser nunca viu)
    entra no conjunto retornado, porque o efeito pretendido (tirar a linha
    `tipo='defesa'` do denominador) e o mesmo de uma falha: o unico jeito
    seguro de ler um status que nao se reconhece e como "nao confirmado
    como ok", nunca como sucesso silencioso. Um arquivo so com cabecalho
    (caminho de skip da regra) devolve conjunto vazio.
    """
    import csv

    failed = set()
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genome = (row.get("genome") or "").strip()
            status = (row.get("status") or "").strip()
            if genome and status != "ok":
                failed.add(genome)
    return failed


def _evaluable(members, completeness, min_completeness, annotated=None):
    return [m for m in members
            if (annotated is None or m in annotated)
            and completeness.get(m, 0.0) >= min_completeness]


def _typed(hit):
    """Normaliza um item de gene_hits para (tipo, nome).

    O chamador pode passar um par (tipo, nome) -- 'defesa'/'amr' -- para
    que uma colisao de nome entre um sistema de defesa e um ARG nao funda
    duas linhas em silencio (a chave interna passa a ser o par, nao so o
    nome). Uma string solta (uso antigo/testes) vira ("", nome).
    """
    if isinstance(hit, tuple) and len(hit) == 2:
        return hit
    return ("", hit)


def build_matrix(clusters, members_by_rep, gene_hits, completeness,
                 annotated=None, min_completeness=MIN_COMPLETENESS,
                 defense_failed=None):
    """Uma linha por (cluster, tipo, gene).

    `annotated`: conjunto opcional de membros que de fato aparecem no
    manifesto de anotacao (prodigal rodou e nao falhou). Quando
    fornecido, um membro fora dele e sempre '?' em TODAS as linhas (sai
    do denominador de defesa e do de amr) -- ver docstring do modulo.

    `defense_failed`: conjunto opcional de membros para os quais o
    DefenseFinder falhou (mas que ESTAO em `annotated`, ou seja, prodigal
    rodou normalmente). Torna '?' e tira do denominador APENAS as linhas
    `tipo='defesa'` desse membro -- as linhas `tipo='amr'` do mesmo
    membro nao sao afetadas, porque AMR vem de outras ferramentas.
    """
    defense_failed = defense_failed or set()
    rows = []
    for rep in clusters:
        members = members_by_rep.get(rep, [])
        base_evaluable = set(_evaluable(members, completeness, min_completeness,
                                        annotated))
        typed_hits = {m: {_typed(h) for h in gene_hits.get(m, set())}
                     for m in members}
        keys = sorted({k for m in members for k in typed_hits[m]})
        for tipo, gene in keys:
            evaluable = {m for m in base_evaluable
                        if not (tipo == "defesa" and m in defense_failed)}
            states, n_present = {}, 0
            for m in members:
                if m not in evaluable:
                    states[m] = "?"
                elif (tipo, gene) in typed_hits[m]:
                    states[m] = "x"
                    n_present += 1
                else:
                    states[m] = "."
            rows.append({
                "representative_id": rep,
                "tipo": tipo,
                "gene": gene,
                "states": states,
                "n_present": n_present,
                "n_evaluable": len(evaluable),
                "freq": f"{n_present}/{len(evaluable)}" if evaluable else "0/0",
            })
    return rows


def summarize_clusters(matrix_rows, members_by_rep, completeness,
                       annotated=None, min_completeness=MIN_COMPLETENESS):
    """Uma linha por cluster: e isto que decide se a fase 2 se justifica.

    O denominador de core/variavel usa `row["n_evaluable"]`, o MESMO valor
    que `build_matrix` ja gravou em cada linha -- nao um `n_eval` calculado
    de novo aqui a partir de `annotated`/`completeness`. Isto e proposital:
    desde que `defense_failed` passou a tirar membros do denominador SO das
    linhas `tipo='defesa'` (ver `build_matrix`), um `n_eval` unico por
    cluster estaria certo para as linhas `amr` e errado para as `defesa` --
    reusar `row["n_evaluable"]` elimina essa segunda fonte de verdade e
    garante que este sumario NUNCA discorda da matriz sobre quantos membros
    contam para cada gene.

    `n_members_avaliaveis`, ao contrario, e um numero por CLUSTER (nao por
    linha): quantos membros passaram no piso geral de completude + prodigal.
    E o denominador "largo" que as linhas `amr` usam; as linhas `defesa` de
    um membro com DefenseFinder falho sao mais estreitas que isto, e essa
    diferenca so aparece na propria matriz (coluna `n_evaluable`), nao aqui.
    """
    by_rep = {}
    for row in matrix_rows:
        by_rep.setdefault(row["representative_id"], []).append(row)

    summaries = []
    for rep, rows in sorted(by_rep.items()):
        members = members_by_rep.get(rep, [])
        n_eval = len(_evaluable(members, completeness, min_completeness,
                                annotated))
        core = variable = singleton = 0
        for row in rows:
            n_eval_row = row["n_evaluable"]
            if not n_eval_row or row["n_present"] == 0:
                continue
            if row["n_present"] >= CORE_FRACTION * n_eval_row:
                core += 1
            else:
                variable += 1
            if row["n_present"] == 1:
                singleton += 1
        summaries.append({
            "representative_id":    rep,
            "n_members":            len(members),
            "n_members_avaliaveis": n_eval,
            "n_genes_core":         core,
            "n_genes_variaveis":    variable,
            "n_genes_singleton":    singleton,
        })
    return summaries
