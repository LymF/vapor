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

'?' cobre duas causas distintas, e quem chama este modulo precisa saber
separa-las no log (nao aqui, que so decide o estado):
  1. completude abaixo do piso (`min_completeness`);
  2. falha de anotacao por membro -- prodigal ou DefenseFinder falhou
     naquele genoma, ou o genoma nem chegou ao pool. Um membro assim nao
     aparece no manifesto do pangenoma, e o chamador passa isso aqui via
     `annotated`: quem nao esta em `annotated` e sempre '?', mesmo que a
     completude esteja ok -- falha de ferramenta nao e ausencia biologica.

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
                 annotated=None, min_completeness=MIN_COMPLETENESS):
    """Uma linha por (cluster, tipo, gene).

    `annotated`: conjunto opcional de membros que de fato aparecem no
    manifesto de anotacao (prodigal + defensefinder/amr rodaram e nao
    falharam). Quando fornecido, um membro fora dele e sempre '?' e sai
    do denominador -- ver docstring do modulo.
    """
    rows = []
    for rep in clusters:
        members = members_by_rep.get(rep, [])
        evaluable = set(_evaluable(members, completeness, min_completeness,
                                   annotated))
        typed_hits = {m: {_typed(h) for h in gene_hits.get(m, set())}
                     for m in members}
        keys = sorted({k for m in members for k in typed_hits[m]})
        for tipo, gene in keys:
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
    """Uma linha por cluster: e isto que decide se a fase 2 se justifica."""
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
            if not n_eval or row["n_present"] == 0:
                continue
            if row["n_present"] >= CORE_FRACTION * n_eval:
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
