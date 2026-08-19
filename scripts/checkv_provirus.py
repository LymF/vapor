"""
Recuperar o contig ORIGINAL a partir de um header de provirus do CheckV.

Por que isto e um modulo, e nao tres linhas inline: as tres linhas inline sao
exatamente o que produziu o bug. `viral_binning.smk`, `coassembly.smk` e o
`lookup_key` do gate de grupo cada um trazia sua propria copia de
`hdr_id.rsplit('|', 1)[0]`, apoiada num comentario que afirmava que o
`proviruses.fna` usa headers "orig_id|start_end". Ele NAO usa.

Verificado contra os dados da Amazonia em 2026-08-19 (7 grupos de co-assembly,
118 proviroses):

    >k141_219139_1 1-13933/18998
    >k141_97527_1 1-2446/3389

O formato real e "{contig_id}_{n}", n = indice da regiao de provirus dentro do
contig, e o intervalo fica na DESCRICAO, depois do espaco -- que todo leitor
descarta com .split()[0]. Como nao ha "|" nenhum, o rsplit devolvia o header
inteiro, a chave nunca batia com o contig original, e o resultado era duplo:
o fragmento aparado era descartado E a sequencia original, com o DNA de
hospedeiro flanqueando o profago, era emitida como "viral".

A resolucao aqui NAO adivinha delimitador. Ela testa candidatos contra o
conjunto de IDs que se sabe existirem (`known`), o que e a unica forma de isto
falhar alto em vez de em silencio se o CheckV mudar o formato de novo.
"""


def resolve_original_id(header_id, known):
    """
    Devolve (orig_id, resolvido) para um header ja reduzido a .split()[0].

    `known` e o conjunto de contig ids originais (chaves do consenso ou do
    quality_summary). `resolvido` e False quando nenhum candidato bate --
    o chamador deve CONTAR esses casos e registra-los, porque um formato
    novo do CheckV aparece exatamente assim.
    """
    if header_id in known:
        return header_id, True

    # "{contig}_{n}" -- formato real do CheckV (>= 1.0).
    base, sep, suffix = header_id.rpartition("_")
    if sep and suffix.isdigit() and base in known:
        return base, True

    # "{contig}|{start}_{end}" -- formato que o comentario antigo afirmava.
    # Mantido porque custa nada e desambigua contra `known`, nao por fe.
    base = header_id.rsplit("|", 1)[0]
    if base != header_id and base in known:
        return base, True

    return header_id, False


def build_trimmed_index(fasta_entries, known):
    """
    fasta_entries: iteravel de (header_line, seq_lines, is_provirus).
    Devolve (index, n_nao_resolvidos), index = {orig_id: [(header, seq), ...]}.
    """
    from collections import defaultdict

    index = defaultdict(list)
    unresolved = 0
    for header_line, seq_lines, is_provirus in fasta_entries:
        hdr_id = header_line[1:].split()[0]
        if is_provirus:
            orig, ok = resolve_original_id(hdr_id, known)
            if not ok:
                unresolved += 1
        else:
            orig = hdr_id
        index[orig].append((header_line, seq_lines))
    return index, unresolved


def inherit_from_original(ids, known):
    """Resolve varios ids de uma vez contra `known`.

    Devolve ({id: orig_id}, {"direct": n, "trimmed": n, "unresolved": n}).
    `direct` conta os ids que ja existem em `known` (sem aparo), `trimmed` os
    que so casaram depois de desfazer o sufixo "{contig}_{n}", e `unresolved`
    os que nao casaram com nada -- estes ultimos TEM de ser contados pelo
    chamador e registrados no log: um formato novo do CheckV aparece
    exatamente como um monte de nao resolvidos, nao como um erro.

    Existe para os chamadores nao repetirem o laco: as tres copias inline de
    `hdr_id.rsplit('|', 1)[0]` foram o que produziu o bug que este modulo
    documenta no topo.
    """
    mapping = {}
    stats = {"direct": 0, "trimmed": 0, "unresolved": 0}
    for i in ids:
        if i in known:
            mapping[i] = i
            stats["direct"] += 1
            continue
        orig, ok = resolve_original_id(i, known)
        if ok:
            mapping[i] = orig
            stats["trimmed"] += 1
        else:
            stats["unresolved"] += 1
    return mapping, stats
