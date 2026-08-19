"""
Catalogo global de MAGs procarioticos — helpers puros.

Mesmo principio do catalogo de vOTU (`docs/ROADMAP_SIMPLIFICACAO.md`, item
"(h)"): computar no representante, herdar no membro. 95% ANI e nivel de
especie para procarioto tambem, entao um MAG da mesma especie recuperado em
tres amostras nao precisa ser anotado tres vezes -- e, mais importante, nao
DEVE, porque nada garantia que as tres recebessem a mesma anotacao.

O NAMESPACE nao e cosmetico. Nomes de bin colidem entre amostras: o Binette
emite "binette_bin1" em toda amostra e o VAMB emite numeros nus ("1", "136").
Juntar os arquivos num diretorio so sem prefixar sobrescreveria bins de
organismos diferentes em silencio -- a mesma familia de bug que este roadmap
inteiro persegue. Aqui o prefixo e o proprio nome do arquivo, porque o galah
recebe caminhos, nao sequencias.
"""
import csv
import os


SEP = "__"


class NamespaceCollision(RuntimeError):
    """Um source_id que contem o separador tornaria o ID ambiguo."""


def namespaced_id(source_id, bin_name):
    """'{source}__{bin}'. Levanta se o source_id contiver o separador.

    Falhar alto aqui e deliberado: um source_id com '__' faria
    `split_namespaced_id` devolver o pedaco errado, e o resultado seria uma
    anotacao atribuida a amostra errada -- silenciosamente, que e o pior modo
    de falha possivel neste codigo.
    """
    if SEP in str(source_id):
        raise NamespaceCollision(
            "source_id %r contem %r, que e o separador de namespace do "
            "catalogo de MAGs. Renomeie a amostra/grupo." % (source_id, SEP)
        )
    return "%s%s%s" % (source_id, SEP, bin_name)


def split_namespaced_id(member_id):
    """Inverso de namespaced_id: (source_id, bin_name).

    Usa split no PRIMEIRO separador, nao no ultimo: nomes de bin podem conter
    '__' (nada os proibe), source_ids nao podem (namespaced_id recusa).
    """
    source, sep, rest = str(member_id).partition(SEP)
    return (source, rest) if sep else ("", str(member_id))


def strip_ext(filename, exts=(".fa", ".fna", ".fasta")):
    """Nome do bin sem a extensao, a mais longa primeiro."""
    base = os.path.basename(filename)
    for e in sorted(exts, key=len, reverse=True):
        if base.endswith(e):
            return base[: -len(e)]
    return base


def collect_bins(sources):
    """[(source_type, source_id, bin_name, caminho)] para cada bin existente.

    sources: [(source_type, source_id, bins_dir, bin_ext)]. Diretorios
    ausentes ou vazios sao pulados sem erro -- uma amostra sem MAG e um
    resultado real, nao uma falha.
    """
    import glob

    found = []
    for source_type, source_id, bins_dir, bin_ext in sources:
        if not bins_dir or not os.path.isdir(bins_dir):
            continue
        for path in sorted(glob.glob(os.path.join(bins_dir, "*" + bin_ext))):
            found.append((source_type, source_id, strip_ext(path), path))
    return found


def build_pool(sources, out_dir, provenance_path):
    """Symlinks com nome namespaced + provenance.tsv.

    Symlink, nao copia: MAGs sao arquivos grandes e o galah so precisa
    le-los. O destino e sempre um caminho absoluto, senao o link quebra
    quando lido de outro diretorio.
    """
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(provenance_path) or ".", exist_ok=True)

    entries = collect_bins(sources)
    n_sources = len({(t, s) for t, s, _, _ in entries})

    with open(provenance_path, "w", newline="") as prov:
        w = csv.writer(prov, delimiter="\t")
        w.writerow(["member_id", "source_type", "source_id", "original_bin_id"])
        for source_type, source_id, bin_name, path in entries:
            member = namespaced_id(source_id, bin_name)
            link = os.path.join(out_dir, member + ".fa")
            if os.path.islink(link) or os.path.exists(link):
                os.remove(link)
            os.symlink(os.path.abspath(path), link)
            w.writerow([member, source_type, source_id, bin_name])

    return {"n_genomes": len(entries), "n_sources": n_sources}


CHECKM2_NAME_COL = "Name"


def merge_checkm2(pairs, out_path):
    """Um quality_report.tsv do CheckM2 para o catalogo inteiro.

    pairs: [(source_id, caminho_do_quality_report)].

    O galah recebe `--checkm2-quality-report` e casa pelo nome do genoma. Como
    o pool renomeia todo bin para o ID namespaced, a coluna Name TEM de ser
    reescrita para o mesmo ID -- concatenar os relatorios crus deixaria o
    galah sem qualidade para nenhum genoma, e ele cairia para escolher o
    representante por outro criterio sem avisar.

    Devolve (n_linhas, n_fontes_lidas).
    """
    header = None
    rows = []
    n_sources = 0
    for source_id, path in pairs:
        if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
            continue
        with open(path, newline="") as fh:
            r = csv.reader(fh, delimiter="\t")
            h = next(r, None)
            if h is None:
                continue
            if CHECKM2_NAME_COL not in h:
                continue
            idx = h.index(CHECKM2_NAME_COL)
            if header is None:
                header = h
            n_sources += 1
            for row in r:
                if idx >= len(row) or not row[idx].strip():
                    continue
                row = list(row)
                row[idx] = namespaced_id(source_id, strip_ext(row[idx].strip()))
                rows.append(row)

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(header or [CHECKM2_NAME_COL, "Completeness", "Contamination"])
        w.writerows(rows)
    return len(rows), n_sources


def parse_galah_clusters(path):
    """galah --output-cluster-definition -> {member_id: representative_id}.

    O galah escreve CAMINHOS de arquivo, nao ids. Converte para o ID
    namespaced (basename sem extensao), que e a chave usada em todo o resto
    do catalogo.
    """
    mapping = {}
    if not path or not os.path.exists(path):
        return mapping
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 2:
                continue
            if row[0].strip().lower() == "representative":
                continue          # cabecalho do fallback
            rep, member = strip_ext(row[0].strip()), strip_ext(row[1].strip())
            if rep and member:
                mapping[member] = rep
    return mapping


def representative_view(clusters, provenance):
    """Linhas (source_id, bin_name, member_id, representative_id).

    E isto que permite "as analises so nas representantes, mas sabendo de qual
    amostra cada MAG veio": a anotacao vive no representante e a vista por
    amostra e uma juncao, exatamente como `viral_taxonomy` virou uma vista
    sobre o catalogo de vOTU.
    """
    out = []
    for member, source_type, source_id, bin_name in provenance:
        out.append((source_id, bin_name, member, clusters.get(member, member)))
    return out


def read_provenance(path):
    """[(member_id, source_type, source_id, original_bin_id)]."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rows.append((
                (row.get("member_id") or "").strip(),
                (row.get("source_type") or "").strip(),
                (row.get("source_id") or "").strip(),
                (row.get("original_bin_id") or "").strip(),
            ))
    return rows
