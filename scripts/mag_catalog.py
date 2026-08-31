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

    Fontes discordam de esquema: uma amostra sem bins escreve so o cabecalho
    de compatibilidade (`Name\tCompleteness\tContamination\tGenome_Size`,
    zero linhas), enquanto uma amostra com bins escreve as 14 colunas reais
    do CheckM2. Ler por NOME de coluna (DictReader/DictWriter), nunca por
    posicao, e essencial: a versao anterior travava o cabecalho da PRIMEIRA
    fonte lida e escrevia as linhas de fontes seguintes por baixo dele sem
    verificar a contagem de colunas -- se essa primeira fonte fosse a de
    esquema curto, linhas de 14 campos iam parar sob um cabecalho de 4, e o
    parser Rust do galah (`checkm-0.3.0`) sofria panic com
    `Option::unwrap() on a None value` ao indexar uma coluna que o cabecalho
    dizia nao existir.
    """
    fieldnames = []
    seen = set()
    rows = []
    n_sources = 0
    for source_id, path in pairs:
        if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
            continue
        with open(path, newline="") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            if not r.fieldnames or CHECKM2_NAME_COL not in r.fieldnames:
                continue
            for col in r.fieldnames:
                if col not in seen:
                    seen.add(col)
                    fieldnames.append(col)
            n_sources += 1
            for row in r:
                name = (row.get(CHECKM2_NAME_COL) or "").strip()
                if not name:
                    continue
                row[CHECKM2_NAME_COL] = namespaced_id(source_id, strip_ext(name))
                rows.append(row)

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", newline="") as out:
        w = csv.DictWriter(
            out, delimiter="\t",
            fieldnames=fieldnames or [CHECKM2_NAME_COL, "Completeness", "Contamination"],
            restval="", extrasaction="ignore",
        )
        w.writeheader()
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


# ── Vistas por fonte ──────────────────────────────────────────────────────
#
# Uma vista implementa HERANCA, nao identidade: todo MAG da fonte recebe a
# linha do SEU representante, e nao apenas os MAGs que por acaso SAO
# representantes. Com 32 amostras, o representante da maioria dos clusters
# pertence a outra amostra -- filtrar a tabela global pelo prefixo da fonte
# devolveria quase nada, em silencio. Foi exatamente esse o bug do
# `viral_taxonomy` em 2026-08-18.


def member_map(membership_rows, source_id):
    """{representative_id: [original_bin_id, ...]} para UMA fonte.

    Lista e nao escalar: duas linhagens da mesma especie no mesmo metagenoma
    dao dois bins com o mesmo representante.
    """
    out = {}
    for row in membership_rows:
        if (row.get("source_id") or "").strip() != source_id:
            continue
        rep = (row.get("representative_id") or "").strip()
        bin_name = (row.get("original_bin_id") or "").strip()
        if rep and bin_name:
            out.setdefault(rep, []).append(bin_name)
    return out


def resolve_prefixed_id(value, representatives):
    """('{rep}__{resto}', {reps}) -> (rep, resto). (None, value) se nao casar.

    NAO se pode cortar no primeiro separador. O ID do catalogo ja e
    '{source}__{bin}', entao uma proteina sai como
    'S1__binette_bin1__k141_1_5' e o corte no primeiro '__' devolveria 'S1'
    -- todo hit de AMR atribuido a AMOSTRA em vez do MAG. Aqui casa-se
    contra os representantes conhecidos, do prefixo mais longo para o mais
    curto, que e a mesma disciplina de `scripts/checkv_provirus.py`: resolver
    contra IDs conhecidos em vez de adivinhar onde fica o delimitador.
    """
    value = (value or "").strip()
    if not value:
        return None, value
    parts = value.split(SEP)
    for n in range(len(parts) - 1, 0, -1):
        cand = SEP.join(parts[:n])
        if cand in representatives:
            return cand, SEP.join(parts[n:])
    return None, value
