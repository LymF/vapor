"""
Recuperar o hospedeiro anotado pelo BANCO a partir dos .sylphmpa por amostra.

Por que este script existe: o `sylph-tax merge` (rule sylph_merge) roda com
`--column relative_abundance`, ou seja, so a coluna de abundancia sobrevive na
tabela mesclada. A coluna "Virus_host (if viral)" -- que e o unico lugar onde o
hospedeiro anotado pelo banco aparece -- fica para tras. Medido nos dados da
Amazonia em 2026-08-19: os .sylphmpa por amostra tem a coluna (107 linhas com
valor em P01_RNG_08_947), a tabela mesclada nao tem a coluna nenhuma.

Consequencia antes disto: `collapse_by_host.py` nunca via um hospedeiro, caia
sempre no ramo "no 'host' column found" e escrevia a saida em nivel de taxon.
O agrupamento por hospedeiro nunca chegou a rodar.

Saida: TSV com clade_name<TAB>host_db, uma linha por clado viral com
hospedeiro conhecido. Clados cujo host e UNKNOWN/NA sao omitidos -- ausencia
declarada e melhor que uma linha dizendo "UNKNOWN", que depois se confunde com
"o banco nao foi consultado".
"""
import csv
import os
import sys


_HOST_COL_HINT = "virus_host"
_NULL = {"", "na", "nan", "unknown", "none", "-"}


def _host_column(header):
    """Indice da coluna de hospedeiro.

    O sylph-tax escreve o cabecalho como "Virus_host (if viral)" -- com espaco
    e parenteses. Casar por prefixo normalizado em vez de string exata para
    nao quebrar se a anotacao entre parenteses mudar.
    """
    for i, col in enumerate(header):
        if col.strip().lower().replace(" ", "_").startswith(_HOST_COL_HINT):
            return i
    return None


def _is_null_host(value):
    """Um host so conta se algum rank for real.

    O IMG/VR emite "UNKNOWN;UNKNOWN;UNKNOWN;..." quando nao atribuiu
    hospedeiro -- string longa, nao vazia, que passaria por qualquer teste
    ingenuo de "tem valor".
    """
    v = (value or "").strip()
    if v.lower() in _NULL:
        return True
    return all(p.strip().lower() in _NULL for p in v.split(";"))


def collect(paths):
    """{clade_name: host_lineage} para os clados com hospedeiro conhecido."""
    hosts = {}
    n_rows = n_with_host = 0
    for path in paths:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            continue
        with open(path, newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = None
            for row in reader:
                if not row:
                    continue
                # a 1a linha do .sylphmpa e um comentario "#SampleID\t..."
                if row[0].startswith("#"):
                    continue
                if header is None:
                    header = row
                    host_i = _host_column(header)
                    continue
                n_rows += 1
                if host_i is None or host_i >= len(row):
                    continue
                clade, host = row[0].strip(), row[host_i].strip()
                if not clade or _is_null_host(host):
                    continue
                n_with_host += 1
                # mesma amostra ou outra: o host do banco e propriedade do
                # clado, nao da amostra -- o primeiro nao-nulo vale.
                hosts.setdefault(clade, host)
    return hosts, n_rows, n_with_host


def main():
    if len(sys.argv) < 3:
        sys.exit("Uso: build_host_map.py <saida.tsv> <arquivo.sylphmpa> [...]")
    out_path, paths = sys.argv[1], sys.argv[2:]
    hosts, n_rows, n_with_host = collect(paths)

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["clade_name", "host_db"])
        for clade in sorted(hosts):
            w.writerow([clade, hosts[clade]])

    sys.stderr.write(
        "[build_host_map] %d arquivos, %d linhas, %d com host nao-nulo, "
        "%d clados unicos -> %s\n"
        % (len(paths), n_rows, n_with_host, len(hosts), out_path)
    )
    if n_rows and not hosts:
        sys.stderr.write(
            "[build_host_map] AVISO: nenhum clado com hospedeiro. Isso e "
            "esperado se o banco nao anota host (ou anota tudo como UNKNOWN, "
            "como o IMG/VR faz para virus sem hospedeiro atribuido) -- nao "
            "confunda com falha de leitura.\n"
        )


if __name__ == "__main__":
    main()
