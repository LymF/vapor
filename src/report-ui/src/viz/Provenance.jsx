// Marca de procedência.
//
// Sob o princípio (h) do ROADMAP_SIMPLIFICACAO.md, a maior parte do que se
// exibe sob uma amostra NÃO foi computado nela: taxonomia GTDB, defesa, AMR,
// plasmídeo, KEGG, CAZy e taxonomia viral saem do REPRESENTANTE do cluster,
// que pode ter vindo de outra amostra ou de um grupo de co-montagem.
//
// Exibir isso sem marca seria fazer uma afirmação que a pipeline nunca fez:
// "o MAG da P01 é Pseudomonadota" quando quem foi classificado foi outro
// genoma. A marca é obrigatória, e o representante viaja nela.
export function Herdado({ de, o = 'valor' }) {
  return (
    <span className="prov prov--herdado" data-inherited-from={de}
          title={`${o} computado em ${de}, o representante do cluster — não neste genoma`}>
      herdado de <code>{de}</code>
    </span>
  );
}

// O par da marca acima. Existe para que "medido aqui" seja uma afirmação
// explícita e não a ausência de outra: QC, montagem, mapeamento, binning,
// CheckM2, GUNC, abundância e diversidade SÃO da amostra.
export function Medido({ o = 'valor' }) {
  return (
    <span className="prov prov--medido" data-measured="true"
          title={`${o} medido nesta amostra`}>
      medido aqui
    </span>
  );
}
