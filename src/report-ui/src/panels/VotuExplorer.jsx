// Explorador de vOTU: um seletor + a trilha genomica em coordenadas reais
// (bp) do vOTU escolhido -- PHROGs, AMGs e anti-defesa nas posicoes que o
// Prodigal deu de graca. Desenhar ordem de gene em vez de posicao real e
// anti-padrao declarado no docs/REPORT_VIZ_GUIDE.md.
import { useState } from 'react';
import { GenomeTrack } from '../charts/GenomeTrack.jsx';

// Mapeia o `kind` cru (vindo de phrog_category/category do pharokka) para
// cor + rotulo em portugues. Um kind desconhecido ainda desenha (fallback do
// proprio GenomeTrack para PAL_MUTED + o nome cru), so nao ganha traducao.
const KINDS = {
  phrog: { label: 'PHROG', color: '#0d9488' },
  amg: { label: 'AMG', color: '#d97706' },
  antidefense: { label: 'Anti-defesa', color: '#7c3aed' },
  'anti-defense': { label: 'Anti-defesa', color: '#7c3aed' },
};

export function VotuExplorer({ vOTUs = [] }) {
  const [selecionado, setSelecionado] = useState(vOTUs[0]?.votu_id ?? '');
  const vOTU = vOTUs.find((v) => v.votu_id === selecionado) ?? vOTUs[0];

  if (!vOTUs.length) {
    return <p className="empty">Sem anotação disponível para esta rodada.</p>;
  }

  return (
    <div>
      <label className="nav__filter">
        <span>vOTU</span>
        <select value={vOTU?.votu_id ?? ''} onChange={(e) => setSelecionado(e.target.value)}>
          {vOTUs.map((v) => (
            <option key={v.votu_id} value={v.votu_id}>
              {v.votu_id} ({v.length.toLocaleString('pt-BR')} bp)
            </option>
          ))}
        </select>
      </label>
      <GenomeTrack length={vOTU?.length} features={vOTU?.features ?? []} kinds={KINDS} />
    </div>
  );
}
