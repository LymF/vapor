// Paleta validada com scripts/validate_palette.js nos modos claro e escuro.
// Oito e o teto: uma nona cor reprova sob CVD (REPORT_VIZ_GUIDE.md §6).
export const PAL = [
  '#0d9488', '#d97706', '#7c3aed', '#0891b2',
  '#16a34a', '#db2777', '#9333ea', '#ef4444',
];

// Neutro para "Other"/"Unknown". Nunca e um slot de identidade categorica.
export const PAL_MUTED = '#64748b';

export function foldOther(counts, max = 8) {
  const ordenado = Object.entries(counts).sort((a, b) => b[1] - a[1]);
  if (ordenado.length <= max) return ordenado;
  const cabeca = ordenado.slice(0, max - 1);
  const cauda = ordenado.slice(max - 1).reduce((soma, [, v]) => soma + v, 0);
  return [...cabeca, ['Other', cauda]];
}
