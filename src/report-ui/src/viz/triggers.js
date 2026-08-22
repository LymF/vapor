// Os gatilhos numericos do REPORT_VIZ_GUIDE.md §4. Vivem DENTRO das formas:
// no report antigo dependiam de o autor do grafico lembrar de chamar o helper
// certo, e o padrao so ficava correto por disciplina.
export const TRIGGERS = {
  manySamples: 12,
  densityMinN: 20,
  manyGroups: 8,
  denseScatter: 500,
  networkNodes: 150,
  tableRows: 200,
};

export function pickDistributionForm(n) {
  // Uma curva de densidade sobre poucos pontos afirma uma distribuicao
  // continua que o dado nao tem: a forma vem da largura de banda.
  return n < TRIGGERS.densityMinN ? 'strip' : 'density';
}

export function pickAxisOrientation(nSamples) {
  return nSamples > TRIGGERS.manySamples ? 'horizontal' : 'vertical';
}
