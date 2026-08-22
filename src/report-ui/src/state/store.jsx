import { createContext, useContext, useMemo, useState } from 'react';

const Ctx = createContext(null);

export const TODAS = '__all__';

export function ReportProvider({ data, children }) {
  const [sample, setSample] = useState(TODAS);
  const [tab, setTab] = useState('overview');
  const samples = data?.run?.samples ?? [];
  const valor = useMemo(
    () => ({ data, sample, setSample, samples, tab, setTab }),
    [data, sample, samples, tab],
  );
  return <Ctx.Provider value={valor}>{children}</Ctx.Provider>;
}

export function useReport() {
  const v = useContext(Ctx);
  if (!v) throw new Error('useReport fora de ReportProvider');
  return v;
}
