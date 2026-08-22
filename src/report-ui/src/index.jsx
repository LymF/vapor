import { createRoot } from 'react-dom/client';

export function App({ data }) {
  return <h1>{data?.run?.title ?? 'VAPOR'}</h1>;
}

const el = typeof document !== 'undefined' && document.getElementById('vapor-root');
if (el) createRoot(el).render(<App data={window.VAPOR_DATA ?? {}} />);
