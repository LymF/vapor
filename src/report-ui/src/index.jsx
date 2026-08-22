import { createRoot } from 'react-dom/client';
import { App } from './App.jsx';
import './styles.css';

export { App };

const el = typeof document !== 'undefined' && document.getElementById('vapor-root');
if (el) createRoot(el).render(<App data={window.VAPOR_DATA ?? {}} />);
