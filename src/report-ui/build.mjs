import * as esbuild from 'esbuild';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const out = path.resolve(here, '../../scripts/report/assets/report-ui.js');

await esbuild.build({
  entryPoints: [path.join(here, 'src/index.jsx')],
  bundle: true,
  format: 'iife',
  minify: true,
  jsx: 'automatic',
  target: ['es2020'],
  define: { 'process.env.NODE_ENV': '"production"' },
  loader: { '.jsx': 'jsx' },
  outfile: out,
});

console.log(`[report-ui] bundle escrito em ${out}`);
