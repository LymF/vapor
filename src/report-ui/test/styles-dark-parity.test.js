import fs from 'node:fs';
import path from 'node:path';

function tokens(bloco) {
  return Object.fromEntries([...bloco.matchAll(/(--[\w-]+):\s*([^;]+);/g)].map(m => [m[1], m[2].trim()]));
}

test('os dois blocos de modo escuro definem os mesmos tokens com os mesmos valores', () => {
  const css = fs.readFileSync(path.join(__dirname, '../src/styles.css'), 'utf8');
  const media = css.match(/@media \(prefers-color-scheme: dark\)\s*\{\s*:root:not\(\[data-theme="light"\]\)\s*\{([\s\S]*?)\}/)[1];
  const toggle = css.match(/:root\[data-theme="dark"\]\s*\{([\s\S]*?)\}/)[1];
  expect(tokens(toggle)).toEqual(tokens(media));
});
