/* Expose the existing archive to supported agent browsers. */
(() => {
  if (window.hljs) window.hljs.highlightAll();
  const index = document.getElementById('journal-index');
  const context = document.modelContext;
  if (!index || !context?.registerTool) return;
  const posts = JSON.parse(index.textContent);
  const topics = [...new Set(posts.map(post => post.topic))];
  const lifecycle = new AbortController();
  window.addEventListener('pagehide', () => lifecycle.abort(), { once: true });
  try {
    Promise.resolve(context.registerTool({
      name: 'list_journal_entries',
      title: 'Read the journal index',
      description: 'List the published journal entries, optionally limited to one of the visible archive topics.',
      inputSchema: { type: 'object', properties: { topic: { type: 'string', enum: topics } }, additionalProperties: false },
      annotations: { readOnlyHint: true, untrustedContentHint: true },
      execute(input) {
        if (!input || typeof input !== 'object' || Array.isArray(input) || Object.keys(input).some(key => key !== 'topic')) throw new Error('Expected an object with an optional topic.');
        if (input.topic !== undefined && !topics.includes(input.topic)) throw new Error('Choose an existing archive topic.');
        const entries = input.topic ? posts.filter(post => post.topic === input.topic) : posts;
        return { count: entries.length, entries };
      }
    }, { signal: lifecycle.signal })).catch(() => {});
  } catch { /* Normal article navigation remains available. */ }
})();
