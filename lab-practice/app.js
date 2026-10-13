(function () {
  'use strict';
  const G = globalThis.ChemGenerators, H = globalThis.ChemCore;
  const $ = id => document.getElementById(id);
  const ids = ['freezing-nonelectrolyte','freezing-electrolyte','ksp-11','ksp-12','ksp-21','ksp-common','solubility-common','titration-strong','titration-weak','titration-diprotic','titration-dihydroxide','molarity','molality','lewis'];
  if (!H || typeof H.sigfig !== 'function' || !G || ids.some(id => typeof G[id]?.generate !== 'function')) {
    $('load-message').textContent = 'The study tool could not finish loading. Check your connection and reload the page to load all 14 topics.';
    return;
  }
  const storageKey = 'genchem-lab-practice-v1';
  let progress = {}, storageAvailable = true, current, attempted = false, revealed = 0, mixed = false, toastTimer, resetPending = false, resetTimer;
  function readProgress() {
    // Keep access errors distinct from damaged JSON so a denied read never
    // discards the progress already held by this session.
    const raw = localStorage.getItem(storageKey);
    let stored;
    try { stored = JSON.parse(raw || '{}'); } catch (_) { stored = {}; }
    if (!stored || typeof stored !== 'object' || Array.isArray(stored)) stored = {};
    const valid = {};
    for (const id of ids) {
      const row = stored[id];
      if (row && Number.isSafeInteger(row.attempted) && row.attempted >= 0 && Number.isSafeInteger(row.firstCorrect) && row.firstCorrect >= 0 && row.firstCorrect <= row.attempted) valid[id] = {attempted:row.attempted, firstCorrect:row.firstCorrect, lastCorrect:typeof row.lastCorrect === 'boolean' ? row.lastCorrect : undefined};
    }
    return valid;
  }
  function refreshProgress() {
    if (storageAvailable) try { progress = readProgress(); } catch (_) { storageAvailable = false; }
  }
  try {
    progress = readProgress();
  } catch (_) { storageAvailable = false; }
  const el = (tag, text, className) => { const node = document.createElement(tag); if (text !== undefined) node.textContent = text; if (className) node.className = className; return node; };
  function toast(message, delay = 3500) { clearTimeout(toastTimer); $('toast').textContent = message; $('toast').hidden = false; toastTimer = setTimeout(() => { $('toast').hidden = true; }, delay); }
  function randomSeed() { if (globalThis.crypto?.getRandomValues) return crypto.getRandomValues(new Uint32Array(1))[0]; return (Date.now() ^ Math.floor(Math.random()*4294967296)) >>> 0; }
  function save() { if (storageAvailable) try { localStorage.setItem(storageKey, JSON.stringify(progress)); } catch (_) { storageAvailable = false; } renderProgress(); }
  function renderProgress() {
    const rows = Object.values(progress);
    $('total-practiced').textContent = rows.reduce((sum, row) => sum + row.attempted, 0);
    $('first-try').textContent = rows.reduce((sum, row) => sum + row.firstCorrect, 0) + ' correct on first try';
    $('storage-note').textContent = storageAvailable ? 'Saved only in this browser.' : 'Storage unavailable; progress lasts this session.';
    document.querySelectorAll('.topic-attempts').forEach(node => { const row = progress[node.dataset.topic]; node.textContent = row?.attempted ? `${row.firstCorrect}/${row.attempted}` : ''; node.title = 'Correct on first try / attempted'; });
  }
  function buildTopics() {
    $('topic-count').textContent = `${ids.length} types`;
    const groups = new Map();
    for (const id of ids) { const g = G[id]; if (!groups.has(g.group)) groups.set(g.group, []); groups.get(g.group).push(id); }
    for (const [group, groupIds] of groups) {
      const section = el('section', undefined, 'topic-section'); section.append(el('h3', group));
      for (const id of groupIds) {
        const button = el('button', undefined, 'topic-button'); button.type = 'button'; button.dataset.topic = id; button.title = G[id].description;
        button.append(el('span', G[id].title)); const count = el('span', '', 'topic-attempts'); count.dataset.topic = id; button.append(count);
        const dot = el('span', '', 'topic-dot'); dot.setAttribute('aria-hidden', 'true'); button.append(dot);
        button.addEventListener('click', () => { mixed = false; newProblem(id, randomSeed(), true); }); section.append(button);
      }
      $('topic-list').append(section);
    }
  }
  function newProblem(id, seed, focus = false) {
    cancelReset();
    if (!G[id]) id = ids[0];
    current = { id, seed: seed >>> 0, ...G[id].generate(H.rng(seed)) }; attempted = false; revealed = 0;
    const params = new URLSearchParams({v:H.version, type:id, seed:String(current.seed)});
    // replaceState avoids creating a browser-history entry for every practice item.
    try { history.replaceState(null, '', `${location.pathname}${location.search}#${params}`); } catch (_) { /* file previews can restrict history */ }
    $('topic-group').textContent = mixed ? `MIXED REVIEW / ${G[id].group}` : G[id].group;
    $('problem-heading').textContent = G[id].title;
    $('problem-code').textContent = `#${current.seed.toString(16).padStart(8,'0')}`;
    $('problem-prompt').innerHTML = current.prompt;
    $('model-content').innerHTML = current.model || ''; $('model-details').hidden = !current.model; $('model-details').open = false;
    $('hint-content').innerHTML = current.hint; $('hint-content').hidden = true; $('hint-button').setAttribute('aria-expanded', 'false');
    $('working').value = ''; $('working').closest('details').open = false;
    $('feedback').hidden = true; $('feedback').replaceChildren();
    $('solution').hidden = true; $('solution-steps').replaceChildren(); $('solution-finish').hidden = true;
    $('reveal-next').disabled = true; $('reveal-next').textContent = 'Show first step'; $('reveal-all').hidden = false;
    $('reveal-next-bottom').hidden = true;
    $('solution-instruction').textContent = 'Submit an answer attempt to unlock the worked solution.';
    $('answer-fields').replaceChildren();
    for (const field of current.fields) {
      const wrap = el('div', undefined, 'field'); const label = el('label', field.label); label.htmlFor = `answer-${field.id}`; wrap.append(label);
      const inputWrap = el('div', undefined, 'input-wrap'); let input;
      if (field.type === 'select') {
        input = el('select'); const placeholder = el('option', 'Choose an answer'); placeholder.value = ''; input.append(placeholder);
        // Shuffle choices independently of generation; a shared problem keeps its order.
        const options = field.options.slice(), r = H.rng(current.seed ^ (field.id.length * 7919));
        for (let i = options.length-1; i>0; i--) { const j = r.int(0,i); [options[i],options[j]]=[options[j],options[i]]; }
        for (const option of options) { const o = el('option', option.label); o.value = option.value; input.append(o); }
      } else {
        input = el('input'); input.type = 'text'; input.inputMode = 'text'; input.autocomplete = 'off'; input.spellcheck = false; input.placeholder = 'Your numerical answer';
      }
      input.id = `answer-${field.id}`; input.name = field.id; input.required = true; input.setAttribute('aria-describedby',`note-${field.id} answer-help${field.unit ? ' unit-'+field.id : ''}`);
      input.addEventListener('input', clearFeedback); input.addEventListener('change', clearFeedback);
      inputWrap.append(input); if (field.unit) { const unit = el('span', field.unit, 'unit'); unit.id = `unit-${field.id}`; inputWrap.append(unit); } wrap.append(inputWrap);
      const note = el('p', '', 'field-note'); note.id = `note-${field.id}`; wrap.append(note); $('answer-fields').append(wrap);
    }
    document.querySelectorAll('.topic-button').forEach(button => button.setAttribute('aria-current', String(button.dataset.topic === id)));
    $('mixed-mode').setAttribute('aria-pressed', String(mixed));
    renderProgress(); if (focus) $('problem-heading').focus({preventScroll:false});
  }
  function clearFeedback() {
    $('feedback').hidden = true;
    for (const field of current.fields) {
      $(`answer-${field.id}`).removeAttribute('aria-invalid');
      $(`note-${field.id}`).textContent = '';
    }
  }
  function chooseMixed() {
    refreshProgress();
    const choices = ids.flatMap(id => {
      const row = progress[id];
      const weight = !row?.attempted ? 3 : row.lastCorrect === false ? 5 : 1;
      return Array(weight).fill(id);
    });
    return choices[H.rng(randomSeed()).int(0,choices.length-1)];
  }
  function loadHash() {
    const params = new URLSearchParams(location.hash.slice(1)), id = params.get('type'), rawSeed = params.get('seed'), version = params.get('v');
    const validSeed = /^\d{1,10}$/.test(rawSeed || '') && Number(rawSeed) <= 4294967295;
    if (version && version !== H.version) { toast('This link uses a different problem version. A fresh problem is shown.', 7000); newProblem(ids[0], randomSeed()); return; }
    if (G[id] && validSeed) newProblem(id, Number(rawSeed));
    else { if (location.hash) toast('That problem link is incomplete. Here is a fresh problem.'); newProblem(ids[0], randomSeed()); }
  }
  $('answer-form').addEventListener('submit', event => {
    event.preventDefault(); clearFeedback(); const grades = current.fields.map(field => ({field, ...H.grade(field, $(`answer-${field.id}`).value)}));
    let invalid = null;
    for (const row of grades) {
      const input = $(`answer-${row.field.id}`); input.setAttribute('aria-invalid', String(!row.valid));
      $(`note-${row.field.id}`).textContent = row.valid ? '' : row.field.type === 'select' ? 'Choose an answer before checking.' : 'Enter a finite number without units (for example, -0.125 or 2.4e-8).';
      if (!row.valid && !invalid) invalid = input;
    }
    if (invalid) { invalid.focus(); toast('Complete each answer so your attempt can be checked.'); return; }
    const correct = grades.every(row => row.correct);
    if (!attempted) {
      // Re-read immediately before this update so an older tab cannot restore
      // counts cleared elsewhere or overwrite a completed update in another tab.
      // localStorage has no atomic increment; simultaneous submissions remain
      // an informal-counter limitation, not suitable for recorded grades.
      refreshProgress();
      const row = progress[current.id] || {attempted:0,firstCorrect:0}; row.attempted++; if (correct) row.firstCorrect++; row.lastCorrect = correct; progress[current.id] = row;
      attempted = true; save();
    }
    $('feedback').className = `feedback${correct ? '' : ' needs-work'}`; $('feedback').replaceChildren();
    $('feedback').append(el('h4', correct ? 'Your answers check out.' : 'You have a starting point. Let’s check the method.'));
    const list = el('ul');
    for (const row of grades) {
      const item = el('li', `${row.field.label}: ${row.correct ? 'correct' : 'needs another look'}.`);
      if (!row.correct && row.field.type !== 'select') {
        const match = current.errorHints?.find(h => h.field === row.field.id && Math.abs(row.value - h.value) <= Math.abs(h.value)*0.01);
        if (match) item.append(document.createTextNode(' ' + match.message));
      }
      list.append(item);
    }
    $('feedback').append(list, el('p', correct ? 'Open the worked steps to compare your reasoning, or try fresh values.' : 'You can revise your answers, ask for a hint, or open the worked steps below.'));
    $('feedback').hidden = false; $('reveal-next').disabled = revealed >= current.steps.length;
    $('solution-instruction').textContent = revealed ? `${revealed} of ${current.steps.length} steps shown.` : 'Your attempt is in. Reveal the solution at your pace.';
  });
  function showSteps(all = false) {
    if (!attempted) return;
    const target = all ? current.steps.length : Math.min(revealed + 1, current.steps.length);
    $('solution').hidden = false;
    while (revealed < target) {
      const step = current.steps[revealed++], li = el('li'); li.append(el('h4', step.title)); const body = el('div', undefined, 'step-body'); body.innerHTML = step.body; li.append(body); $('solution-steps').append(li);
    }
    $('solution-instruction').textContent = `${revealed} of ${current.steps.length} steps shown.`;
    $('reveal-next').textContent = revealed >= current.steps.length ? 'All steps shown' : 'Show next step'; $('reveal-next').disabled = revealed >= current.steps.length;
    $('reveal-next-bottom').hidden = revealed >= current.steps.length;
    if (revealed === current.steps.length) {
      $('reveal-all').hidden = true; $('solution-finish').hidden = false; $('answer-summary').replaceChildren(el('strong','Final answers'));
      for (const field of current.fields) {
        const value = field.type === 'select' ? field.options.find(o => o.value === field.answer).label : `${field.tolerance === 0 ? String(field.answer) : H.sigfig(field.answer, 3)}${field.unit ? ' ' + field.unit : ''}`;
        $('answer-summary').append(el('p', `${field.label}: ${value}`));
      }
    }
    const latestStep = $('solution-steps').lastElementChild;
    latestStep.tabIndex = -1; latestStep.focus({preventScroll:true});
    latestStep.scrollIntoView({behavior:window.matchMedia?.('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth',block:'nearest'});
  }
  $('reveal-next').addEventListener('click', () => showSteps()); $('reveal-all').addEventListener('click', () => showSteps(true));
  $('reveal-next-bottom').addEventListener('click', () => showSteps());
  $('new-problem').addEventListener('click', () => newProblem(mixed ? chooseMixed() : current.id, randomSeed(), true));
  $('another-like').addEventListener('click', () => { mixed = false; newProblem(current.id, randomSeed(), true); });
  $('mixed-mode').addEventListener('click', () => { mixed = true; newProblem(chooseMixed(), randomSeed(), true); });
  $('hint-button').addEventListener('click', () => { const show = $('hint-content').hidden; $('hint-content').hidden = !show; $('hint-button').setAttribute('aria-expanded',String(show)); });
  function toggleReference(show) { $('reference').hidden = !show; $('reference-toggle').setAttribute('aria-expanded',String(show)); }
  $('reference-toggle').addEventListener('click', () => toggleReference($('reference').hidden)); $('reference-close').addEventListener('click', () => { toggleReference(false); $('reference-toggle').focus(); });
  $('share-problem').addEventListener('click', async () => {
    try { await navigator.clipboard.writeText(location.href); toast('Problem link copied. The link contains no answers or working.'); }
    catch (_) { toast('Copy this problem’s URL from your browser address bar.', 7000); }
  });
  document.querySelector('.skip-link').addEventListener('click', event => { event.preventDefault(); $('problem-heading').focus(); });
  function cancelReset() { clearTimeout(resetTimer); resetPending = false; $('reset-progress').textContent = 'Reset'; $('cancel-reset').hidden = true; }
  $('cancel-reset').addEventListener('click', cancelReset);
  $('reset-progress').addEventListener('click', () => {
    if (!resetPending) { resetPending = true; $('reset-progress').textContent = 'Confirm reset'; $('cancel-reset').hidden = false; resetTimer = setTimeout(cancelReset,15000); toast('Click Confirm reset to clear this browser’s practice counts, or Cancel to keep them.'); return; }
    progress = {}; cancelReset(); save(); toast('Practice counts cleared.');
  });
  window.addEventListener('hashchange', () => { mixed = false; loadHash(); });
  window.addEventListener('storage', event => {
    if (event.key === storageKey || event.key === null) { refreshProgress(); renderProgress(); }
  });
  buildTopics(); loadHash();
  $('load-status').hidden = true; $('practice-workspace').hidden = false;
})();
