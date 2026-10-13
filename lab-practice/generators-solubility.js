(function () {
  'use strict';

  const G = globalThis.ChemGenerators = globalThis.ChemGenerators || {};
  const H = globalThis.ChemCore;
  const group = 'Solubility & Ksp';
  const f = (value, sig = 3) => H.sigfig(value, sig);
  const ag = 'Ag<sup>+</sup>';
  const cl = 'Cl<sup>−</sup>';
  const agclReaction = 'AgCl(s) ⇌ Ag<sup>+</sup>(aq) + Cl<sup>−</sup>(aq)';
  const baseModel = '<p>Use a saturated solution at equilibrium with excess solid at a fixed temperature. Treat the dissolved ions as ideal and neglect side reactions. The activity of a pure solid is 1, so the solid is absent from the solubility-product expression. Strictly, <var>K</var><sub>sp</sub> is dimensionless: each dissolved-ion concentration is divided by the standard concentration <var>c</var><sup>°</sup> = 1 mol L<sup>−1</sup>. The shorthand calculations below use the numerical concentrations in mol L<sup>−1</sup> (M).</p><p>The generated measurements and constants are hypothetical practice data, not published values for these salts.</p>';

  function iceTable(species, initial, changes, equilibrium) {
    const row = (label, cells) => '<tr><th scope="row">' + label + '</th>' + cells.map(cell => '<td>' + cell + '</td>').join('') + '</tr>';
    return '<div class="table-scroll"><table><caption>ICE table: dissolved-ion concentrations in mol L<sup>−1</sup> (M)</caption><thead><tr><th scope="col">Stage</th>' + species.map(name => '<th scope="col">' + name + '</th>').join('') + '</tr></thead><tbody>' + row('Initial', initial) + row('Change', changes) + row('Equilibrium', equilibrium) + '</tbody></table></div><p><var>s</var> is the molar solubility: moles of dissolved salt formula units per litre of solution. The solid is present throughout but is omitted from this concentration table.</p>';
  }

  function registerPure(id, title, description, slides, salt, reaction, species, coefficients, makeS, expression, substitution, extraModel) {
    G[id] = {
      id, group, title, description, slides,
      generate(rng) {
        const s = makeS(rng);
        const concentrations = coefficients.map(coefficient => coefficient * s);
        const Ksp = concentrations.reduce((product, concentration, index) => product * Math.pow(concentration, coefficients[index]), 1);
        const symbols = coefficients.map(coefficient => coefficient === 1 ? '<var>s</var>' : coefficient + '<var>s</var>');
        return {
          prompt: '<p>In a hypothetical measurement, the molar solubility of ' + salt + ' in initially pure water is <strong>' + f(s) + ' M</strong>. The solution is saturated and remains in contact with excess solid at equilibrium. Calculate <var>K</var><sub>sp</sub> using the stated ideal model.</p>',
          fields: [H.num('Ksp', 'Solubility product, Ksp', Ksp, 'dimensionless')],
          hint: 'First balance the dissolution reaction. Multiply the molar solubility by each ion’s stoichiometric coefficient, then raise each ion concentration to that coefficient in Ksp.',
          model: baseModel + (extraModel || ''),
          steps: [
            { title: 'Write the balanced dissolution reaction', body: '<p class="equation">' + reaction + '</p><p>Each mole of dissolved salt produces ' + coefficients.map((coefficient, index) => coefficient + ' mol of ' + species[index]).join(' and ') + '.</p>' },
            { title: 'Build the ICE table', body: iceTable(species, ['0', '0'], symbols.map(symbol => '+' + symbol), symbols) },
            { title: 'Use the equilibrium ion concentrations', body: '<p>' + species.map((name, index) => '[' + name + '] = ' + symbols[index] + ' = <strong>' + f(concentrations[index], 4) + ' M</strong>').join('; ') + '.</p><p>The molar solubility refers to dissolved formula units; it is not automatically the concentration of every ion.</p>' },
            { title: 'Form and evaluate the solubility product', body: '<p class="equation">' + expression + '</p><p class="equation">' + substitution(s) + ' = <strong>' + f(Ksp) + '</strong>.</p><p>The coefficients in the balanced reaction become exponents. Retain extra digits while calculating, then report the result to 3 significant figures.</p>' }
          ],
          data: { kind: id, s, c0: 0, Ksp, stoichiometry: coefficients.slice(), concentrations }
        };
      }
    };
  }

  registerPure(
    'ksp-11', 'Ksp from solubility: 1:1 salt', 'Use a saturated AgCl solution to connect molar solubility and Ksp.', [16, 17, 18, 19],
    'AgCl', agclReaction, [ag, cl], [1, 1], rng => rng.int(100, 300) * 1e-7,
    '<var>K</var><sub>sp</sub> = [Ag<sup>+</sup>][Cl<sup>−</sup>] = <var>s</var><sup>2</sup>',
    s => '<var>K</var><sub>sp</sub> = (' + f(s) + ')<sup>2</sup>'
  );

  registerPure(
    'ksp-12', 'Ksp from solubility: 1:2 salt', 'Account for the two fluoride ions produced by CaF2.', [20, 21, 22, 23],
    'CaF<sub>2</sub>', 'CaF<sub>2</sub>(s) ⇌ Ca<sup>2+</sup>(aq) + 2F<sup>−</sup>(aq)', ['Ca<sup>2+</sup>', 'F<sup>−</sup>'], [1, 2], rng => rng.int(100, 350) * 1e-6,
    '<var>K</var><sub>sp</sub> = [Ca<sup>2+</sup>][F<sup>−</sup>]<sup>2</sup> = (<var>s</var>)(2<var>s</var>)<sup>2</sup> = 4<var>s</var><sup>3</sup>',
    s => '<var>K</var><sub>sp</sub> = 4(' + f(s) + ')<sup>3</sup>'
  );

  registerPure(
    'ksp-21', 'Ksp from solubility: 2:1 salt', 'Account for the two silver ions in an idealized Ag2CO3 model.', [24, 25, 26, 27],
    'Ag<sub>2</sub>CO<sub>3</sub>', 'Ag<sub>2</sub>CO<sub>3</sub>(s) ⇌ 2Ag<sup>+</sup>(aq) + CO<sub>3</sub><sup>2−</sup>(aq)', [ag, 'CO<sub>3</sub><sup>2−</sup>'], [2, 1], rng => rng.int(100, 300) * 1e-6,
    '<var>K</var><sub>sp</sub> = [Ag<sup>+</sup>]<sup>2</sup>[CO<sub>3</sub><sup>2−</sup>] = (2<var>s</var>)<sup>2</sup>(<var>s</var>) = 4<var>s</var><sup>3</sup>',
    s => '<var>K</var><sub>sp</sub> = 4(' + f(s) + ')<sup>3</sup>',
    '<p><strong>Carbonate simplification:</strong> neglect carbonate protonation, carbonate hydrolysis, and silver complex formation. Therefore this exercise sets [CO<sub>3</sub><sup>2−</sup>] = <var>s</var>. In real water, total dissolved carbonate need not equal the free carbonate concentration used in <var>K</var><sub>sp</sub>.</p>'
  );

  G['ksp-common'] = {
    id: 'ksp-common', group,
    title: 'Ksp with an initial common ion',
    description: 'Include the initial silver concentration when calculating Ksp from a measured solubility.',
    slides: [28, 29, 30],
    generate(rng) {
      const c0 = rng.int(100, 400) * rng.pick([1e-7, 1e-6]);
      const s = rng.int(100, 500) * 1e-8;
      const silver = c0 + s;
      const Ksp = silver * s;
      return {
        prompt: '<p>A hypothetical solution initially contains <strong>' + f(c0) + ' M Ag<sup>+</sup></strong> from AgNO<sub>3</sub> and no Cl<sup>−</sup>. Excess AgCl(s) is added. After saturation and equilibrium, its measured molar solubility is <strong>' + f(s) + ' M</strong>. Calculate <var>K</var><sub>sp</sub> using the ideal model. Include the additional Ag<sup>+</sup> released by the dissolving AgCl.</p>',
        fields: [H.num('Ksp', 'Solubility product, Ksp', Ksp, 'dimensionless')],
        hint: 'The starting silver concentration is not zero. At equilibrium, [Ag+] = c₀ + s and [Cl−] = s.',
        model: baseModel + '<p>AgNO<sub>3</sub> is treated as fully dissociated, nitrate is a spectator ion, and dilution from adding solid is neglected.</p>',
        steps: [
          { title: 'Write the dissolution reaction', body: '<p class="equation">' + agclReaction + '</p><p>Dissolving <var>s</var> M of AgCl increases each ion concentration by <var>s</var> M.</p>' },
          { title: 'Include the initial common ion in the ICE table', body: iceTable([ag, cl], ['<var>c</var><sub>0</sub> = ' + f(c0), '0'], ['+<var>s</var>', '+<var>s</var>'], ['<var>c</var><sub>0</sub> + <var>s</var>', '<var>s</var>']) },
          { title: 'Find both equilibrium concentrations', body: '<p>[Ag<sup>+</sup>] = ' + f(c0) + ' + ' + f(s) + ' = <strong>' + f(silver, 5) + ' M</strong>.</p><p>[Cl<sup>−</sup>] = <var>s</var> = <strong>' + f(s) + ' M</strong>.</p>' },
          { title: 'Calculate Ksp without dropping the solubility term', body: '<p class="equation"><var>K</var><sub>sp</sub> = [Ag<sup>+</sup>][Cl<sup>−</sup>] = (<var>c</var><sub>0</sub> + <var>s</var>)<var>s</var></p><p class="equation"><var>K</var><sub>sp</sub> = (' + f(c0) + ' + ' + f(s) + ')(' + f(s) + ') = <strong>' + f(Ksp) + '</strong>.</p><p>The common-ion term is added before multiplication. Using <var>s</var><sup>2</sup> would incorrectly discard the silver already present.</p>' }
        ],
        data: { kind: 'ksp-common', s, c0, Ksp, stoichiometry: [1, 1], concentrations: [silver, s] }
      };
    }
  };

  G['solubility-common'] = {
    id: 'solubility-common', group,
    title: 'Solubility with a common ion',
    description: 'Solve for AgCl solubility and check a common-ion approximation.',
    slides: [31, 32, 33, 34],
    generate(rng) {
      const Ksp = rng.int(100, 400) * 1e-12;
      const c0 = rng.int(100, 400) * rng.pick([1e-8, 1e-7, 1e-6, 1e-5]);
      // Rationalized positive root avoids subtracting nearly equal numbers.
      const s = 2 * Ksp / (c0 + Math.sqrt(c0 * c0 + 4 * Ksp));
      const approximateS = Ksp / c0;
      const fraction = s / c0;
      const approximationValid = fraction < 0.05;
      const percent = fraction * 100;
      return {
        prompt: '<p>For this hypothetical practice dataset, AgCl has <strong><var>K</var><sub>sp</sub> = ' + f(Ksp) + '</strong>. A solution initially contains <strong>' + f(c0) + ' M Ag<sup>+</sup></strong> from AgNO<sub>3</sub> and no Cl<sup>−</sup>. Excess AgCl(s) is added and allowed to reach equilibrium.</p><p>Find the molar solubility <var>s</var> (equal to the equilibrium [Cl<sup>−</sup>]) using the exact quadratic solution. Then decide whether dropping <var>s</var> from <var>c</var><sub>0</sub> + <var>s</var> meets the criterion <strong><var>s</var>/<var>c</var><sub>0</sub> &lt; 0.05</strong>, using your exact result.</p>',
        fields: [
          H.num('solubility', 'Molar solubility, s', s, 'mol/L'),
          H.choice('approximation', 'Does s/c₀ < 0.05?', approximationValid ? 'yes' : 'no', [{ value: 'yes', label: 'Yes — the criterion is met' }, { value: 'no', label: 'No — keep the solubility term' }])
        ],
        hint: 'Set Ksp = (c₀ + s)s. Rearrange to s² + c₀s − Ksp = 0 and keep the positive root. Check s/c₀ against 0.05.',
        model: baseModel + '<p>AgNO<sub>3</sub> is fully dissociated and nitrate is a spectator ion. The 5% criterion is a stated approximation rule for this exercise, not a physical boundary.</p>',
        steps: [
          { title: 'Write the reaction and ICE table', body: '<p class="equation">' + agclReaction + '</p>' + iceTable([ag, cl], ['<var>c</var><sub>0</sub> = ' + f(c0), '0'], ['+<var>s</var>', '+<var>s</var>'], ['<var>c</var><sub>0</sub> + <var>s</var>', '<var>s</var>']) },
          { title: 'Set up the solubility-product equation', body: '<p class="equation"><var>K</var><sub>sp</sub> = [Ag<sup>+</sup>][Cl<sup>−</sup>] = (<var>c</var><sub>0</sub> + <var>s</var>)<var>s</var></p><p class="equation">' + f(Ksp) + ' = (' + f(c0) + ' + <var>s</var>)<var>s</var></p><p class="equation"><var>s</var><sup>2</sup> + (' + f(c0) + ')<var>s</var> − ' + f(Ksp) + ' = 0</p><p>Here the symbols in the arithmetic denote numerical concentration values in M.</p>' },
          { title: 'Take the positive quadratic root', body: '<p class="equation"><var>s</var> = [−<var>c</var><sub>0</sub> + √(<var>c</var><sub>0</sub><sup>2</sup> + 4<var>K</var><sub>sp</sub>)]/2</p><p>The equivalent form below avoids subtracting nearly equal values:</p><p class="equation"><var>s</var> = 2<var>K</var><sub>sp</sub> / [<var>c</var><sub>0</sub> + √(<var>c</var><sub>0</sub><sup>2</sup> + 4<var>K</var><sub>sp</sub>)]</p><p class="equation"><var>s</var> = 2(' + f(Ksp) + ') / [' + f(c0) + ' + √((' + f(c0) + ')<sup>2</sup> + 4(' + f(Ksp) + '))] = <strong>' + f(s) + ' M</strong>.</p><p>Reject the negative root because a molar solubility cannot be negative.</p>' },
          { title: 'Check the common-ion approximation', body: '<p>If <var>s</var> is small relative to <var>c</var><sub>0</sub>, then <var>K</var><sub>sp</sub> ≈ <var>c</var><sub>0</sub><var>s</var>, giving <var>s</var><sub>approx</sub> = <var>K</var><sub>sp</sub>/<var>c</var><sub>0</sub> = ' + f(Ksp) + ' / ' + f(c0) + ' = <strong>' + f(approximateS) + ' M</strong>.</p><p>Using the exact result, 100 × <var>s</var>/<var>c</var><sub>0</sub> = <strong>' + f(percent) + '%</strong>. This is ' + (approximationValid ? 'less than 5%, so the stated criterion is met.' : 'at least 5%, so the stated criterion is not met; retain the exact result.') + '</p><p>For this 1:1 model, the approximation overestimates <var>s</var>; its relative error compared with the exact value is also <var>s</var>/<var>c</var><sub>0</sub>. The exact quadratic answer is valid in either case.</p>' },
          { title: 'Check the equilibrium product', body: '<p>[Ag<sup>+</sup>] = <var>c</var><sub>0</sub> + <var>s</var> = ' + f(c0 + s, 5) + ' M; [Cl<sup>−</sup>] = ' + f(s, 5) + ' M.</p><p>Using unrounded values, (<var>c</var><sub>0</sub> + <var>s</var>)<var>s</var> = ' + f((c0 + s) * s, 4) + ', recovering the given <var>K</var><sub>sp</sub>.</p>' }
        ],
        data: { kind: 'solubility-common', s, c0, Ksp, stoichiometry: [1, 1], concentrations: [c0 + s, s], approximateS, fraction, approximationValid }
      };
    }
  };
})();
