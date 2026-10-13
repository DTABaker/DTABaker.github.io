(function () {
  'use strict';
  const G = globalThis.ChemGenerators;
  const H = globalThis.ChemCore;

  const nonelectrolytes = [
    { name: 'urea', formula: 'CH₄N₂O', mm: 60.06, i: 1 },
    { name: 'glucose', formula: 'C₆H₁₂O₆', mm: 180.16, i: 1 }
  ];
  const salts = [
    { name: 'sodium chloride', formula: 'NaCl', mm: 58.44, i: 2, dissolution: 'NaCl(s) → Na⁺(aq) + Cl⁻(aq)' },
    { name: 'magnesium chloride', formula: 'MgCl₂', mm: 95.21, i: 3, dissolution: 'MgCl₂(s) → Mg²⁺(aq) + 2 Cl⁻(aq)' },
    { name: 'calcium chloride', formula: 'CaCl₂', mm: 110.98, i: 3, dissolution: 'CaCl₂(s) → Ca²⁺(aq) + 2 Cl⁻(aq)' }
  ];

  function freezing(rng, electrolyte) {
    const solute = rng.pick(electrolyte ? salts : nonelectrolytes);
    const soluteMass = rng.int(electrolyte ? 50 : 100, electrolyte ? 200 : 450) / 100;
    const solventMass = rng.int(20, 50) * 10;
    const kg = solventMass / 1000;
    const Kf = 1.86;
    const moles = soluteMass / solute.mm;
    const molality = moles / kg;
    const delta = solute.i * Kf * molality;
    const tf = -delta;
    const soluteText = soluteMass.toFixed(3);
    const waterText = solventMass.toFixed(1);
    const particleStep = electrolyte
      ? `<p>${solute.dissolution}</p><p>Each formula unit gives ${solute.i} ions, so the ideal van ’t Hoff factor is <strong>i = ${solute.i}</strong>.</p>`
      : `<p>${solute.name} dissolves as intact molecules. It is a nonelectrolyte, so <strong>i = 1</strong>.</p>`;
    return {
      prompt: `<p>Dissolve <strong>${soluteText} g ${solute.name} (${solute.formula})</strong> in <strong>${waterText} g water</strong>. The solute molar mass is <strong>${solute.mm.toFixed(2)} g/mol</strong>. For water, <strong>K<sub>f</sub> = 1.86 °C·kg/mol</strong> and the pure-solvent freezing point is <strong>0.00 °C</strong> at 1 atm.</p><p>Using the ideal dilute-solution model${electrolyte ? ' and complete dissociation' : ''}, calculate the <strong>positive freezing-point depression ΔT<sub>f</sub></strong> and the <strong>solution freezing point T<sub>f</sub></strong>.</p>`,
      fields: [
        H.num('delta', 'Freezing-point depression ΔTf (positive)', delta, '°C'),
        H.num('tf', 'Solution freezing point Tf', tf, '°C')
      ],
      hint: 'Use ΔTf = i Kf m, where m = moles of solute / kilograms of solvent. Then Tf = 0.00 °C − ΔTf.',
      model: '<p>Ideal dilute aqueous solution; solute is nonvolatile. The calculation uses the stated molar mass and ideal particle count. Real electrolyte solutions can deviate from the ideal van ’t Hoff factor.</p>',
      steps: [
        { title: 'Convert the solute mass to moles', body: `<p>n = mass / molar mass</p><p>n = ${soluteText} g / (${solute.mm.toFixed(2)} g/mol) ≈ <strong>${H.fmt(moles)} mol</strong>.</p>` },
        { title: 'Convert the solvent mass to kilograms', body: `<p>${waterText} g H₂O × (1 kg / 1000 g) = <strong>${kg.toFixed(3)} kg H₂O</strong>.</p><p>Molality uses the mass of the solvent alone.</p>` },
        { title: 'Find the molality', body: `<p>m = n / mass of solvent in kg</p><p>m = [${soluteText} / ${solute.mm.toFixed(2)}] mol / ${kg.toFixed(3)} kg ≈ <strong>${H.fmt(molality)} mol/kg</strong>.</p>` },
        { title: 'Count the dissolved particles', body: particleStep },
        { title: 'Calculate the positive depression', body: `<p>ΔT<sub>f</sub> = i K<sub>f</sub> m</p><p>ΔT<sub>f</sub> = ${solute.i} × (1.86 °C·kg/mol) × [${soluteText} / ${solute.mm.toFixed(2)} / ${kg.toFixed(3)}] mol/kg ≈ <strong>${H.sigfig(delta)} °C</strong>.</p><p>The kg and mol units cancel, leaving °C. ΔT<sub>f</sub> is the positive amount by which the freezing point falls.</p>` },
        { title: 'Find the solution freezing point', body: `<p>T<sub>f,solution</sub> = T<sub>f,pure</sub> − ΔT<sub>f</sub></p><p>T<sub>f,solution</sub> = 0.00 °C − ${H.sigfig(delta)} °C ≈ <strong>${H.sigfig(tf)} °C</strong>.</p><p>Keep unrounded values until the final calculation.</p>` }
      ],
      errorHints: [
        { field: 'delta', value: -delta, message: 'ΔTf is defined here as the positive depression. The solution freezing point is the negative temperature.' },
        { field: 'tf', value: delta, message: 'Freezing is depressed below 0 °C: subtract the positive ΔTf from 0.00 °C.' },
        { field: 'delta', value: delta / 1000, message: 'Molality needs kilograms of solvent. Convert the water mass from grams to kilograms before dividing.' },
        ...(electrolyte ? [{ field: 'delta', value: delta / solute.i, message: `Include all dissolved ions: the ideal van ’t Hoff factor for ${solute.formula} is ${solute.i}.` }] : [])
      ],
      data: { kind: 'freezing', electrolyte, solute: solute.formula, soluteMass, solventMass, masses: { solute: soluteMass, solvent: solventMass }, mm: solute.mm, kg, i: solute.i, Kf, moles, molality, delta, tf }
    };
  }

  G['freezing-nonelectrolyte'] = {
    id: 'freezing-nonelectrolyte',
    group: 'Freezing-point depression',
    title: 'Nonelectrolyte freezing point',
    description: 'Convert mass to molality and calculate ΔTf for urea or glucose in water.',
    slides: [3, 4, 5, 6, 11, 12, 13, 14],
    generate(rng) { return freezing(rng, false); }
  };
  G['freezing-electrolyte'] = {
    id: 'freezing-electrolyte',
    group: 'Freezing-point depression',
    title: 'Electrolyte freezing point',
    description: 'Count ions, calculate the depression, and find the solution freezing point.',
    slides: [7, 8, 9, 10],
    generate(rng) { return freezing(rng, true); }
  };

  const titrations = [
    {
      id: 'titration-strong', title: 'KOH titrated with HCl', slides: [36, 37, 38, 39],
      sample: 'KOH', titrant: 'HCl', sampleCoeff: 1, titrantCoeff: 1,
      reaction: 'KOH(aq) + HCl(aq) → KCl(aq) + H₂O(l)',
      explanation: 'One mole of HCl neutralizes one mole of KOH.',
      model: 'Complete neutralization of a strong monoprotic acid and a strong base; the measured endpoint is treated as the equivalence point.'
    },
    {
      id: 'titration-weak', title: 'Acetic acid titrated with NaOH', slides: [40, 41, 42],
      sample: 'CH₃COOH', sampleName: 'acetic acid (CH₃COOH)', titrant: 'NaOH', sampleCoeff: 1, titrantCoeff: 1,
      reaction: 'CH₃COOH(aq) + NaOH(aq) → CH₃COONa(aq) + H₂O(l)',
      explanation: 'Acetic acid is monoprotic: one mole reacts with one mole of NaOH. Its being a weak acid does not change the neutralization ratio.',
      model: 'Complete neutralization of acetic acid by NaOH; the endpoint is treated as the equivalence point. This determines the initial analytical concentration, not the equilibrium concentration of undissociated acid. The equivalence-point pH need not be 7.'
    },
    {
      id: 'titration-diprotic', title: 'H₂SO₄ titrated with KOH', slides: [43, 44, 45],
      sample: 'H₂SO₄', titrant: 'KOH', sampleCoeff: 1, titrantCoeff: 2,
      reaction: 'H₂SO₄(aq) + 2 KOH(aq) → K₂SO₄(aq) + 2 H₂O(l)',
      explanation: 'Complete neutralization removes both acidic protons: two moles of KOH react with one mole of H₂SO₄.',
      endpoint: 'the complete-neutralization equivalence point (both acidic protons neutralized)',
      model: 'The titration reaches complete neutralization of both protons of H₂SO₄. Calculate the original analytical molarity of H₂SO₄ using the 1:2 stoichiometric ratio.'
    },
    {
      id: 'titration-dihydroxide', title: 'Ba(OH)₂ titrated with HCl', slides: [46, 47, 48],
      sample: 'Ba(OH)₂', titrant: 'HCl', sampleCoeff: 1, titrantCoeff: 2,
      reaction: 'Ba(OH)₂(aq) + 2 HCl(aq) → BaCl₂(aq) + 2 H₂O(l)',
      explanation: 'Each formula unit of Ba(OH)₂ supplies two hydroxide ions, so two moles of HCl neutralize one mole of Ba(OH)₂.',
      model: 'Dissolved Ba(OH)₂ is completely neutralized by HCl; the endpoint is treated as the equivalence point. The requested molarity refers to Ba(OH)₂ formula units, not hydroxide ions.'
    }
  ];

  function titration(rng, spec) {
    const diluteBarium = spec.id === 'titration-dihydroxide';
    const sampleML = rng.pick(diluteBarium ? [25, 30, 40] : [10, 15, 20, 25, 30, 40]);
    const titrantML = rng.int(120, diluteBarium ? 300 : 450) / 10;
    const titrantM = rng.int(diluteBarium ? 25 : 50, diluteBarium ? 75 : 250) / 1000;
    const titrantL = titrantML / 1000;
    const sampleL = sampleML / 1000;
    const ratio = spec.sampleCoeff / spec.titrantCoeff;
    const titrantMoles = titrantM * titrantL;
    const sampleMoles = titrantMoles * ratio;
    const sampleM = sampleMoles / sampleL;
    const sampleName = spec.sampleName || spec.sample;
    return {
      prompt: `<p>A <strong>${sampleML.toFixed(1)} mL sample of ${sampleName}</strong> is titrated with <strong>${titrantM.toPrecision(3)} M ${spec.titrant}</strong>. It takes <strong>${titrantML.toFixed(1)} mL</strong> of the titrant to reach ${spec.endpoint || 'the equivalence point'}.</p><p>Find the <strong>original molar concentration of ${spec.sample}</strong> in the sample.</p>`,
      fields: [H.num('concentration', `Original ${spec.sample} concentration`, sampleM, 'M')],
      hint: `Balance the reaction first. The ratio is ${spec.sampleCoeff} mol ${spec.sample} : ${spec.titrantCoeff} mol ${spec.titrant}. Find the titrant moles, convert to sample moles, then divide by the original sample volume in liters.`,
      model: `<p>${spec.model}</p>`,
      steps: [
        { title: 'Balance the neutralization reaction', body: `<p>${spec.reaction}</p><p>${spec.explanation}</p>` },
        { title: 'Convert both volumes to liters', body: `<p>V<sub>titrant</sub> = ${titrantML.toFixed(1)} mL × (1 L / 1000 mL) = <strong>${titrantL.toFixed(4)} L</strong>.</p><p>V<sub>sample</sub> = ${sampleML.toFixed(1)} mL × (1 L / 1000 mL) = <strong>${sampleL.toFixed(4)} L</strong>.</p>` },
        { title: 'Calculate moles of titrant', body: `<p>n<sub>titrant</sub> = M<sub>titrant</sub>V<sub>titrant</sub></p><p>n(${spec.titrant}) = (${titrantM.toPrecision(3)} mol/L)(${titrantL.toFixed(4)} L) ≈ <strong>${H.fmt(titrantMoles)} mol</strong>.</p>` },
        { title: 'Use the mole ratio to find sample moles', body: `<p>n(${spec.sample}) = n(${spec.titrant}) × (${spec.sampleCoeff} mol ${spec.sample} / ${spec.titrantCoeff} mol ${spec.titrant})</p><p>n(${spec.sample}) = [${titrantM.toPrecision(3)} × ${titrantL.toFixed(4)}] × (${spec.sampleCoeff}/${spec.titrantCoeff}) mol ≈ <strong>${H.fmt(sampleMoles)} mol</strong>.</p>` },
        { title: 'Divide by the original sample volume', body: `<p>M<sub>sample</sub> = n<sub>sample</sub> / V<sub>sample</sub></p><p>M(${spec.sample}) = [${titrantM.toPrecision(3)} × ${titrantL.toFixed(4)} × (${spec.sampleCoeff}/${spec.titrantCoeff})] mol / ${sampleL.toFixed(4)} L ≈ <strong>${H.sigfig(sampleM)} mol/L = ${H.sigfig(sampleM)} M</strong>.</p><p>Use the original ${sampleML.toFixed(1)} mL sample volume because the question asks for its concentration before titration. Keep intermediate values unrounded.</p>` }
      ],
      errorHints: [
        ...(spec.titrantCoeff === 2 ? [
          { field: 'concentration', value: sampleM * 2, message: `The reaction consumes 2 mol ${spec.titrant} per 1 mol ${spec.sample}. Divide the titrant moles by 2.` },
          { field: 'concentration', value: sampleM * 4, message: 'The sample has half as many moles as the titrant at equivalence. Divide by the stoichiometric factor of 2 instead of multiplying by it.' }
        ] : []),
        { field: 'concentration', value: sampleMoles / (sampleL + titrantL), message: 'Use the original sample volume to find its original concentration. Adding the titrant volume gives a concentration in the final mixture instead.' },
        { field: 'concentration', value: sampleM / 1000, message: 'Molarity is mol/L. Convert the original sample volume from mL to L before dividing the moles.' }
      ],
      data: { kind: 'titration', sample: spec.sample, titrant: spec.titrant, titrantM, titrantML, sampleML, ratio, sampleCoeff: spec.sampleCoeff, titrantCoeff: spec.titrantCoeff, titrantMoles, sampleMoles, sampleM }
    };
  }

  for (const spec of titrations) {
    G[spec.id] = {
      id: spec.id,
      group: 'Acid–base titrations',
      title: spec.title,
      description: `Find the original sample concentration using a ${spec.sampleCoeff}:${spec.titrantCoeff} sample-to-titrant mole ratio.`,
      slides: spec.slides,
      generate(rng) { return titration(rng, spec); }
    };
  }
})();
