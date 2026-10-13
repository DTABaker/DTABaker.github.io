(function () {
  "use strict";
  const G = globalThis.ChemGenerators;
  const H = globalThis.ChemCore;
  const f = (value, sig = 4) => H.fmt(value, sig);
  const shuffle = (items, rng) => {
    const result = items.slice();
    for (let i = result.length - 1; i > 0; i -= 1) {
      const j = rng.int(0, i);
      [result[i], result[j]] = [result[j], result[i]];
    }
    return result;
  };
  const choose = (id, label, answer, options, rng) => H.choice(id, label, answer, shuffle(options, rng));
  const option = (value, label) => ({ value, label });
  const solutes = [
    { name: "sodium chloride", formula: "NaCl", mm: 58.44 },
    { name: "potassium chloride", formula: "KCl", mm: 74.55 },
    { name: "glucose", formula: "C₆H₁₂O₆", mm: 180.16 },
    { name: "urea", formula: "CH₄N₂O", mm: 60.06 },
    { name: "sucrose", formula: "C₁₂H₂₂O₁₁", mm: 342.30 }
  ];
  const concentrationOptions = [
    option("molarity", "Moles of solute per liter of solution"),
    option("molality", "Moles of solute per kilogram of solvent"),
    option("solvent-volume", "Moles of solute per liter of solvent"),
    option("solution-mass", "Moles of solute per kilogram of solution")
  ];

  function concentrationDefinition(kind, rng) {
    const isMolarity = kind === "molarity";
    const name = isMolarity ? "Molarity" : "Molality";
    const denominator = isMolarity ? "liters of the final solution" : "kilograms of the solvent";
    return {
      prompt: `<p><strong>${name}</strong> describes a solution's concentration. Select its definition, then identify the quantity in the denominator.</p>`,
      fields: [
        choose("definition", `Which definition describes ${kind}?`, kind, concentrationOptions, rng),
        choose("denominator", "What belongs in the denominator?", isMolarity ? "solution-volume" : "solvent-mass", [
          option("solution-volume", "Volume of the final solution (L)"),
          option("solvent-volume", "Volume of solvent before mixing (L)"),
          option("solvent-mass", "Mass of solvent alone (kg)"),
          option("solution-mass", "Mass of solute plus solvent (kg)")
        ], rng)
      ],
      steps: [
        { title: "Identify the numerator", body: "<p>Both concentration units use the amount of solute in moles.</p>" },
        { title: "Choose the denominator", body: `<p>${name} divides moles of solute by <strong>${denominator}</strong>.</p><p>${isMolarity ? "M = n(solute) / V(solution), with V in L. The final solution includes solute and solvent." : "m = n(solute) / mass(solvent), with mass in kg. Exclude the solute mass from this denominator."}</p>` },
        { title: "Keep the units distinct", body: "<p>Molarity uses mol/L of solution (M); molality uses mol/kg of solvent (m). Their values cannot generally be converted into each other without additional information, such as solution density and composition.</p>" }
      ],
      hint: "Both start with moles of solute. Decide whether the denominator measures solution volume or solvent mass.",
      data: { variant: "definition", kind, answer: kind }
    };
  }

  G.molarity = {
    id: "molarity", group: "General questions", title: "Molarity", slides: [49],
    description: "Define molarity and calculate concentration from solute mass and final solution volume.",
    generate(rng) {
      if (rng.int(0, 3) === 0) return concentrationDefinition("molarity", rng);
      const solute = rng.pick(solutes);
      const mass = rng.int(100, 1200) / 100;
      const volumeML = rng.int(10, 100) * 10;
      const volumeL = volumeML / 1000;
      const n = mass / solute.mm;
      const answer = n / volumeL;
      return {
        prompt: `<p>Dissolve <strong>${mass.toFixed(2)} g ${solute.name} (${solute.formula})</strong> in water and dilute to a <strong>final solution volume of ${volumeML.toFixed(1)} mL</strong>. Use a molar mass of <strong>${solute.mm.toFixed(2)} g/mol</strong>.</p><p>Calculate the molarity of the dissolved solute.</p>`,
        fields: [H.num("molarity", "Molarity", answer, "mol/L")],
        steps: [
          { title: "Convert solute mass to moles", body: `<p>n = mass / molar mass = (${mass.toFixed(2)} g) / (${solute.mm.toFixed(2)} g/mol) = <strong>${f(n, 6)} mol</strong>.</p>` },
          { title: "Convert the final solution volume", body: `<p>V = ${volumeML.toFixed(1)} mL × (1 L / 1000 mL) = <strong>${f(volumeL)} L</strong>.</p>` },
          { title: "Divide moles by liters of solution", body: `<p>M = n / V = (${f(n, 6)} mol) / (${f(volumeL)} L) ≈ <strong>${H.sigfig(answer, 3)} mol/L</strong> (3 significant figures).</p><p>The denominator is the final solution volume. Use the unrounded mole value in this calculation.</p>` }
        ],
        hint: "Convert grams of solute to moles, convert the final solution volume to liters, and divide.",
        model: "<p>Complete dissolution is assumed. This is solute formula-unit concentration; do not multiply it by the number of ions formed.</p>",
        data: { variant: "numeric", mass, mm: solute.mm, volumeL, answer }
      };
    }
  };

  G.molality = {
    id: "molality", group: "General questions", title: "Molality", slides: [49],
    description: "Define molality and calculate concentration using the mass of solvent alone.",
    generate(rng) {
      if (rng.int(0, 3) === 0) return concentrationDefinition("molality", rng);
      const solute = rng.pick(solutes);
      const mass = rng.int(200, 2000) / 100;
      const solventG = rng.int(1000, 8000) / 10;
      const kg = solventG / 1000;
      const n = mass / solute.mm;
      const answer = n / kg;
      return {
        prompt: `<p>Dissolve <strong>${mass.toFixed(2)} g ${solute.name} (${solute.formula})</strong> in <strong>${solventG.toFixed(1)} g water</strong>. Use a molar mass of <strong>${solute.mm.toFixed(2)} g/mol</strong>.</p><p>Calculate the molality of the dissolved solute.</p>`,
        fields: [H.num("molality", "Molality", answer, "mol/kg")],
        steps: [
          { title: "Convert solute mass to moles", body: `<p>n = mass / molar mass = (${mass.toFixed(2)} g) / (${solute.mm.toFixed(2)} g/mol) = <strong>${f(n, 6)} mol</strong>.</p>` },
          { title: "Use the solvent mass in kilograms", body: `<p>mass(solvent) = ${solventG.toFixed(1)} g × (1 kg / 1000 g) = <strong>${f(kg)} kg</strong>.</p><p>The ${mass.toFixed(2)} g of solute is excluded from this denominator.</p>` },
          { title: "Divide moles by kilograms of solvent", body: `<p>m = n / mass(solvent) = (${f(n, 6)} mol) / (${f(kg)} kg) ≈ <strong>${H.sigfig(answer, 3)} mol/kg</strong> (3 significant figures).</p><p>Use the unrounded mole value. Molality is not moles per kilogram of the entire solution.</p>` }
        ],
        hint: "Water is the solvent. Convert only its mass to kilograms; convert solute grams to moles separately.",
        model: "<p>Complete dissolution is assumed. This is solute formula-unit molality; an ion-count factor belongs in a colligative-property calculation, not in this definition.</p>",
        data: { variant: "numeric", mass, mm: solute.mm, kg, answer }
      };
    }
  };

  const geometryOptions = [
    option("bent", "Bent"), option("linear", "Linear"),
    option("trigonal-pyramidal", "Trigonal pyramidal"),
    option("trigonal-planar", "Trigonal planar"), option("tetrahedral", "Tetrahedral")
  ];
  const molecules = [
    {
      id: "H2O", formula: "H₂O", name: "water", valence: 8, count: "2(1) + 6 = 8",
      bonds: 2, lonePairs: 2, loneAtom: "O", structure: "H—O(··)(··)—H",
      connectivity: "Two O–H single bonds and two lone pairs on oxygen. Each H has one bond.",
      geometry: "bent", geometryLabel: "Molecular geometry around O", geometryText: "O has four electron domains: two bonds and two lone pairs. Its electron-domain geometry is tetrahedral; the molecular geometry is bent.",
      polarBond: "O–H", negative: "O", positive: "H", polarityText: "Each O–H bond is polarized toward O: O δ− and H δ+. Water is bent, so the bond dipoles do not cancel; the molecule is polar.",
      electronegativity: "O > H", formalAtom: "O", formalValence: 6, formalNonbonding: 4, formalBonds: 2
    },
    {
      id: "HCl", formula: "HCl", name: "hydrogen chloride", valence: 8, count: "1 + 7 = 8",
      bonds: 1, lonePairs: 3, loneAtom: "Cl", structure: "H—Cl(··)(··)(··)",
      connectivity: "One H–Cl single bond and three lone pairs on chlorine.",
      geometry: "linear", geometryLabel: "Shape of the two-atom molecule", geometryText: "A two-atom molecule is linear. There is no central atom for a three-atom bond angle.",
      polarBond: "H–Cl", negative: "Cl", positive: "H", polarityText: "The H–Cl bond is polarized toward Cl: H δ+ and Cl δ−. An isolated HCl molecule is polar.",
      electronegativity: "Cl > H", formalAtom: "Cl", formalValence: 7, formalNonbonding: 6, formalBonds: 1
    },
    {
      id: "NH3", formula: "NH₃", name: "ammonia", valence: 8, count: "5 + 3(1) = 8",
      bonds: 3, lonePairs: 1, loneAtom: "N", structure: "    H\n    |\nH — N(··) — H",
      connectivity: "Three N–H single bonds and one lone pair on nitrogen.",
      geometry: "trigonal-pyramidal", geometryLabel: "Molecular geometry around N", geometryText: "N has four electron domains: three bonds and one lone pair. Its electron-domain geometry is tetrahedral; the molecular geometry is trigonal pyramidal.",
      polarBond: "N–H", negative: "N", positive: "H", polarityText: "Each N–H bond is polarized toward N: N δ− and H δ+. The trigonal pyramidal molecule has a net dipole and is polar.",
      electronegativity: "N > H", formalAtom: "N", formalValence: 5, formalNonbonding: 2, formalBonds: 3
    },
    {
      id: "HBr", formula: "HBr", name: "hydrogen bromide", valence: 8, count: "1 + 7 = 8",
      bonds: 1, lonePairs: 3, loneAtom: "Br", structure: "H—Br(··)(··)(··)",
      connectivity: "One H–Br single bond and three lone pairs on bromine.",
      geometry: "linear", geometryLabel: "Shape of the two-atom molecule", geometryText: "A two-atom molecule is linear. There is no central atom for a three-atom bond angle.",
      polarBond: "H–Br", negative: "Br", positive: "H", polarityText: "The H–Br bond is polarized toward Br: H δ+ and Br δ−. An isolated HBr molecule is polar.",
      electronegativity: "Br > H", formalAtom: "Br", formalValence: 7, formalNonbonding: 6, formalBonds: 1
    },
    {
      id: "CH3CH2MgBr", formula: "CH₃CH₂MgBr", name: "ethylmagnesium bromide", valence: 22, count: "2(4) + 5(1) + 2 + 7 = 22",
      bonds: 8, lonePairs: 3, loneAtom: "Br", structure: "    H   H\n    |   |\nH — C — C — Mg — Br(··)(··)(··)\n    |   |\n    H   H",
      connectivity: "Use five C–H, one C–C, one C–Mg, and one Mg–Br single bond. Put three lone pairs on Br. No other lone pairs occur in this specified covalent bookkeeping model.",
      geometry: "tetrahedral", geometryLabel: "Local geometry around the carbon bonded to Mg, in the specified model", geometryText: "The CH₂ carbon has four single bonds and no lone pairs, giving a local tetrahedral arrangement in this model. No prediction of the actual Mg coordination geometry is made.",
      polarBond: "C–Mg", negative: "C", positive: "Mg", polarityText: "The C–Mg bond is polarized toward carbon: C δ− and Mg δ+. The Mg–Br bond is polarized toward Br. These bond-polarity labels are distinct from integer formal charges.",
      electronegativity: "C > Mg (and Br > Mg)", formalAtom: "Mg", formalValence: 2, formalNonbonding: 0, formalBonds: 2,
      caveat: "For this exercise only, use the simplified covalent monomer skeleton CH₃–CH₂–Mg–Br. It is an electron-bookkeeping model, with an incomplete octet at Mg. Real Grignard solutions contain solvent-coordinated species and equilibria; this drawing does not specify their structures or proportions."
    },
    {
      id: "CH3OH", formula: "CH₃OH", name: "methanol", valence: 14, count: "4 + 4(1) + 6 = 14",
      bonds: 5, lonePairs: 2, loneAtom: "O", structure: "    H\n    |\nH — C — O(··)(··) — H\n    |\n    H",
      connectivity: "Use three C–H, one C–O, and one O–H single bond. Put two lone pairs on oxygen; carbon and hydrogen have no lone pairs.",
      geometry: "bent", geometryLabel: "Local molecular geometry around O", geometryText: "O has two bonds and two lone pairs: a bent local molecular geometry, with tetrahedral electron-domain geometry. Carbon has four bonds and a local tetrahedral geometry.",
      polarBond: "C–O", negative: "O", positive: "C", polarityText: "The C–O bond is polarized toward O: C δ+ and O δ− for that bond. The O–H bond is also polarized toward O, with its H δ+. Methanol is polar. C–H bonds are usually treated as nearly nonpolar in this introductory model.",
      electronegativity: "O > C and O > H", formalAtom: "O", formalValence: 6, formalNonbonding: 4, formalBonds: 2
    }
  ];

  G.lewis = {
    id: "lewis", group: "General questions", title: "Lewis structures & partial charges", slides: [49],
    description: "Practice the six review molecules, with electron counts, lone pairs, partial charges, and geometry.",
    generate(rng) {
      const molecule = rng.pick(molecules);
      const askFormal = rng.int(0, 1) === 1;
      const fields = [
        { ...H.num("valence", "Total valence electrons in the neutral formula", molecule.valence, "electrons"), tolerance: 0 },
        { ...H.num("lone-pairs", `Number of lone pairs on ${molecule.loneAtom}`, molecule.lonePairs, "pairs"), tolerance: 0 },
        choose("negative", `Which end of the ${molecule.polarBond} bond is δ−?`, molecule.negative, [
          option(molecule.negative, molecule.negative), option(molecule.positive, molecule.positive),
          option("neither", "Neither; the bond has no polarity")
        ], rng)
      ];
      if (askFormal) {
        fields.push(choose("formal-charge", `Formal charge on ${molecule.formalAtom} in the specified covalent Lewis structure`, "0", [
          option("0", "0"), option("-1", "−1"), option("+1", "+1"), option("+2", "+2")
        ], rng));
      } else {
        fields.push(choose("geometry", molecule.geometryLabel, molecule.geometry, geometryOptions, rng));
      }
      const formalText = `<p>Formal charge on ${molecule.formalAtom} = valence electrons − nonbonding electrons − ½(bonding electrons) = ${molecule.formalValence} − ${molecule.formalNonbonding} − ${molecule.formalBonds} = <strong>0</strong>.</p><p>A formal charge assigns bonding electrons equally. A partial charge (δ+ or δ−) describes unequal sharing, so an atom can have formal charge 0 and still carry a partial charge.</p>`;
      return {
        prompt: `<p>Draw the Lewis structure of <strong>${molecule.formula} (${molecule.name})</strong>, showing every lone pair and the relevant partial charges. Use paper or the scratchpad, then answer the checks below.</p>${molecule.caveat ? `<p>${molecule.caveat}</p>` : "<p>Consider the neutral molecular formula, rather than its ions after reaction with water.</p>"}<p><strong>Only the structured answers are checked automatically.</strong> Compare your sketch with the revealed solution yourself.</p>`,
        fields: shuffle(fields, rng),
        steps: [
          { title: "Count the valence electrons", body: `<p>Use H = 1, C = 4, N = 5, O = 6, Mg = 2, and Cl/Br = 7 valence electrons per neutral atom, as needed.</p><p>${molecule.formula}: ${molecule.count} valence electrons.</p>` },
          { title: "Connect atoms and place the remaining electrons", body: `<p>${molecule.connectivity}</p><p>${molecule.bonds} single bonds use ${2 * molecule.bonds} electrons. The remaining ${molecule.valence - 2 * molecule.bonds} electrons form ${molecule.lonePairs} lone pair${molecule.lonePairs === 1 ? "" : "s"} on ${molecule.loneAtom}.</p><pre class="lewis-diagram" style="overflow-x:auto">${molecule.structure}</pre><p>Each (··) represents one lone pair. Each line represents one electron-pair bond. This connectivity drawing does not show three-dimensional bond angles.</p>` },
          { title: "Read the local geometry", body: `<p>${molecule.geometryText}</p>` },
          { title: "Assign partial charges from bond polarity", body: `<p>The relevant electronegativity order is ${molecule.electronegativity}. The more electronegative end attracts the bonding electrons more strongly.</p><p>${molecule.polarityText}</p>` },
          { title: "Distinguish formal and partial charges", body: formalText },
          { title: "Compare your drawing", body: `<p>Check every atom, every single bond, the ${molecule.lonePairs} lone pair${molecule.lonePairs === 1 ? "" : "s"} on ${molecule.loneAtom}, and the direction of each marked bond polarity. A correct set of structured answers does not verify the drawing.</p>` }
        ],
        hint: `Count each atom's valence electrons first. A single bond uses two electrons, and δ− belongs at the more electronegative end of a bond. ${molecule.caveat ? "Use only the specified monomer bookkeeping model." : "Hydrogen needs a duet, not an octet."}`,
        model: molecule.caveat ? `<p>${molecule.caveat}</p>` : "<p>Lewis structures are electron-bookkeeping models. Geometry is interpreted with introductory VSEPR; partial-charge labels indicate bond polarity, not numerical atomic charges.</p>",
        data: { variant: askFormal ? "formal-charge" : "geometry", species: molecule.id, valence: molecule.valence, bonds: molecule.bonds, lonePairs: molecule.lonePairs, loneAtom: molecule.loneAtom, valenceCount: molecule.count, geometry: molecule.geometry, negative: molecule.negative, formalAtom: molecule.formalAtom, formalCharge: 0 }
      };
    }
  };
}());
