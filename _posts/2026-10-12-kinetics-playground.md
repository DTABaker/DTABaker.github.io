---
layout: post
title: "Kinetics playground: explore reaction rates and rate laws"
short_title: "A playground for reaction kinetics"
date: 2026-10-12 12:00:00 -0500
permalink: /kinetics/
topic: "Chemistry & teaching"
tags: ["Teaching & outreach"]
description: "An interactive lab for exploring reaction order, concentration, half-life, temperature, and catalysts, with an initial-rate detective challenge."
---

What happens to a reaction rate when you double the starting concentration? Does doubling the temperature double the rate? Explore those questions by changing one variable at a time, watching the reaction progress, and comparing the graphs.

<p><a href="{{ '/assets/kinetics-playground.html' | relative_url }}" target="_blank" rel="noopener">Open the kinetics playground in a full page</a></p>

<iframe
  id="kinetics-playground-embed"
  src="{{ '/assets/kinetics-playground.html' | relative_url }}"
  title="Interactive kinetics playground: reaction rates and rate laws"
  width="100%"
  height="1150"
  style="display:block; width:100%; border:0;"
  loading="lazy">
</iframe>

<script>
(() => {
  const frame = document.getElementById('kinetics-playground-embed');
  if (!frame) return;
  const origin = new URL(frame.src, document.baseURI).origin;
  window.addEventListener('message', event => {
    if (event.origin !== origin || event.source !== frame.contentWindow) return;
    const data = event.data;
    if (data?.type !== 'kinetics-playground:resize') return;
    if (!Number.isFinite(data.height) || data.height < 200 || data.height > 10000) return;
    frame.style.height = Math.ceil(data.height) + 'px';
  });
})();
</script>

<noscript><p>This playground needs JavaScript. The learning prompts and model description below remain available without it.</p></noscript>

## Try these experiments

1. **Predict, then test.** Select zero, first, or second order. Double the starting concentration while keeping temperature and the rate constant fixed. How does the initial rate change?
2. **Follow a half-life.** Scrub through time and compare how long it takes to consume half the reactant. Change the starting concentration and repeat. The time window rescales to four initial half-lives for each setting.
3. **Find the straight line.** Switch to the linearized graphs and compare concentration, logarithmic, and reciprocal plots. Which transformation is linear for each order?
4. **Solve a mystery.** Open the rate-law detective, collect initial rates at different concentrations, and infer the reaction order before revealing the answer.

Temperature and catalyst controls offer a second investigation: predict which way the rate constant and half-life will move, then compare your prediction with the model.

## What the model represents

The lab uses ideal, irreversible **A → P** kinetics with 1:1 stoichiometry in a closed, well-mixed system at constant volume and temperature during each run. Reaction order is an empirical property of the rate law; it does not follow automatically from the reaction equation. The particle animation illustrates progress, rather than simulating molecular collisions.

The temperature control uses the Arrhenius equation. The hypothetical catalyst lowers the activation energy from 50 to 45 kJ/mol while keeping the prefactor and reaction order fixed. The detective's measurements are exact synthetic data, without experimental noise. Expand the model's assumptions for the equations, units, and further limitations.

The model runs in your browser and remembers settings locally when browser storage is available. It needs no account or server connection after loading.

For the underlying chemistry, see OpenStax *Chemistry 2e*: [Rate Laws](https://openstax.org/books/chemistry-2e/pages/12-3-rate-laws), [Integrated Rate Laws](https://openstax.org/books/chemistry-2e/pages/12-4-integrated-rate-laws), [Collision Theory](https://openstax.org/books/chemistry-2e/pages/12-5-collision-theory), and [Catalysis](https://openstax.org/books/chemistry-2e/pages/12-7-catalysis). [Third-party software notices]({{ '/assets/kinetics-THIRD-PARTY-LICENSES.txt' | relative_url }}) accompany the bundled model.
