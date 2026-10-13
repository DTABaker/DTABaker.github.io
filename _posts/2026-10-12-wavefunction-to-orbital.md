---
layout: post
title: "From wavefunction to orbital: explore the mathematics of shape"
short_title: "From wavefunction to orbital"
date: 2026-10-12 17:00:00 -0500
permalink: /wavefunction-to-orbital/
topic: "Chemistry & teaching"
tags: ["Teaching & outreach"]
description: "Follow one point through four linked views of a hydrogenic orbital: radial function, angular function, wavefunction cross-section, and three-dimensional probability surface."
---

How does a wavefunction become the orbital shapes we draw in chemistry? This companion to the [Orbital Observatory]({{ '/2026/09/07/orbital-observatory.html' | relative_url }}) connects the mathematics to the picture. Move one shared probe and watch its distance, direction, wavefunction amplitude, and probability density update together.

The model includes all **16 Observatory orbitals**—1s, three 2p, five 3d, and seven 4f—plus **2s, 3s, and 3p_z** to explore radial nodes. These 19 states form a selected teaching set, rather than every orbital through n = 4.

<p><a href="{{ '/assets/wavefunction-to-orbital.html' | relative_url }}" target="_blank" rel="noopener">Open From Wavefunction to Orbital in a full window ↗</a></p>

<iframe
  id="wavefunction-orbital-embed"
  src="{{ '/assets/wavefunction-to-orbital.html' | relative_url }}"
  title="From Wavefunction to Orbital: four linked views of quantum mathematics"
  width="100%" height="1500"
  style="display:block; width:100%; border:1px solid #d8dfdf; border-radius:10px;"
  loading="lazy">
</iframe>

<noscript><p>The interactive model needs JavaScript, and its 3D views require WebGL. The explanation below remains available without JavaScript.</p></noscript>

The model scrolls within the frame. For more room, especially on a phone, use the full-window link above.

## Follow a point through four views

1. **Distance:** the radial function R(r) shows how amplitude changes with distance from the nucleus. The second graph shows the shell probability density P(r) = r²|R(r)|², including the radial volume weighting.
2. **Direction:** the angular function Y(θ, φ) colors a unit sphere. Rotate it or change the angles to explore signs, orientations, and angular nodes.
3. **Product:** the cross-section shows ψ = RY or the probability density |ψ|². A node in either factor makes the wavefunction zero.
4. **Shape:** the three-dimensional surface shows a region containing approximately 90% of the orbital's probability before cutting. Move the cutting plane to inspect the interior; its position stays linked to the probe and cross-section.

Blue and coral indicate opposite signs of the real spatial wavefunction, not opposite electrical charges. Squaring removes the sign but preserves the nodes. The squared amplitude |ψ|² at a point is a probability density; a probability comes from integrating over a region.

## Try these experiments

- Select **1s** and compare the radial function with the shell probability graph. Why is the most probable radial distance away from the nucleus even though the density is greatest at the nucleus?
- Use **2s** to find a radial node, then **2p_z** to find an angular node. Try **3p_z** to see both in one state.
- Change the azimuth of **2p_x** or **3d_xy**. Follow a sign change in Y into ψ, then switch to |ψ|².
- Choose **4f_xyz** and inspect its three nodal planes. A central coordinate slice can lie entirely on a node; move the cut to reveal the surrounding lobes.

The guided-example buttons set up these investigations. The equation panel displays the selected normalized functions and evaluates them at the shared probe.

## What the model represents

These are normalized, nonrelativistic **hydrogenic wavefunctions with one electron and Z = 1**, treating the nucleus as fixed. Electron–electron interactions, spin dynamics, and relativistic effects are omitted, so the model does not give quantitative orbital sizes or energies for neutral many-electron atoms.

Lengths are measured in Bohr radii, a₀. The physical plotting window is labeled and changes with the selected orbital; this viewer is not a common-screen-scale comparison of orbital sizes. The angular sphere has unit radius and represents directions, not atomic size.

Most of the familiar real p, d, and f shapes combine complex states with +m and −m. They remain eigenstates of energy and total orbital angular momentum squared, L², but generally do not have a single L_z eigenvalue. The model identifies these real combinations and displays their full angular formulas.

The approximately 90% surfaces are numerical isodensity boundaries, not hard edges or electron paths. Cutting away part of a surface removes part of the displayed probability region. The radial panel's cumulative probability instead refers to a sphere centered on the nucleus; it describes a different region.

The model runs in your browser without an account. All runtime code and geometry are bundled in the standalone page. Its mathematical functions, mesh probabilities, linked controls, and desktop/mobile layouts were reviewed before publication; browser checks covered Chromium rather than every browser or device.

For the underlying mathematics, see [Oregon State University's hydrogen radial functions and spherical harmonics](https://books.physics.oregonstate.edu/GMM/hydroform.html), Richard Fitzpatrick's [Hydrogen Atom](https://farside.ph.utexas.edu/teaching/qmech/lectures/node82.html), and [Spherical Harmonics](https://farside.ph.utexas.edu/teaching/qmech/lectures/node76.html). The bundled Three.js renderer retains its [MIT license]({{ '/assets/THREE-LICENSE.txt' | relative_url }}).
