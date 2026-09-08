# Dillon Baker — My PhD Journey

A custom editorial redesign of the original Jekyll blog. All 13 original articles and their public URL paths are preserved. Article content remains in `_posts/`.

Two September 2026 entries cover Dillon's Orbital Observatory teaching model and election as TigerAI president. The orbital article embeds the existing angular-shape viewer and the render of the separate printable design. The self-contained viewer includes keyboard rotation, zoom/reset buttons, and a WebGL failure message; its Three.js MIT notice is retained. The GitHub download destination can be added to the orbital article once the new repository is available.

## Local development

Use Ruby 3.3, then run:

```sh
bundle install
bundle exec jekyll serve
```

The local site is served at http://127.0.0.1:4000. The standard production build is `JEKYLL_ENV=production bundle exec jekyll build`; public output is written to `dist/`.

For the private Sites copy, build with `--config _config.yml,_config.sites.yml`. The default configuration retains https://dtabaker.github.io as the canonical origin for GitHub Pages.

## Writing

Add Markdown files to `_posts/` with the usual YYYY-MM-DD-slug.md name and front matter:

```yaml
---
layout: post
title: "The full article title"
short_title: "A shorter title for lists"
date: 2026-09-07
topic: Research
tags: [Research]
description: "A one-sentence summary."
---
```

The homepage shows the newest post automatically. The `feature` values `main`, `quantum`, and `coding` select the three highlighted articles. Keep one article per featured slot.

Existing permalinks are explicit to preserve every original address. New articles use the same date-based structure. The archive groups posts using `tags`.

## Design and preservation

- Original chemistry-inspired molecular art, responsive editorial layouts, and reduced-motion support.
- Semantic navigation, keyboard focus, a skip link, readable article typography, RSS, and sitemap.
- Shared contact and research link includes keep the homepage, article sidebars, and every footer consistent. Article headers show the author and date without reading-time estimates.
- Original article words, photos, photograph credit, student blog links, and course link preserved. Redundant opening titles moved into front matter; nonstandard paragraph elements repaired.
- Social-sharing buttons and unused theme social icons are excluded. Email and profile links point directly to Dillon’s contact details.
- No analytics or external tracking added. Google Fonts supplies the typefaces, with local system fallbacks.
- An optional read-only WebMCP journal index is feature-detected on the archive. Runtime WebMCP verification was unavailable in this environment; ordinary navigation is independent of it.

The live GitHub Pages repository was not changed by the private Sites publication.

