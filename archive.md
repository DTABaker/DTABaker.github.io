---
layout: default
title: The Journal
description: "Every entry from Dillon Baker’s PhD journey, organized by research, teaching and outreach, and academic life."
permalink: /archive.html
---
<section class="archive-header section-wrap">
<span class="eyebrow section-kicker">THE COMPLETE JOURNAL / {{ site.posts.size }} ENTRIES</span>
<h1>The notebook.</h1>
<p>Research, teaching, and everything learned along the way.</p>
<nav class="archive-topics" aria-label="Browse journal topics">{% for tag in site.tags %}<a href="#{{ tag[0] | slugify }}">{{ tag[0] }} <span>{{ tag[1].size }}</span></a>{% endfor %}</nav>
</section>
<div class="archive-content section-wrap">{% for tag in site.tags %}<section class="archive-group" id="{{ tag[0] | slugify }}" aria-labelledby="heading-{{ tag[0] | slugify }}"><div class="archive-group-title"><span class="eyebrow">0{{ forloop.index }} / {{ tag[1].size }} ENTRIES</span><h2 id="heading-{{ tag[0] | slugify }}">{{ tag[0] }}</h2></div><div class="archive-entries">{% for post in tag[1] %}<a class="archive-entry" href="{{ post.url | relative_url }}"><div class="archive-entry-meta"><time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: '%d %b %Y' }}</time><span>{{ post.topic }}</span></div><h3>{{ post.short_title | default: post.title }}</h3><p>{{ post.description }}</p><span class="row-arrow" aria-hidden="true">↗</span></a>{% endfor %}</div></section>{% endfor %}</div>
<script type="application/json" id="journal-index">[{% for post in site.posts %}{"title":{{ post.title | jsonify }},"date":{{ post.date | date: '%Y-%m-%d' | jsonify }},"topic":{{ post.tags[0] | jsonify }},"url":{{ post.url | relative_url | jsonify }}}{% unless forloop.last %},{% endunless %}{% endfor %}]</script>
