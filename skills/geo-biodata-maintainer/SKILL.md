---
name: geo-biodata-maintainer
description: Maintain the GEO_biodata repository. Use when editing skill structure, README/docs, dependency profiles, validation fixtures, smoke tests, GitHub Actions, release checks, or repository hygiene. Do not run heavy scientific validation unless explicitly requested.
---

# GEO_biodata Maintainer

Use this for repository changes, not ordinary analysis.

## Defaults

- Keep public-facing entry points short.
- Prefer narrow skills over one large workflow skill.
- Keep existing script paths compatible unless a migration is explicitly planned.
- Run quick local validation before pushing.
- Keep heavy Bioconductor/scRNA validation manual or path-scoped.

## Validation Tiers

- Quick: parse scripts, validate manifests, download guards, skill frontmatter.
- Bulk: bulk/microarray driver fixtures and negative scale tests.
- scRNA: object-intake fixture and dangerous dense-conversion checks.
- Real data: manual workflow dispatch with accession-specific evidence.

