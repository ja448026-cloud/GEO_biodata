# Agent operating contract

> FREEZE NOTICE (2026-07-29): No new gates, contracts, YAMLs, schemas,
> decision rules, or smoke-check assertions until Slim-Down Phases 1-5
> are complete. Only deletions/consolidations are allowed.
> smoke checks size cap: 34,459 bytes.

This repository is a reusable GEO biodata workflow for coding agents (Claude Code, Codex, etc.). Read `skills/geo-biodata-workflow/SKILL.md` before working on any GEO accession.

## For Claude Code

This project works as a Claude Code skill. To use it:

1. Copy or symlink `skills/geo-biodata-workflow/` into `.claude/skills/` in your project or home directory.
2. Invoke via `/geo-biodata-workflow <GSE accession>` or by asking Claude to inspect a GSE.

The skill is self-contained: it brings its own R scripts and reference documents. No additional installation beyond the R packages listed in `check_environment.R`.

## Repository rules

Keep the repository generic. Do not introduce private data, local-project assumptions, disease-specific logic, local absolute paths, generated data, or unpublished results.

For each accession:

1. discover resources before downloading;
2. inspect metadata and author supplements before choosing an analysis route;
3. preserve raw files and author labels;
4. record every mapping or inference;
5. generate diagnostic figures from saved tables or objects;
6. stop explicitly when inputs, metadata, or compute are insufficient.

Do not install packages silently. Report missing dependencies and provide reproducible installation instructions. Do not treat cells as biological replicates, infer missing experimental groups, or scrape paywalled resources.
