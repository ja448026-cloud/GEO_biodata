# Maintainer Architecture Notes

The public skill surface is intentionally small:

- `geo-biodata`
- `geo-biodata-intake`
- `geo-biodata-bulk`
- `geo-biodata-enrichment`
- `geo-biodata-scrna-intake`
- `geo-biodata-figure`
- `geo-biodata-workflow` as a deprecated compatibility alias

Stable user-facing executors live under `core/R`. The legacy workflow skill is retained for compatibility with old prompts; it is not the primary execution surface.

Do not promote a handoff into a formal skill until it has a deterministic executor, validation fixture, output contract, and route maturity entry.
