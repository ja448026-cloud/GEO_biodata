# Maintainer Architecture Notes

The public skill surface is intentionally small:

- `geo-biodata`
- `geo-biodata-intake`
- `geo-biodata-bulk`
- `geo-biodata-enrichment`
- `geo-biodata-scrna`
- `geo-biodata-figure`
- `geo-biodata-workflow` as a deprecated compatibility alias

Stable helper scripts live under `core/R`. The legacy workflow skill is retained for compatibility with old prompts; it should stay as thin shims over current helpers where possible.

Prefer short skill instructions over duplicated YAML rules. Add code only when it removes manual GEO busywork or prevents a common wrong analysis.
