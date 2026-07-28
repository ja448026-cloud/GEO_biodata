# Maintainer Architecture Notes

The public skill surface is intentionally small:

- `geo-biodata`
- `geo-biodata-intake`
- `geo-biodata-bulk`
- `geo-biodata-enrichment`
- `geo-biodata-scrna-intake`
- `geo-biodata-figure`
- `geo-biodata-workflow` as a deprecated compatibility alias

Stable user-facing executors live under `core/R`. During the transition period, wrappers forward to the compatibility implementation under `skills/geo-biodata-workflow/scripts` so validated behavior is not duplicated.

Do not promote a handoff into a formal skill until it has a deterministic executor, validation fixture, output contract, and route maturity entry.
