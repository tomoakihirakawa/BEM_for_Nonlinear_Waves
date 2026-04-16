# Public Literature Benchmark Notes

The actual public benchmark case directories are now located directly under
`gui/examples/`, not under this directory. This directory only keeps the shared
case template and checklist.

The intended use is:

- extract cable, wind, rain, damper, and response data from published sources;
- record source quality and missing quantities explicitly;
- convert usable cases into solver or GUI input files;
- keep confidential project notes out of the benchmark dataset.

## Candidate Cases

| Directory | Case | Primary use |
|---|---|---|
| `../dongting_lake_2007/` | Dongting Lake Bridge RWIV field measurements | Main full-scale rain-wind benchmark |
| `../stavanger_city_2021/` | Stavanger City Bridge wet/dry cable observations | Open-access full-scale RWIV comparison |
| `../east_china_sea_viv_2019/` | Long-span bridge VIV field investigation | High-mode and multimode VIV benchmark |
| `../nrc_rwiv_2023/` | NRC Canada large-scale RWIV test | Controlled parametric RWIV benchmark |
| `../fhwa_dry_inclined_2007/` | FHWA/NRC dry inclined cable tests | Damping and Scruton-number design reference |

Existing sibling directories such as `../fred_hartman_1995`,
`../stonecutters_2009`, and `../sutong_2008` contain older case notes. They
can either stay as standalone examples or be summarized here later.

## Per-Case Checklist

Each case should record:

- source references and DOI/URL;
- bridge or test-rig geometry;
- cable length, diameter, mass per unit length, inclination, and tension if
  available;
- measured natural frequencies and damping ratios;
- damper type, coefficient, location, and support condition if available;
- wind speed, direction, rain rate, turbulence, and reduced velocity;
- response amplitude, mode number, time segment, and whether the value is
  peak, peak-to-peak, RMS, or modal amplitude;
- confidence level: direct table value, figure digitization, inferred value,
  or missing.

## File Convention

Suggested files inside each case directory:

- `README.md`: human-readable case summary.
- `sources.md`: bibliographic notes, links, access status.
- `extracted_data.csv`: tabular values taken from tables or digitized figures.
- `assumptions.md`: inferred values and unit conversions.
- `settings.json`: solver/GUI input when a reproducible case is ready.
