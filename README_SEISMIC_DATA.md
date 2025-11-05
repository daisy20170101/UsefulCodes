# Seismic Moment Release Data for Major Earthquakes

This directory contains scripts and data for seismic moment release rates of major earthquakes (M 7.5-9.0).

## Data Files

### 1. `major_earthquakes_moment_data.csv`
Complete dataset with all calculated parameters:
- `time`: Event date and time
- `latitude`, `longitude`: Epicenter coordinates
- `depth_km`: Focal depth in kilometers
- `magnitude`: Moment magnitude (Mw)
- `seismic_moment_Nm`: Seismic moment in Newton-meters (N·m)
- `cumulative_moment_Nm`: Cumulative moment release up to that event
- `moment_release_rate_Nm_per_year`: Moment release rate (1-year moving window)
- `location`: Event location description
- `years_since_start`: Time in years from first event

### 2. `major_earthquakes_simple.csv`
Simplified version with essential parameters only

### 3. `gutenberg_richter_distribution.csv`
Frequency-magnitude distribution:
- `magnitude_threshold`: Magnitude threshold
- `cumulative_count`: Number of events ≥ threshold
- `annual_rate`: Annual rate of events ≥ threshold

## Scripts

### 1. `create_seismic_moment_database.py`
Creates the seismic moment database from curated historical earthquake records.

**Usage:**
```bash
python3 create_seismic_moment_database.py
```

**Features:**
- Curated database of 77+ major earthquakes (M 7.5-9.0) from 1868-2024
- Calculates seismic moment from magnitude using: Mw = (2/3)×log₁₀(M₀) - 10.7
- Computes cumulative moment release
- Calculates moment release rates (1-year moving window)
- Generates Gutenberg-Richter frequency-magnitude distribution
- Comprehensive summary statistics

### 2. `download_seismic_moment_data.py`
Script to download data from USGS earthquake catalog (requires internet access)

### 3. `download_global_cmt_data.py`
Alternative script to download from Global CMT catalog

## Key Statistics (Current Dataset)

- **Total earthquakes**: 77 events
- **Time span**: 1868-2024 (155+ years)
- **Magnitude range**: M 6.2 - 9.5
- **Total moment released**: 6.095×10²³ N·m
- **Average moment rate**: 3.923×10²¹ N·m/year

### Largest Earthquakes

1. **M 9.5** - 1960-05-22 - Valdivia, Chile (1.995×10²³ N·m)
2. **M 9.2** - 1964-03-28 - Prince William Sound, Alaska (7.079×10²² N·m)
3. **M 9.1** - 2004-12-26 - Sumatra-Andaman Islands (5.012×10²² N·m)
4. **M 9.1** - 2011-03-11 - Tohoku, Japan (5.012×10²² N·m)
5. **M 9.0** - 1952-11-04 - Kamchatka, Russia (3.548×10²² N·m)

### Magnitude Distribution

- M 9.0-10.0: 5 events (6.5%)
- M 8.5-9.0: 14 events (18.2%)
- M 8.0-8.5: 18 events (23.4%)
- M 7.5-8.0: 34 events (44.2%)
- M < 7.5: 6 events (7.8%)

## Seismic Moment Calculation

The seismic moment M₀ is calculated from moment magnitude Mw using:

```
Mw = (2/3) × log₁₀(M₀) - 10.7
```

where M₀ is in dyne·cm. We convert to N·m using:
```
M₀ (N·m) = M₀ (dyne·cm) × 10⁻⁷
```

This gives:
```
M₀ (N·m) = 10^(1.5×Mw + 16.05) × 10⁻⁷
```

## Data Sources

The curated earthquake database is compiled from:
- USGS Earthquake Catalog
- Global CMT Catalog
- Published scientific literature

## Notes

- The database focuses on well-documented major earthquakes
- Historical earthquakes (pre-1900) may have magnitude uncertainties
- Seismic moment values are calculated from reported magnitudes
- For the most up-to-date catalog, use the download scripts with internet access

## Requirements

```bash
pip install pandas numpy requests
```

## References

1. USGS Earthquake Catalog: https://earthquake.usgs.gov/
2. Global CMT Catalog: https://www.globalcmt.org/
3. Hanks & Kanamori (1979). A moment magnitude scale. JGR, 84(B5), 2348-2350.
4. Dziewonski et al. (1981). Determination of earthquake source parameters. JGR, 86(B4), 2825-2852.

## Author

Generated for seismic moment analysis and research purposes.

Last updated: 2025-11-05
