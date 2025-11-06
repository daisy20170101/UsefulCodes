# GeoNet AWS S3 Data Download Tools

Tools for downloading seismic waveform data from GeoNet's AWS S3 bucket.

## Data Source

GeoNet provides open access to New Zealand seismic data via AWS S3:
- **Bucket**: `s3://geonet-open-data/waveforms/miniseed/`
- **Format**: MiniSEED (standard seismology format)
- **Update frequency**: Near real-time
- **Access**: Public, no AWS account required

## Data Structure

```
s3://geonet-open-data/waveforms/miniseed/
└── YEAR/
    └── YEAR.DOY/                    (DOY = Day of Year, e.g., 310)
        └── NETWORK.STATION.LOCATION.CHANNEL.D.YEAR.DOY
```

Example:
```
s3://geonet-open-data/waveforms/miniseed/2025/2025.310/NZ.WEL.10.HHZ.D.2025.310
```

## Requirements

### AWS CLI
```bash
# Ubuntu/Debian
sudo apt-get install awscli

# macOS
brew install awscli

# pip
pip install awscli
```

### Python packages (for Python script)
```bash
pip install obspy boto3
```

## Available Tools

### 1. Bash Script (Simple and Fast)

**File**: `download_geonet_aws.sh`

Quick downloads using AWS CLI directly.

**Usage:**
```bash
# Default: Download data for 2025-11-06
./download_geonet_aws.sh

# Specify date
./download_geonet_aws.sh 2025-11-06

# Specify date and stations
./download_geonet_aws.sh 2025-11-06 WEL,BFZ,PUZ,MXZ

# Custom everything
DATE=2025-11-05 STATIONS=WEL,BFZ CHANNELS=HHZ,HHN ./download_geonet_aws.sh
```

**Features:**
- Fast parallel-like downloads
- Color-coded output
- Automatic skip if file exists
- Download statistics

### 2. Python Script (Advanced Features)

**File**: `download_geonet_aws_data.py`

More features including file merging and ObsPy integration.

**Usage:**
```bash
# List S3 bucket structure
python download_geonet_aws_data.py --list

# List data for specific date
python download_geonet_aws_data.py --list --date 2025-11-06

# Download data for default stations
python download_geonet_aws_data.py --date 2025-11-06

# Download for specific stations
python download_geonet_aws_data.py --date 2025-11-06 --stations WEL,BFZ,PUZ

# Download with time range (multi-day)
python download_geonet_aws_data.py --date 2025-11-06 --hours-before 12 --hours-after 36

# Download and merge files by station
python download_geonet_aws_data.py --date 2025-11-06 --merge

# Custom output directory
python download_geonet_aws_data.py --date 2025-11-06 --output ./my_data
```

**Features:**
- Multi-day download support
- Automatic file merging (requires ObsPy)
- Station metadata support
- Progress tracking
- Error handling

## Station Codes

### Major GeoNet Stations

| Code | Name | Location | Description |
|------|------|----------|-------------|
| WEL | Wellington | North Island | Capital city |
| BFZ | Birch Farm | North Island | Near Wellington |
| PUZ | Putaroa | North Island | Central |
| MXZ | Mangaotuku | North Island | Central |
| THZ | Top House | South Island | Northern |
| KNZ | Kaikoura | South Island | Northeast coast |
| ODZ | Otago University | South Island | Dunedin |
| DCZ | Denniston | South Island | West coast |
| WVZ | White Island | North Island | Volcanic |
| RPZ | Raoul Island | Kermadec Islands | Remote |

### Network Codes
- **NZ**: GeoNet permanent network (most common)
- **NZ.**: GeoNet locations

### Channel Codes

| Code | Description | Sample Rate |
|------|-------------|-------------|
| HH* | High broadband | 100 Hz |
| BH* | Broadband | 20-80 Hz |
| EH* | Extremely short period | 100+ Hz |
| HN* | High-gain | 100 Hz |

Component suffixes:
- **Z**: Vertical
- **N**: North
- **E**: East
- **1**: Horizontal 1
- **2**: Horizontal 2

## Example Workflows

### 1. Download New Zealand earthquake data
```bash
# For 2025-11-06 earthquake
./download_geonet_aws.sh 2025-11-06 WEL,BFZ,PUZ,MXZ,THZ,KNZ

# Multi-day download with Python
python download_geonet_aws_data.py \
    --date 2025-11-06 \
    --hours-before 6 \
    --hours-after 48 \
    --stations WEL,BFZ,PUZ \
    --merge
```

### 2. View downloaded data with ObsPy
```python
from obspy import read

# Read single file
st = read('geonet_data/2025.310/NZ.WEL.10.HHZ.D.2025.310')
print(st)
st.plot()

# Read all files for a station
st = read('geonet_data/2025.310/NZ.WEL.*.D.2025.310')
st.plot()

# Process data
st.detrend('linear')
st.taper(max_percentage=0.05)
st.filter('bandpass', freqmin=0.5, freqmax=20.0)
st.plot()
```

### 3. Merge multiple days
```python
from obspy import read, Stream

# Read multiple days
st = Stream()
st += read('geonet_data/2025.309/NZ.WEL.10.HHZ.D.2025.309')
st += read('geonet_data/2025.310/NZ.WEL.10.HHZ.D.2025.310')
st += read('geonet_data/2025.311/NZ.WEL.10.HHZ.D.2025.311')

# Merge with gap handling
st.merge(method=1, fill_value='interpolate')

print(f"Total duration: {st[0].stats.endtime - st[0].stats.starttime} seconds")
st.plot()
```

## Output Structure

```
geonet_data/
├── 2025.309/
│   ├── NZ.WEL.10.HHZ.D.2025.309
│   ├── NZ.WEL.10.HHN.D.2025.309
│   ├── NZ.WEL.10.HHE.D.2025.309
│   ├── NZ.BFZ.10.HHZ.D.2025.309
│   └── ...
└── 2025.310/
    ├── NZ.WEL.10.HHZ.D.2025.310
    └── ...

geonet_merged/          (if using --merge)
├── NZ.WEL.10.HHZ.mseed
├── NZ.WEL.10.HHN.mseed
└── ...
```

## Troubleshooting

### AWS CLI not found
```bash
# Install AWS CLI
sudo apt-get install awscli
# or
brew install awscli
# or
pip install awscli
```

### Date conversion errors (bash script)
The bash script uses GNU date. On macOS, install GNU date:
```bash
brew install coreutils
# Then use: gdate instead of date in the script
```

### File not found in S3
- Check that the date is correct (format: YYYY-MM-DD)
- Verify station code is correct (case-sensitive)
- Some stations may not have data for all time periods
- Try listing the bucket: `aws s3 ls --no-sign-request s3://geonet-open-data/waveforms/miniseed/2025/2025.310/`

### Cannot merge files
Ensure ObsPy is installed:
```bash
pip install obspy
```

## Additional Resources

- **GeoNet website**: https://www.geonet.org.nz/
- **GeoNet data services**: https://www.geonet.org.nz/data
- **AWS Open Data Registry**: https://registry.opendata.aws/geonet/
- **ObsPy documentation**: https://docs.obspy.org/

## Notes

- Data is organized by day (UTC)
- Files are typically ~1-10 MB per station-channel-day
- Download only what you need to save time and bandwidth
- GeoNet data is provided by GNS Science and funded by EQC
- No AWS account or authentication required
