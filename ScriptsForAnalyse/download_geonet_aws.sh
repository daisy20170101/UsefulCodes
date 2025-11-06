#!/bin/bash
#
# Download GeoNet seismic data from AWS S3
#
# GeoNet provides open access to New Zealand seismic waveform data via AWS S3
# Data structure: s3://geonet-open-data/waveforms/miniseed/YEAR/YEAR.DOY/
#
# Requirements:
#   - AWS CLI (apt-get install awscli or brew install awscli)
#
# Usage:
#   ./download_geonet_aws.sh
#   ./download_geonet_aws.sh 2025-11-06 WEL,BFZ,PUZ
#

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Default parameters
DATE=${1:-"2025-11-06"}
STATIONS=${2:-"WEL,BFZ,PUZ,MXZ,THZ,KNZ"}
NETWORK="NZ"
LOCATION="10"
CHANNELS="HHZ,HHN,HHE"
OUTPUT_DIR="./geonet_data"

echo "======================================================================"
echo "GEONET AWS S3 DATA DOWNLOADER"
echo "======================================================================"
echo "Date: $DATE"
echo "Stations: $STATIONS"
echo "Channels: $CHANNELS"
echo "Output: $OUTPUT_DIR"
echo "======================================================================"
echo ""

# Check if AWS CLI is installed
if ! command -v aws &> /dev/null; then
    echo -e "${RED}✗ AWS CLI not found${NC}"
    echo ""
    echo "Please install AWS CLI:"
    echo "  Ubuntu/Debian: sudo apt-get install awscli"
    echo "  macOS: brew install awscli"
    echo "  pip: pip install awscli"
    exit 1
fi

echo -e "${GREEN}✓ AWS CLI found: $(aws --version)${NC}"

# Convert date to year and day of year
YEAR=$(date -d "$DATE" +%Y)
DOY=$(date -d "$DATE" +%j)

echo "Year: $YEAR, Day of Year: $DOY"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR/$YEAR.$DOY"

# S3 base path
S3_BASE="s3://geonet-open-data/waveforms/miniseed/$YEAR/$YEAR.$DOY"

echo "S3 path: $S3_BASE"
echo ""

# Split stations and channels into arrays
IFS=',' read -ra STATION_ARRAY <<< "$STATIONS"
IFS=',' read -ra CHANNEL_ARRAY <<< "$CHANNELS"

# Download counters
TOTAL=0
SUCCESS=0
FAILED=0

echo "Starting download..."
echo "======================================================================"

# Download data for each station and channel
for STATION in "${STATION_ARRAY[@]}"; do
    echo ""
    echo "Station: $STATION"

    for CHANNEL in "${CHANNEL_ARRAY[@]}"; do
        TOTAL=$((TOTAL + 1))

        # Construct filename: NETWORK.STATION.LOCATION.CHANNEL.D.YEAR.DOY
        FILENAME="${NETWORK}.${STATION}.${LOCATION}.${CHANNEL}.D.${YEAR}.${DOY}"
        S3_PATH="${S3_BASE}/${FILENAME}"
        LOCAL_PATH="${OUTPUT_DIR}/${YEAR}.${DOY}/${FILENAME}"

        # Check if file already exists
        if [ -f "$LOCAL_PATH" ]; then
            SIZE=$(du -h "$LOCAL_PATH" | cut -f1)
            echo -e "  ${YELLOW}⊙${NC} $CHANNEL - already exists ($SIZE)"
            SUCCESS=$((SUCCESS + 1))
            continue
        fi

        # Download file
        if aws s3 cp --no-sign-request "$S3_PATH" "$LOCAL_PATH" 2>/dev/null; then
            SIZE=$(du -h "$LOCAL_PATH" | cut -f1)
            echo -e "  ${GREEN}✓${NC} $CHANNEL - downloaded ($SIZE)"
            SUCCESS=$((SUCCESS + 1))
        else
            echo -e "  ${RED}✗${NC} $CHANNEL - not available"
            FAILED=$((FAILED + 1))
        fi
    done
done

echo ""
echo "======================================================================"
echo "DOWNLOAD SUMMARY"
echo "======================================================================"
echo "Total attempts: $TOTAL"
echo -e "${GREEN}Successful: $SUCCESS${NC}"
echo -e "${RED}Failed: $FAILED${NC}"
echo "Success rate: $(awk "BEGIN {printf \"%.1f\", ($SUCCESS/$TOTAL)*100}")%"
echo "Output directory: $OUTPUT_DIR"
echo "======================================================================"

# List downloaded files
echo ""
echo "Downloaded files:"
find "$OUTPUT_DIR" -type f -name "*.D.$YEAR.$DOY" -exec ls -lh {} \; | awk '{print "  " $9 " (" $5 ")"}'

echo ""
echo "To view with ObsPy:"
echo "  from obspy import read"
echo "  st = read('$OUTPUT_DIR/$YEAR.$DOY/${NETWORK}.${STATION_ARRAY[0]}.${LOCATION}.${CHANNEL_ARRAY[0]}.D.$YEAR.$DOY')"
echo "  st.plot()"
