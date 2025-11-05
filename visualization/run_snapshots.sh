#!/bin/bash
#
# Wrapper script for fault snapshot generation
# Automatically finds pvpython and runs the visualization scripts
#

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STATE_FILE="$SCRIPT_DIR/rakeNew.pvsm"
FALLBACK_STATE_FILE="/Users/DuoL/Documents/NSHM/Central/Paraviews/rakeNew.pvsm"
OUTPUT_DIR="./snapshots"
RESOLUTION="1920 1080"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Function to find pvpython
find_pvpython() {
    # Check if pvpython is in PATH
    if command -v pvpython &> /dev/null; then
        echo "pvpython"
        return 0
    fi

    # Check common installation paths
    local common_paths=(
        "/Applications/ParaView-*.app/Contents/bin/pvpython"
        "/usr/local/bin/pvpython"
        "/opt/paraview/bin/pvpython"
        "$HOME/ParaView-*/bin/pvpython"
    )

    for path in "${common_paths[@]}"; do
        if ls $path &> /dev/null; then
            echo "$path"
            return 0
        fi
    done

    return 1
}

# Function to show usage
show_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Generate fault slip rate snapshots using ParaView

Options:
    -h, --help              Show this help message
    -i, --input FILE        Input XDMF file
    -m, --mode MODE         Mode: single|batch|timesteps (default: single)
    -d, --data-dir DIR      Data directory (for batch mode)
    -s, --state FILE        ParaView state file (optional)
    -o, --output DIR        Output directory (default: $OUTPUT_DIR)
    -r, --resolution WxH    Resolution (default: $RESOLUTION)
    -p, --pattern PATTERN   File pattern for batch mode (default: *.xdmf)
    -f, --field FIELD       Field to visualize: SRs|SRd|all (default: all)

Modes:
    single      Generate snapshots from XDMF file or state file (default)
    batch       Process all files in data directory
    timesteps   Process all timesteps from state file

Examples:
    # Generate snapshot from XDMF file
    $0 -i fault_output.xdmf

    # Generate snapshot with custom state file
    $0 -i fault_output.xdmf -s rakeNew.pvsm

    # Generate single snapshot from state file only (legacy)
    $0 -m single -s /path/to/state.pvsm

    # Batch process data files
    $0 -m batch -d /path/to/data -p "fault_*.xdmf"

    # Process all timesteps
    $0 -m timesteps -s /path/to/state.pvsm

    # High resolution output
    $0 -i data.xdmf -r 3840x2160

EOF
}

# Parse command line arguments
MODE="single"
INPUT_FILE=""
DATA_DIR=""
PATTERN="*.xdmf"
FIELD="all"

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_usage
            exit 0
            ;;
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -m|--mode)
            MODE="$2"
            shift 2
            ;;
        -d|--data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -s|--state)
            STATE_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -r|--resolution)
            RESOLUTION="${2/x/ }"
            shift 2
            ;;
        -p|--pattern)
            PATTERN="$2"
            shift 2
            ;;
        -f|--field)
            FIELD="$2"
            shift 2
            ;;
        *)
            print_error "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Main execution
main() {
    echo "========================================="
    echo "Fault Snapshot Generator"
    echo "========================================="

    # Find pvpython
    print_info "Looking for pvpython..."
    PVPYTHON=$(find_pvpython)

    if [ $? -ne 0 ]; then
        print_error "pvpython not found!"
        print_error "Please install ParaView and ensure pvpython is in your PATH"
        exit 1
    fi

    print_info "Found pvpython: $PVPYTHON"

    # Check for state file - use fallback if local not found
    if [ ! -f "$STATE_FILE" ] && [ -f "$FALLBACK_STATE_FILE" ]; then
        STATE_FILE="$FALLBACK_STATE_FILE"
    fi

    # Execute based on mode
    case $MODE in
        single)
            print_info "Mode: Single snapshot generation"

            if [ -n "$INPUT_FILE" ]; then
                # Mode: Load from XDMF file
                if [ ! -f "$INPUT_FILE" ]; then
                    print_error "Input file not found: $INPUT_FILE"
                    exit 1
                fi

                print_info "Input file: $INPUT_FILE"
                [ -n "$STATE_FILE" ] && print_info "State file: $STATE_FILE"
                print_info "Output directory: $OUTPUT_DIR"
                print_info "Resolution: $RESOLUTION"
                print_info "Field: $FIELD"

                ARGS="-i \"$INPUT_FILE\" -o \"$OUTPUT_DIR\" -r $RESOLUTION -f \"$FIELD\""
                [ -n "$STATE_FILE" ] && ARGS="$ARGS -s \"$STATE_FILE\""

                eval $PVPYTHON "$SCRIPT_DIR/generate_fault_snapshots.py" $ARGS

            else
                # Mode: Load from state file (legacy)
                if [ ! -f "$STATE_FILE" ]; then
                    print_error "State file not found: $STATE_FILE"
                    print_info "Either provide --input (XDMF file) or --state (state file)"
                    exit 1
                fi

                print_info "State file: $STATE_FILE"
                print_info "Output directory: $OUTPUT_DIR"
                print_info "Resolution: $RESOLUTION"
                print_info "Field: $FIELD"

                $PVPYTHON "$SCRIPT_DIR/generate_fault_snapshots.py" \
                    -s "$STATE_FILE" \
                    -o "$OUTPUT_DIR" \
                    -r $RESOLUTION \
                    -f "$FIELD"
            fi
            ;;

        batch)
            print_info "Mode: Batch processing"

            if [ -n "$INPUT_FILE" ]; then
                # Batch process single XDMF file (all timesteps)
                if [ ! -f "$INPUT_FILE" ]; then
                    print_error "Input file not found: $INPUT_FILE"
                    exit 1
                fi

                print_info "Input file: $INPUT_FILE"
                print_info "Output directory: $OUTPUT_DIR"
                print_info "Resolution: $RESOLUTION"

                ARGS="-i \"$INPUT_FILE\" -o \"$OUTPUT_DIR\" -r $RESOLUTION"
                [ -n "$STATE_FILE" ] && ARGS="$ARGS -s \"$STATE_FILE\""

                eval $PVPYTHON "$SCRIPT_DIR/batch_fault_snapshots.py" $ARGS

            elif [ -n "$DATA_DIR" ]; then
                # Batch process directory
                if [ ! -d "$DATA_DIR" ]; then
                    print_error "Data directory not found: $DATA_DIR"
                    exit 1
                fi

                print_info "Data directory: $DATA_DIR"
                print_info "File pattern: $PATTERN"
                print_info "Output directory: $OUTPUT_DIR"
                print_info "Resolution: $RESOLUTION"

                $PVPYTHON "$SCRIPT_DIR/batch_fault_snapshots.py" \
                    -d "$DATA_DIR" \
                    -p "$PATTERN" \
                    -o "$OUTPUT_DIR" \
                    -r $RESOLUTION

            else
                print_error "Either --input or --data-dir required for batch mode"
                exit 1
            fi
            ;;

        timesteps)
            print_info "Mode: Timestep processing"
            print_info "State file: $STATE_FILE"
            print_info "Output directory: $OUTPUT_DIR"
            print_info "Resolution: $RESOLUTION"

            $PVPYTHON "$SCRIPT_DIR/batch_fault_snapshots.py" \
                --state-mode \
                -s "$STATE_FILE" \
                -o "$OUTPUT_DIR" \
                -r $RESOLUTION
            ;;

        *)
            print_error "Invalid mode: $MODE"
            print_error "Valid modes: single, batch, timesteps"
            exit 1
            ;;
    esac

    if [ $? -eq 0 ]; then
        print_info "Snapshots generated successfully!"
        print_info "Output location: $OUTPUT_DIR"
    else
        print_error "Snapshot generation failed!"
        exit 1
    fi
}

# Run main function
main
