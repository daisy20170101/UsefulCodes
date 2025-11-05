#!/usr/bin/env python3
"""
Visualize seismic moment release data for major earthquakes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

def load_data(filename='major_earthquakes_moment_data.csv'):
    """Load earthquake data"""
    df = pd.read_csv(filename)
    df['time'] = pd.to_datetime(df['time'])
    return df

def plot_magnitude_time_series(df, output_file='magnitude_time_series.png'):
    """Plot magnitude vs time"""
    fig, ax = plt.subplots(figsize=(14, 6))

    # Color code by magnitude
    colors = []
    for mag in df['magnitude']:
        if mag >= 9.0:
            colors.append('red')
        elif mag >= 8.5:
            colors.append('orange')
        elif mag >= 8.0:
            colors.append('gold')
        else:
            colors.append('blue')

    # Scatter plot with size proportional to moment
    scatter = ax.scatter(df['time'], df['magnitude'],
                        c=colors, s=df['seismic_moment_Nm']/1e20,
                        alpha=0.7, edgecolors='black', linewidth=0.5)

    ax.set_xlabel('Year', fontsize=12, fontweight='bold')
    ax.set_ylabel('Magnitude (Mw)', fontsize=12, fontweight='bold')
    ax.set_title('Major Earthquakes (M 7.5-9.5) Timeline', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(7.0, 10.0)

    # Format x-axis
    ax.xaxis.set_major_locator(mdates.YearLocator(20))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    plt.xticks(rotation=45)

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
               markersize=10, label='M ≥ 9.0'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='orange',
               markersize=10, label='8.5 ≤ M < 9.0'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gold',
               markersize=10, label='8.0 ≤ M < 8.5'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue',
               markersize=10, label='M < 8.0')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def plot_cumulative_moment(df, output_file='cumulative_moment.png'):
    """Plot cumulative moment release"""
    fig, ax = plt.subplots(figsize=(14, 6))

    ax.plot(df['time'], df['cumulative_moment_Nm'], 'b-', linewidth=2)

    # Mark major events (M ≥ 8.5)
    major = df[df['magnitude'] >= 8.5]
    ax.scatter(major['time'], major['cumulative_moment_Nm'],
              c='red', s=100, zorder=5, edgecolors='black', linewidth=1)

    ax.set_xlabel('Year', fontsize=12, fontweight='bold')
    ax.set_ylabel('Cumulative Seismic Moment (N·m)', fontsize=12, fontweight='bold')
    ax.set_title('Cumulative Seismic Moment Release (M 7.5+)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

    # Format x-axis
    ax.xaxis.set_major_locator(mdates.YearLocator(20))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def plot_moment_release_rate(df, output_file='moment_release_rate.png'):
    """Plot moment release rate"""
    fig, ax = plt.subplots(figsize=(14, 6))

    # Remove NaN values
    df_clean = df.dropna(subset=['moment_release_rate_Nm_per_year'])

    ax.plot(df_clean['time'], df_clean['moment_release_rate_Nm_per_year'],
           'g-', linewidth=1.5, alpha=0.7)

    # Mark major events
    major = df_clean[df_clean['magnitude'] >= 8.5]
    ax.scatter(major['time'], major['moment_release_rate_Nm_per_year'],
              c='red', s=100, zorder=5, edgecolors='black', linewidth=1,
              label='M ≥ 8.5 events')

    ax.set_xlabel('Year', fontsize=12, fontweight='bold')
    ax.set_ylabel('Moment Release Rate (N·m/year)', fontsize=12, fontweight='bold')
    ax.set_title('Seismic Moment Release Rate (1-year moving window)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
    ax.legend(fontsize=10)

    # Format x-axis
    ax.xaxis.set_major_locator(mdates.YearLocator(20))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def plot_gutenberg_richter(output_file='gutenberg_richter.png'):
    """Plot Gutenberg-Richter distribution"""
    gr_df = pd.read_csv('gutenberg_richter_distribution.csv')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Linear plot
    ax1.semilogy(gr_df['magnitude_threshold'], gr_df['cumulative_count'],
                'bo-', linewidth=2, markersize=6)
    ax1.set_xlabel('Magnitude (Mw)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Cumulative Number (N)', fontsize=12, fontweight='bold')
    ax1.set_title('Gutenberg-Richter Distribution', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, which='both')

    # Log-log plot with annual rate
    ax2.loglog(gr_df['magnitude_threshold'], gr_df['annual_rate'],
              'ro-', linewidth=2, markersize=6)
    ax2.set_xlabel('Magnitude (Mw)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Annual Rate (events/year)', fontsize=12, fontweight='bold')
    ax2.set_title('Annual Rate vs Magnitude', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def plot_world_map(df, output_file='earthquake_map.png'):
    """Plot earthquake locations on world map"""
    fig, ax = plt.subplots(figsize=(16, 9))

    # Color and size by magnitude
    sizes = (df['magnitude'] - 7.0) ** 3 * 5  # Scale for visibility
    colors = df['magnitude']

    scatter = ax.scatter(df['longitude'], df['latitude'],
                        c=colors, s=sizes, alpha=0.6,
                        cmap='hot_r', edgecolors='black', linewidth=0.5)

    ax.set_xlabel('Longitude', fontsize=12, fontweight='bold')
    ax.set_ylabel('Latitude', fontsize=12, fontweight='bold')
    ax.set_title('Epicenter Locations of Major Earthquakes (M 7.5-9.5)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Magnitude (Mw)', fontsize=11, fontweight='bold')

    # Add plate boundaries (simplified)
    # Pacific Ring of Fire
    ax.plot([-180, -120, -70, -70], [0, 40, -55, -55], 'k--', alpha=0.3, linewidth=1)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def create_summary_figure(df, output_file='summary_figure.png'):
    """Create a comprehensive summary figure"""
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

    # 1. Magnitude time series
    ax1 = fig.add_subplot(gs[0, :])
    colors = ['red' if m >= 9.0 else 'orange' if m >= 8.5 else 'gold' if m >= 8.0 else 'blue'
              for m in df['magnitude']]
    ax1.scatter(df['time'], df['magnitude'], c=colors, s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
    ax1.set_ylabel('Magnitude (Mw)', fontweight='bold')
    ax1.set_title('A) Magnitude Timeline', fontweight='bold', loc='left')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(7.0, 10.0)

    # 2. Cumulative moment
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(df['time'], df['cumulative_moment_Nm'], 'b-', linewidth=2)
    ax2.set_ylabel('Cumulative Moment (N·m)', fontweight='bold')
    ax2.set_title('B) Cumulative Moment Release', fontweight='bold', loc='left')
    ax2.grid(True, alpha=0.3)
    ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

    # 3. Moment release rate
    ax3 = fig.add_subplot(gs[1, 1])
    df_clean = df.dropna(subset=['moment_release_rate_Nm_per_year'])
    ax3.plot(df_clean['time'], df_clean['moment_release_rate_Nm_per_year'], 'g-', linewidth=1.5)
    ax3.set_ylabel('Release Rate (N·m/year)', fontweight='bold')
    ax3.set_title('C) Moment Release Rate', fontweight='bold', loc='left')
    ax3.grid(True, alpha=0.3)
    ax3.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

    # 4. World map
    ax4 = fig.add_subplot(gs[2, 0])
    sizes = (df['magnitude'] - 7.0) ** 2 * 10
    scatter = ax4.scatter(df['longitude'], df['latitude'],
                         c=df['magnitude'], s=sizes, alpha=0.6,
                         cmap='hot_r', edgecolors='black', linewidth=0.5)
    ax4.set_xlabel('Longitude', fontweight='bold')
    ax4.set_ylabel('Latitude', fontweight='bold')
    ax4.set_title('D) Epicenter Locations', fontweight='bold', loc='left')
    ax4.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax4, label='Magnitude')

    # 5. Magnitude histogram
    ax5 = fig.add_subplot(gs[2, 1])
    bins = np.arange(7.0, 10.0, 0.2)
    ax5.hist(df['magnitude'], bins=bins, color='steelblue', edgecolor='black', alpha=0.7)
    ax5.set_xlabel('Magnitude (Mw)', fontweight='bold')
    ax5.set_ylabel('Frequency', fontweight='bold')
    ax5.set_title('E) Magnitude Distribution', fontweight='bold', loc='left')
    ax5.grid(True, alpha=0.3, axis='y')

    fig.suptitle('Seismic Moment Release Analysis: Major Earthquakes (M 7.5-9.5)',
                fontsize=16, fontweight='bold', y=0.995)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def main():
    """Main function"""
    print("="*70)
    print("SEISMIC MOMENT DATA VISUALIZATION")
    print("="*70)

    # Load data
    print("\nLoading earthquake data...")
    df = load_data('major_earthquakes_moment_data.csv')
    print(f"Loaded {len(df)} earthquakes")

    # Create plots
    print("\nGenerating plots...")
    plot_magnitude_time_series(df)
    plot_cumulative_moment(df)
    plot_moment_release_rate(df)
    plot_gutenberg_richter()
    plot_world_map(df)
    create_summary_figure(df)

    print("\n" + "="*70)
    print("DONE! Generated 6 figures:")
    print("  1. magnitude_time_series.png")
    print("  2. cumulative_moment.png")
    print("  3. moment_release_rate.png")
    print("  4. gutenberg_richter.png")
    print("  5. earthquake_map.png")
    print("  6. summary_figure.png")
    print("="*70)

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
        print("\nNote: Matplotlib is required for visualization.")
        print("Install with: pip install matplotlib")
