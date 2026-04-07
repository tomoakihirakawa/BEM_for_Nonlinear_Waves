#!/usr/bin/env python3
"""
Wave spectrum visualization tool.

Uses the C++ wave spectrum implementation (via pybind11) directly.

Usage:
  python plot_spectrum.py --type jonswap --Hm0 7.4 --Tp 12.0 --gamma 3.3 --h 200
  python plot_spectrum.py --type bretschneider --H13 2.0 --T13 8.0 --h 50
  python plot_spectrum.py --type jonswap --Hm0 7.4 --Tp 12.0 --gamma 3.3 --h 200 --N 50
"""

import argparse
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

import wave_spectrum as ws


def create_wave(args):
    """Create RandomWaterWaveTheory from command-line arguments."""
    if args.type == "jonswap":
        if args.Hm0 is not None and args.Tp is not None:
            return ws.RandomWaterWaveTheory.JONSWAP_Hm0_Tp(
                args.Hm0, args.Tp, args.gamma, args.h, args.bottom_z)
        elif args.H13 is not None and args.T13 is not None:
            return ws.RandomWaterWaveTheory.JONSWAP_H13_T13(
                args.H13, args.T13, args.gamma, args.h, args.bottom_z)
        else:
            sys.exit("JONSWAP requires (--Hm0, --Tp) or (--H13, --T13)")
    elif args.type == "bretschneider":
        if args.H13 is not None and args.T13 is not None:
            return ws.RandomWaterWaveTheory.Bretschneider(
                args.H13, args.T13, args.h, args.bottom_z)
        else:
            sys.exit("Bretschneider requires --H13 and --T13")
    else:
        sys.exit(f"Unknown spectrum type: {args.type}")


def plot_spectrum_and_components(wave, N_display=None):
    """Plot continuous spectrum and discretized components."""
    comp = wave.get_components()
    f_comp = np.array(comp["f"])
    A_comp = np.array(comp["A"])

    # Continuous spectrum curve (fine resolution)
    f_fine = np.linspace(wave.f_min * 0.8, wave.f_max * 1.1, 2000)
    f_fine = f_fine[f_fine > 0]
    S_fine = wave.spectrum_array(f_fine)

    # Discrete spectrum: S_i = A_i^2 / (2 * df)
    S_comp = A_comp**2 / (2 * wave.df)

    # Hm0 recovery check
    m0_continuous = np.trapezoid(S_fine, f_fine)
    m0_discrete = np.sum(A_comp**2 / 2)
    Hm0_continuous = 4 * np.sqrt(m0_continuous)
    Hm0_discrete = 4 * np.sqrt(m0_discrete)

    # Actual input wave height
    Hm0_input = wave.Hm0 if wave.Hm0 > 0 else wave.H13

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"{wave.spectrum_type.name}  |  mode={wave.mode.name}\n"
                 f"gamma={wave.gamma:.1f}  h={wave.h:.0f} m  "
                 f"N={wave.N} components", fontsize=13)

    # ---- (1) Spectrum S(f) ----
    ax = axes[0, 0]
    ax.plot(f_fine, S_fine, "k-", linewidth=1.5, label="S(f) continuous")
    ax.bar(f_comp, S_comp, width=wave.df * 0.8, alpha=0.4, color="steelblue",
           label=f"Discretized (N={wave.N})")
    ax.set_xlabel("f [Hz]")
    ax.set_ylabel("S(f) [m²/Hz]")
    ax.set_title("Spectrum")
    ax.legend(fontsize=9)
    ax.set_xlim(0, wave.f_max * 1.15)

    # Hm0 annotation
    ax.annotate(
        f"Hm0 (input): {Hm0_input:.3f} m\n"
        f"Hm0 (continuous): {Hm0_continuous:.3f} m\n"
        f"Hm0 (discrete N={wave.N}): {Hm0_discrete:.3f} m",
        xy=(0.97, 0.95), xycoords="axes fraction",
        ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8))

    # ---- (2) Component amplitudes ----
    ax = axes[0, 1]
    ax.fill_between(f_comp, 0, A_comp, alpha=0.4, color="steelblue", label="A_i (C++)")
    ax.plot(f_comp, A_comp, "steelblue", linewidth=0.8)
    # Overlay theoretical: A = sqrt(2 * S(f) * df)
    A_theory = np.sqrt(2 * S_fine * wave.df)
    ax.plot(f_fine, A_theory, "k--", linewidth=1.0, alpha=0.7, label=r"$\sqrt{2\,S(f)\,\Delta f}$")
    ax.set_xlabel("f [Hz]")
    ax.set_ylabel("A [m]")
    ax.set_title("Component amplitudes")
    ax.legend(fontsize=9)
    ax.set_xlim(0, wave.f_max * 1.15)

    # ---- (3) Coarser discretization comparison ----
    ax = axes[1, 0]
    if N_display is not None and N_display < wave.N:
        # Sub-sample: group components into N_display bins
        step = wave.N // N_display
        f_coarse = []
        A_coarse = []
        for i in range(0, wave.N, step):
            group = slice(i, min(i + step, wave.N))
            f_group = f_comp[group]
            A_group = A_comp[group]
            f_coarse.append(np.mean(f_group))
            # Energy-preserving: A_combined = sqrt(sum(A_i^2))
            A_coarse.append(np.sqrt(np.sum(A_group**2)))
        f_coarse = np.array(f_coarse)
        A_coarse = np.array(A_coarse)
        df_coarse = (wave.f_max - wave.f_min) / N_display
        S_coarse = A_coarse**2 / (2 * df_coarse)

        m0_coarse = np.sum(A_coarse**2 / 2)
        Hm0_coarse = 4 * np.sqrt(m0_coarse)

        ax.plot(f_fine, S_fine, "k-", linewidth=1.0, alpha=0.5, label="S(f)")
        ax.bar(f_coarse, S_coarse, width=df_coarse * 0.8, alpha=0.5,
               color="coral", label=f"N={N_display}")
        ax.bar(f_comp, S_comp, width=wave.df * 0.8, alpha=0.2,
               color="steelblue", label=f"N={wave.N}")
        ax.set_xlabel("f [Hz]")
        ax.set_ylabel("S(f) [m²/Hz]")
        ax.set_title(f"Discretization comparison: N={N_display} vs N={wave.N}")
        ax.legend(fontsize=9)
        ax.set_xlim(0, wave.f_max * 1.15)
        ax.annotate(f"Hm0 (N={N_display}): {Hm0_coarse:.3f} m",
                    xy=(0.97, 0.95), xycoords="axes fraction",
                    ha="right", va="top", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8))
    else:
        ax.text(0.5, 0.5, "Use --N <number> to compare\ncoarser discretization",
                transform=ax.transAxes, ha="center", va="center", fontsize=12, color="gray")
        ax.set_title("Discretization comparison")

    # ---- (4) Time series sample ----
    ax = axes[1, 1]
    T_rep = 1.0 / f_comp[np.argmax(A_comp)]  # peak period from components
    n_waves = 200  # enough waves for statistical convergence
    t_end = n_waves * T_rep
    dt = T_rep / 20  # ~20 points per peak period
    t = np.arange(0, t_end, dt)
    eta = wave.eta_array([0, 0, 0], t)
    eta_surface = eta - wave.h - wave.bottom_z  # relative to still water

    ax.plot(t, eta_surface, "b-", linewidth=0.3)
    ax.axhline(0, color="k", linewidth=0.3)
    ax.set_xlabel("t [s]")
    ax.set_ylabel("eta [m]")
    ax.set_title(f"Surface elevation at x=0 ({n_waves} waves)")

    # Check statistics
    std_eta = np.std(eta_surface)
    Hm0_timeseries = 4 * std_eta
    ax.annotate(f"Hm0 (input): {Hm0_input:.3f} m\n"
                f"Hm0 (discrete): {Hm0_discrete:.3f} m\n"
                f"4*std(eta): {Hm0_timeseries:.3f} m",
                xy=(0.97, 0.95), xycoords="axes fraction",
                ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8))

    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Wave spectrum visualization (C++ backend)")
    parser.add_argument("--type", choices=["jonswap", "bretschneider"],
                        required=True, help="Spectrum type")
    parser.add_argument("--H13", type=float, help="Significant wave height H1/3 [m]")
    parser.add_argument("--Hm0", type=float, help="Spectral significant wave height Hm0 [m]")
    parser.add_argument("--T13", type=float, help="Significant wave period T1/3 [s]")
    parser.add_argument("--Tp", type=float, help="Peak period Tp [s]")
    parser.add_argument("--gamma", type=float, default=3.3, help="JONSWAP gamma (default: 3.3)")
    parser.add_argument("--h", type=float, required=True, help="Water depth [m]")
    parser.add_argument("--bottom_z", type=float, default=0.0, help="Bottom z-coordinate (default: 0)")
    parser.add_argument("--N", type=int, default=None,
                        help="Coarser discretization count for comparison plot")
    args = parser.parse_args()

    wave = create_wave(args)
    print(wave)
    plot_spectrum_and_components(wave, N_display=args.N)


if __name__ == "__main__":
    main()
