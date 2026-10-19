#!/usr/bin/env python3
"""Run the paper's analyses in the order the paper presents them.

The script numbers are historical citation keys, not a sequence; this runner
supplies the sequence. Every step is a separate process, so a failure in one
does not stop the others.
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent

GROUPS = {
    "section2": [("Fig1a-d", "code/cellline1_bifurcation_suite/code2_bifurcation_analysis_suite.py", []),
                 ("Fig1f",   "code/cellline1_bifurcation_suite/code3_discrimination_diagram.py", [])],
    "section3": [("Eq3",    "code/code52_TrackC2_discrimination_diagram.py",
                  ["results/code52_TrackC2_results.json"])],
    "section4": [("Fisher", "code/code51r3_p53_fisher_audit.py",
                  ["results/code51r3_p53_Fisher_audit_results.json"]),
                 ("spectra","code/code51r1_p53_fisher_audit.py",
                  ["results/code51_p53_Fisher_audit_results.json"]),
                 ("spectra2","code/code51r2_p53_fisher_audit.py",
                  ["results/code51r2_p53_Fisher_audit_results.json"]),
                 ("wrapup-A", "code/code51r4_p53_fisher_wrapup.py",
                  ["results/code51r4A_noise_calibration_robustness.json"], ["A"]),
                 ("wrapup-fig", "code/code51r4_p53_fisher_wrapup.py",
                  ["results/code51r4_wrapup_fourpanel.png"], ["FIG"]),
                 ("cascade","code/code81_p53_full_cascade_audit.py",
                  ["results/code81_cascade_audit_results.json"])],
    "section5": [("mechanism","code/code82_p53_mechanism_model_audit.py",
                  ["results/code82_mechanism_audit_results.json"]),
                 ("period",  "code/code89_p53_molecular_constant_decomposition.py",
                  ["results/molecular_layer/code89_molecular_decomposition_results.json"]),
                 ("links",   "code/code90_p53_full_cascade_consistency_audit.py",
                  ["results/molecular_layer/code90_cascade_audit_results.json"]),
                 ("width",   "code/code91_K2_pulse_width_decomposition.py",
                  ["results/molecular_layer/code91_pulse_width_results.json"]),
                 ("closedform","code/code92_belgacem_closedform_crossvalidation.py",
                  ["results/molecular_layer/code92_closedform_crossvalidation_results.json"]),
                 ("tcs",     "code/code93_tcs_logistic_correspondence.py",
                  ["results/molecular_layer/code93_tcs_logistic_correspondence_results.json"])],
    "extra": [("stratify","code/code13_p53_fisher_stratification_S7.py",
               ["results/p53_channel_census.csv"]),
              ("variant",  "code/code83_why_digital_variational_competition.py",
               ["results/code83_variational_competition_results.json"]),
              ("discovery","code/code84_discovery_cascade_model.py",
               ["results/code84_discovery_cascade_results.json"]),
              ("heterog",  "code/code85_heterogeneity_recheck.py",
               ["results/code85_heterogeneity_recheck_results.json"]),
              ("foci",     "code/code86_foci_persistence_trigger.py",
               ["results/code86_foci_persistence_trigger_results.json"]),
              ("clock",    "code/code87_clock_dose_independent_architecture.py",
               ["results/code87_clock_dose_independent_results.json"]),
              ("terminate","code/code88_foci_termination_channel.py",
               ["results/code88_foci_termination_channel_results.json"]),
              ("pop",      "code/cellline1_bifurcation_suite/code1_p53_stochastic_population_simulation.py", [])],
    "support": [("selfcheck","code/code50_p53_cascade_framework_selfcheck.py", [])],
}

ORDER = ["section2", "section3", "section4", "section5", "extra", "support"]


def run(step):
    label, script, outputs = step[0], step[1], step[2]
    argv = step[3] if len(step) > 3 else []
    path = REPO / script
    if not path.exists():
        return label, "MISSING SCRIPT", 0.0
    t0 = time.time()
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUTF8="1")
    proc = subprocess.run([sys.executable, str(path)] + list(argv), cwd=str(REPO),
                          capture_output=True, text=True,
                          encoding="utf-8", errors="replace", env=env)
    dt = time.time() - t0
    ok = proc.returncode == 0
    missing = [o for o in outputs if not (REPO / o).exists()]
    if ok and missing:
        ok = False
        status = "exit 0 but missing " + ", ".join(missing)
    elif ok:
        status = "ok"
    else:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-1:]
        status = "FAILED " + (tail[0] if tail else "")
    return label, status, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("groups", nargs="*", default=None,
                    help="subset of: " + ", ".join(ORDER))
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for g in ORDER:
            print(f"{g:10s} " + "; ".join(s[0] for s in GROUPS[g]))
        return 0
    groups = args.groups or ORDER
    bad = 0
    for g in groups:
        if g not in GROUPS:
            print(f"unknown group {g}; try --list")
            return 2
        print(f"\n== {g} ==")
        for step in GROUPS[g]:
            label, status, dt = run(step)
            print(f"  {label:14s} {dt:7.1f}s  {status}")
            if status != "ok":
                bad += 1
    print(f"\n{bad} step(s) not ok")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
