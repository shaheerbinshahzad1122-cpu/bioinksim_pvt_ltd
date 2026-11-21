"""
BioInkSim: A Biochemical and Rheological Simulator for Polymer-Based Tissue Engineering and 3D Bioprinting
===========================================================================================================
A comprehensive, scientifically defensible simulator modeling the entire tissue engineering workflow:
Polymer Selection -> Bioink Formulation -> Rheology -> Crosslinking -> Syringe/Nozzle Extrusion
-> Cell Hydrodynamic Shear Safety -> 3D Scaffold Architecture -> Degradation Kinetics -> Tissue Suitability Scoring
-> Inverse Design Optimization.

Run directly with:
    python Bioinksim.py               # Starts web server and opens browser
    python Bioinksim.py --test        # Runs scientific unit and sanity verification suite
    python Bioinksim.py --cli         # Runs standalone CLI bioink simulation
    python Bioinksim.py --benchmark   # Runs literature benchmark parity report
    python Bioinksim.py --pdf [file]  # Generates publication-grade PDF laboratory report
"""

import sys
import os
import argparse
import webbrowser
import threading
import time

# Ensure project root is in sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bioinksim.data import (
    POLYMERS_DB,
    CELLS_DB,
    GROWTH_FACTORS_DB,
    NOZZLES_DB,
    CROSSLINKERS_DB,
    TISSUES_DB,
    BENCHMARKS_DB,
)
from bioinksim.engine import BioInkEngine
from bioinksim.calibration import RheologyCalibrator
from bioinksim.uq import UncertaintyEngine
from bioinksim.sensitivity import SensitivityEngine
from bioinksim.report_pdf import generate_pdf_report
from bioinksim.server import BioInkRequestHandler, ThreadingHTTPServer
from cli import run_test_suite, run_cli_simulation, run_benchmark_report, generate_cli_pdf_report


def main():
    parser = argparse.ArgumentParser(
        description="BioInkSim: Biochemical & Rheological Simulator for 3D Bioprinting"
    )
    parser.add_argument("--test", action="store_true", help="Run automated test suite and exit")
    parser.add_argument("--cli", action="store_true", help="Run standalone simulation in terminal and exit")
    parser.add_argument("--benchmark", action="store_true", help="Run literature parity benchmarks and exit")
    parser.add_argument("--pdf", nargs="?", const="BioInkSim_Laboratory_Report.pdf", help="Generate PDF lab report and exit")
    parser.add_argument("--port", type=int, default=8501, help="Port to bind HTTP server (default: 8501)")
    parser.add_argument("--no-browser", action="store_true", help="Do not automatically open web browser")
    args = parser.parse_args()

    if args.test:
        run_test_suite()
        sys.exit(0)

    if args.cli:
        run_cli_simulation()
        sys.exit(0)

    if args.benchmark:
        run_benchmark_report()
        sys.exit(0)

    if args.pdf:
        generate_cli_pdf_report(filename=args.pdf)
        sys.exit(0)

    # Start multi-threaded HTTP Server
    port = args.port
    server = None

    for attempt_port in range(port, port + 10):
        try:
            server = ThreadingHTTPServer(("127.0.0.1", attempt_port), BioInkRequestHandler)
            port = attempt_port
            break
        except OSError:
            continue

    if not server:
        print(f"[ERROR] Could not bind to port {port} or next 10 ports. Please check running processes.")
        sys.exit(1)

    url = f"http://localhost:{port}"
    print("=================================================================")
    print("   BIOINKSIM v2.5: BIOCHEMICAL & RHEOLOGICAL SIMULATOR           ")
    print("=================================================================")
    print(f" [OK] Localhost server active at: {url}")
    print(" [OK] Modular calculation engine, calibration studio & PDF export active.")
    print(" [OK] Press Ctrl+C in this terminal window to stop the server.")
    print("=================================================================")

    if not args.no_browser:
        def open_browser():
            time.sleep(1.0)
            try:
                webbrowser.open(url)
            except Exception:
                pass
        threading.Thread(target=open_browser, daemon=True).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n[INFO] Shutting down BioInkSim server. Goodbye!")
        server.shutdown()
        server.server_close()


if __name__ == "__main__":
    main()
