"""Single external entry point for the active HTL-DTA program."""

from __future__ import annotations

import importlib
import sys

from scripts import preprocess


def blocked_affinity_command() -> None:
    raise SystemExit(
        "HTL-DTA training is blocked: historical H0-S/H1-S remain "
        "DATA_NOT_READY; complete the certified 24.1/27/31/37 identity "
        "presence audit and the closed historical roster first."
    )


def main() -> None:
    module_commands = {
        "audit": "scripts.audit",
        "topology-audit": "scripts.htl_topology",
        "interaction-audit": "research.shared.interaction_identifiability_audit",
        "h0-inventory": "research.a2s.a2s_h0_metadata_inventory",
        "historical-project": "research.a2s.a2s_historical_projection",
        "historical-presence": "research.a2s.a2s_historical_presence",
        "dataset-run": "research.shared.dataset_processing",
        "natural-tail-audit": "research.a2s.a2s_natural_tail_audit",
        "d0r-roster": "research.a2s.a2s_d0r",
        "a2s-bir": "research.a2s.a2s_bir",
        "a2s-scao": "research.a2s.a2s_scao",
        "prepare-cmal-data": "research.a2s.a2s_cmal_data",
        "a2s-cmal": "research.a2s.a2s_cmal",
        "a2s-trace-q1": "research.a2s.a2s_trace_stratum",
        "a2s-trace": "research.a2s.a2s_trace",
        "a2s-trace-headroom": "research.a2s.a2s_trace_headroom",
        "a2s-mode-gates": "research.a2s.a2s_mode_gates",
        "a2s-mode-generalization": "research.a2s.a2s_mode_generalization",
        "a2s-rip-r0": "research.a2s.a2s_rip_gate",
        "a2s-transfer-object": "research.a2s.a2s_transfer_object_gate",
        "a2s-nea-preconditions": "research.a2s.a2s_nea_preconditions",
    }
    blocked_commands = {"a2s-baseline", "a2s-router", "train"}
    command = sys.argv[1] if len(sys.argv) > 1 else ""
    if command not in module_commands and command not in blocked_commands and command != "prepare-a2s-validation":
        print("usage: main.py {topology-audit,h0-inventory,historical-project,historical-presence,prepare-a2s-validation,prepare-cmal-data,dataset-run,natural-tail-audit,audit,interaction-audit,a2s-bir,a2s-scao,a2s-cmal,a2s-trace-q1,a2s-trace,a2s-trace-headroom,a2s-mode-gates,a2s-mode-generalization,a2s-rip-r0,a2s-transfer-object,a2s-nea-preconditions,a2s-baseline,a2s-router,train} [options]")
        return
    if command in {"h0-inventory", "historical-project", "historical-presence", "prepare-a2s-validation"}:
        importlib.import_module("research.shared.dataset_processing").require_run_context(sys.argv)
    sys.argv = [command, *sys.argv[2:]]
    if command == "prepare-a2s-validation":
        preprocess.innovation_main()
    elif command in blocked_commands:
        blocked_affinity_command()
    else:
        importlib.import_module(module_commands[command]).main()


if __name__ == "__main__":
    main()
