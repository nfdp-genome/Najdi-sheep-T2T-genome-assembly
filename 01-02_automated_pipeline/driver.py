#!/usr/bin/env python3
"""
Najdi assembly pipeline driver: config-driven SLURM (IBEX) or local execution.
See docs/command_inventory.md and config.example.yaml.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import shlex
import subprocess
import sys
import tempfile
import textwrap
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import yaml

PIPELINE_ROOT = Path(__file__).resolve().parent
STEPS_DIR = PIPELINE_ROOT / "steps"


@dataclass
class Job:
    job_id: str
    step_name: str
    script: Path
    env: Dict[str, str]
    resources_key: str
    outputs: List[Path] = field(default_factory=list)
    depends_on: List[str] = field(default_factory=list)


def load_config(path: Path) -> dict:
    with path.open() as f:
        return yaml.safe_load(f)


def _mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _ts() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def resources_for(cfg: dict, key: str) -> dict:
    r = dict(cfg.get("resources", {}).get("defaults", {}))
    overrides = cfg.get("resources", {}).get("per_stage", {})
    if key in overrides:
        r.update(overrides[key])
    return r


def topological_waves(jobs: List[Job]) -> List[List[Job]]:
    by_id = {j.job_id: j for j in jobs}
    pending: Set[str] = set(by_id)
    rev_deps: Dict[str, Set[str]] = defaultdict(set)
    dep_count: Dict[str, int] = {}
    for j in jobs:
        dep_count[j.job_id] = len(j.depends_on)
        for d in j.depends_on:
            rev_deps[d].add(j.job_id)

    waves: List[List[Job]] = []
    while pending:
        ready = [by_id[jid] for jid in pending if dep_count[jid] == 0]
        if not ready:
            raise RuntimeError("Cycle or missing dependency in job graph")
        waves.append(sorted(ready, key=lambda x: x.job_id))
        for j in ready:
            pending.remove(j.job_id)
        for j in ready:
            for child in rev_deps[j.job_id]:
                dep_count[child] -= 1
    return waves


def write_failure_row(
    failures_tsv: Path,
    row: dict,
) -> None:
    failures_tsv.parent.mkdir(parents=True, exist_ok=True)
    exists = failures_tsv.exists()
    with failures_tsv.open("a", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "timestamp",
                "job_id",
                "assembly_id",
                "step",
                "exit_code",
                "slurm_job_id",
                "log_path",
                "tail_stderr",
            ],
        )
        if not exists:
            w.writeheader()
        w.writerow(row)


def append_manifest(manifest_path: Path, entry: dict) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if manifest_path.exists() else "w"
    with manifest_path.open(mode) as f:
        f.write(json.dumps(entry) + "\n")


def load_state(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def save_state(path: Path, state: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2))


def outputs_ready(outputs: List[Path]) -> bool:
    return all(p.is_file() and p.stat().st_size > 0 for p in outputs)


def slurm_wrap(
    job: Job,
    cfg: dict,
    log_path: Path,
    modules: Optional[str],
) -> str:
    res = resources_for(cfg, job.resources_key)
    cpus = int(res.get("cpus", 1))
    mem = res.get("memory", "8G")
    wall = res.get("walltime", "01:00:00")
    nodes = int(res.get("nodes", 1))
    ex = cfg.get("execution", {})
    slurm = ex.get("slurm", {})
    partition = slurm.get("partition", "batch")
    account = slurm.get("account")
    constraint = slurm.get("constraint")

    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=najdi_{job.step_name[:40]}",
        f"#SBATCH --output={log_path}",
        f"#SBATCH --error={log_path}",
        f"#SBATCH --time={wall}",
        f"#SBATCH --mem={mem}",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --nodes={nodes}",
        f"#SBATCH --partition={partition}",
    ]
    if constraint:
        lines.append(f"#SBATCH --constraint={constraint}")
    if account:
        lines.append(f"#SBATCH --account={account}")
    for extra in slurm.get("sbatch_extra") or []:
        lines.append(f"#SBATCH {extra}")

    lines.append("set -euo pipefail")
    if modules:
        lines.append("module purge 2>/dev/null || true")
        for m in modules.split():
            lines.append(f"module load {shlex.quote(m)} || true")

    env_exports = []
    for k, v in job.env.items():
        if v is None:
            continue
        env_exports.append(f"export {k}={shlex.quote(str(v))}")
    lines.append("\n".join(env_exports))
    lines.append(f"bash {shlex.quote(str(job.script))}")
    return "\n".join(lines) + "\n"


def run_local(
    job: Job,
    log_path: Path,
    use_docker: bool,
    docker_image: Optional[str],
    docker_binds: List[str],
    cpus: int,
    memory: str,
) -> int:
    _mkdir(log_path.parent)
    env = os.environ.copy()
    env.update({k: str(v) for k, v in job.env.items() if v is not None})

    if use_docker and docker_image:
        binds = []
        for b in docker_binds:
            binds.extend(["-v", f"{b}:{b}"])
        mem_arg = memory.rstrip("Gg") if memory else "8"
        cmd = [
            "docker",
            "run",
            "--rm",
            f"--cpus={cpus}",
            f"--memory={memory}",
            *binds,
            "-w",
            str(PIPELINE_ROOT),
            "-e",
            "HOME=/tmp",
            docker_image,
            "bash",
            "-c",
            "set -euo pipefail; " + " ; ".join(f"export {k}={shlex.quote(str(v))}" for k, v in job.env.items() if v is not None) + f" ; bash {shlex.quote(str(job.script))}",
        ]
        proc = subprocess.run(cmd, stdout=log_path.open("w"), stderr=subprocess.STDOUT)
        return proc.returncode

    with log_path.open("w") as logf:
        proc = subprocess.run(
            ["bash", str(job.script)],
            env=env,
            stdout=logf,
            stderr=subprocess.STDOUT,
        )
    return proc.returncode


def tail_file(path: Path, n: int) -> str:
    if not path.exists():
        return ""
    lines = path.read_text(errors="replace").splitlines()
    return "\n".join(lines[-n:])


def append_qc_bundle(
    jobs: List[Job],
    cfg: dict,
    workdir: Path,
    *,
    job_id: str,
    label: str,
    fasta: Path,
    depends_on: List[str],
    qc_subtag: str,
) -> None:
    stages = cfg.get("stages", {})
    if not stages.get("intermediate_qc", True) or not stages.get("qc_waves", True):
        return
    log_root = workdir / (cfg.get("logging", {}).get("dir") or "logs") / "qc" / qc_subtag
    busco_root = workdir / "qc" / qc_subtag / "busco"
    quast_root = workdir / "qc" / qc_subtag / "quast"
    jobs.append(
        Job(
            job_id=job_id,
            step_name="qc_bundle",
            script=STEPS_DIR / "step_qc_bundle.sh",
            env=_qc_env(cfg, str(fasta), label, busco_root, log_root, quast_root),
            resources_key="qc",
            outputs=[],
            depends_on=depends_on.copy(),
        )
    )

HIFIASM_CONTIG_PHASES = ("primary", "hap1", "hap2")


def expand_contig_phases(cfg: dict) -> List[str]:
    h = cfg.get("hifiasm") or {}
    raw = h.get("contig_phases", "primary")
    if isinstance(raw, str) and raw.strip().lower() == "all":
        return list(HIFIASM_CONTIG_PHASES)
    if isinstance(raw, list):
        out: List[str] = []
        for x in raw:
            s = str(x).strip().lower()
            if s == "all":
                return list(HIFIASM_CONTIG_PHASES)
            out.append(str(x).strip())
        return out or ["primary"]
    s = str(raw).strip().lower()
    if s == "all":
        return list(HIFIASM_CONTIG_PHASES)
    if s in HIFIASM_CONTIG_PHASES:
        return [str(raw).strip()]
    return ["primary"]


def phase_job_suffix(phase: str, phases: List[str]) -> str:
    if len(phases) == 1 and phase == "primary":
        return ""
    return f"_{phase}"


def phase_root(workdir: Path, sample: str, phase: str, phases: List[str]) -> Path:
    if len(phases) == 1 and phase == "primary":
        return workdir
    return workdir / "phases" / f"{sample}_{phase}"


def hifiasm_gfa_for_phase(hifiasm_out: Path, prefix: str, phase: str) -> Path:
    stem = f"{prefix}.asm"
    if phase == "primary":
        return hifiasm_out / f"{stem}.bp.p_ctg.gfa"
    return hifiasm_out / f"{stem}.bp.{phase}.p_ctg.gfa"


def hic_bam_for_phase(cfg: dict, phase: str) -> str:
    inputs = cfg.get("inputs") or {}
    m = inputs.get("hic_sorted_bam_by_phase")
    if isinstance(m, dict):
        p = m.get(phase) or m.get(str(phase))
        if p:
            return str(Path(p).resolve())
    return str(Path(inputs["hic_sorted_bam"]).resolve())


def existing_contigs_for_phase(cfg: dict, phase: str) -> Optional[Path]:
    inputs = cfg.get("inputs") or {}
    ex = inputs.get("existing_contigs_fasta")
    if isinstance(ex, dict):
        p = ex.get(phase)
        if p:
            return Path(p).resolve()
    if phase == "primary" and inputs.get("existing_primary_fasta"):
        return Path(inputs["existing_primary_fasta"]).resolve()
    return None


def gfa_job_id(sample: str, phase: str, phases: List[str]) -> str:
    if len(phases) == 1 and phase == "primary":
        return "gfa_primary"
    return f"gfa_{phase}_{sample}"


def qc_first_contigs_job_id(sample: str, phase: str, phases: List[str]) -> str:
    if len(phases) == 1 and phase == "primary":
        return "qc_primary"
    return f"qc_contigs_{phase}_{sample}"


def scaffold_job_id(sc: str, sample: str, phase: str, phases: List[str]) -> str:
    ps = phase_job_suffix(phase, phases)
    return f"scaffold_{sc}_{sample}{ps}"


def build_jobs(cfg: dict) -> List[Job]:
    jobs: List[Job] = []
    proj = cfg["project"]
    workdir = Path(proj["workdir"]).resolve()
    sample = str(proj["sample_id"])
    stages = cfg.get("stages", {})
    inputs = cfg["inputs"]
    run_mode = cfg.get("run_mode", "full_optimization")
    hh = cfg.get("hifiasm") or {}
    prefix = hh.get("output_prefix", f"{sample}_hifiasm")
    phases = expand_contig_phases(cfg)

    if run_mode == "selective":
        for i, spec in enumerate(cfg.get("selective", {}).get("jobs", [])):
            step = spec["step"]
            script = STEPS_DIR / f"step_{step}.sh"
            if not script.exists():
                script = STEPS_DIR / f"{step}.sh"
            env = {k: str(v) for k, v in spec.get("env", {}).items()}
            if "FA" not in env and spec.get("assembly_fasta"):
                env["FA"] = spec["assembly_fasta"]
            jobs.append(
                Job(
                    job_id=f"selective_{i}_{step}",
                    step_name=step,
                    script=script,
                    env=env,
                    resources_key=spec.get("resources_key", "qc"),
                    outputs=[],
                    depends_on=[],
                )
            )
        return jobs

    hifiasm_out = workdir / "01_hifiasm" / sample
    ex_pri = inputs.get("existing_primary_fasta")
    skip_hifiasm_asm = (not stages.get("hifiasm", True)) or (
        bool(ex_pri) and len(phases) == 1 and phases[0] == "primary"
    )

    hifiasm_jid = f"hifiasm_{sample}"
    hifiasm_dep: List[str] = []
    if not skip_hifiasm_asm:
        gfa_outs = [hifiasm_gfa_for_phase(hifiasm_out, prefix, ph) for ph in phases]
        jobs.append(
            Job(
                job_id=hifiasm_jid,
                step_name="hifiasm",
                script=STEPS_DIR / "step_hifiasm.sh",
                env={
                    "HIFIASM_OUTDIR": str(hifiasm_out),
                    "HIFIASM_PREFIX": prefix,
                    "THREADS": str(resources_for(cfg, "hifiasm").get("cpus", 40)),
                    "HIFIASM_HIFI": str(Path(inputs["hifi_reads"]).resolve()),
                    "HIFIASM_ONT": str(Path(inputs["ont_reads"]).resolve()),
                },
                resources_key="hifiasm",
                outputs=gfa_outs,
                depends_on=[],
            )
        )
        hifiasm_dep = [hifiasm_jid]

    all_clean_job_ids: List[str] = []

    for phase in phases:
        phase_jobs_mark = len(jobs)
        ps = phase_job_suffix(phase, phases)
        root = phase_root(workdir, sample, phase, phases)
        gfa_path = hifiasm_gfa_for_phase(hifiasm_out, prefix, phase)
        contig_fa = root / "assemblies" / f"{sample}_hifiasm_{phase}_contigs.fasta"
        existing_fa = existing_contigs_for_phase(cfg, phase)
        skip_gfa = bool(existing_fa) or not stages.get("gfa_to_primary", True)
        contig_fa_p: Path = existing_fa if existing_fa else contig_fa

        deps_contig: List[str] = []
        if not skip_gfa:
            gid = gfa_job_id(sample, phase, phases)
            jobs.append(
                Job(
                    job_id=gid,
                    step_name="gfa_primary",
                    script=STEPS_DIR / "step_gfa_primary.sh",
                    env={"GFA_IN": str(gfa_path), "FASTA_OUT": str(contig_fa)},
                    resources_key="gfa_to_fa",
                    outputs=[contig_fa],
                    depends_on=hifiasm_dep.copy(),
                )
            )
            deps_contig = [gid]
        elif not existing_fa and hifiasm_dep:
            deps_contig = hifiasm_dep.copy()

        if existing_fa:
            primary_graph_deps: List[str] = []
        elif not skip_gfa:
            primary_graph_deps = [gfa_job_id(sample, phase, phases)]
        elif hifiasm_dep:
            primary_graph_deps = hifiasm_dep.copy()
        else:
            primary_graph_deps = []

        hic_bam = hic_bam_for_phase(cfg, phase)

        if stages.get("qc_waves", True):
            log_root = root / (cfg.get("logging", {}).get("dir") or "logs")
            qc_dir = root / "qc" / f"post_contigs_{phase}"
            qjid = qc_first_contigs_job_id(sample, phase, phases)
            jobs.append(
                Job(
                    job_id=qjid,
                    step_name="qc_bundle",
                    script=STEPS_DIR / "step_qc_bundle.sh",
                    env=_qc_env(
                        cfg,
                        str(contig_fa_p),
                        f"{sample}_{phase}_contigs",
                        qc_dir / "busco",
                        log_root / f"qc_contigs_{phase}",
                        qc_dir / "quast",
                    ),
                    resources_key="qc",
                    outputs=[],
                    depends_on=primary_graph_deps.copy(),
                )
            )

        scaff_paths: Dict[str, str] = {}
        scaffolders = stages.get("scaffolders") or []

        for sc in scaffolders:
            sid = scaffold_job_id(sc, sample, phase, phases)
            if sc == "yahs":
                out_fa = root / "assemblies" / f"{sample}_hifiasm_{phase}_hic_scaffold_yahs.fasta"
                ydir = root / "02_scaffold" / "yahs" / sample
                jobs.append(
                    Job(
                        job_id=sid,
                        step_name="yahs",
                        script=STEPS_DIR / "step_yahs.sh",
                        env={
                            "YAHS_WORKDIR": str(ydir),
                            "ASM_FASTA": str(contig_fa_p),
                            "HIC_BAM": hic_bam,
                            "YAHS_PREFIX": f"{sample}_{phase}_yahs",
                            "FINAL_SCAFFOLD_FA": str(out_fa),
                        },
                        resources_key="scaffold",
                        outputs=[out_fa],
                        depends_on=primary_graph_deps.copy(),
                    )
                )
                scaff_paths["yahs"] = str(out_fa)
            elif sc == "salsa":
                out_fa = root / "assemblies" / f"{sample}_hifiasm_{phase}_hic_scaffold_SALSA.fasta"
                sdir = root / "02_scaffold" / "salsa" / sample / phase
                jobs.append(
                    Job(
                        job_id=sid,
                        step_name="salsa",
                        script=STEPS_DIR / "step_salsa.sh",
                        env={
                            "SALSA_WORKDIR": str(sdir),
                            "REF_FASTA": str(contig_fa_p),
                            "SORTED_BAM": hic_bam,
                            "PREFIX": f"{sample}_hifiasm_{phase}_hic",
                            "SALSA_RUN": str(Path(inputs["salsa_run_pipeline"]).resolve()),
                            "FINAL_SCAFFOLD_FA": str(out_fa),
                        },
                        resources_key="scaffold",
                        outputs=[out_fa],
                        depends_on=primary_graph_deps.copy(),
                    )
                )
                scaff_paths["salsa"] = str(out_fa)

        for sc, spa in scaff_paths.items():
            append_qc_bundle(
                jobs,
                cfg,
                root,
                job_id=f"qc_scaffold_{sc}_{sample}{ps}",
                label=f"{sample}_{phase}_scaffold_{sc}",
                fasta=Path(spa),
                depends_on=[scaffold_job_id(sc, sample, phase, phases)],
                qc_subtag=f"post_scaffold_{phase}_{sc}",
            )

        polishers = stages.get("polishers") or []
        polish_outputs: List[Tuple[str, str, str, Path, str]] = []

        for sc in scaffolders:
            if sc not in scaff_paths:
                continue
            asm_scaff = scaff_paths[sc]
            tag = f"{sample}_hifiasm_{phase}_hic_scaffold_{sc}"
            scaff_dep = [scaffold_job_id(sc, sample, phase, phases)]
            for pol in polishers:
                if pol == "racon":
                    jid = f"polish_racon_{sc}_{sample}{ps}"
                    wdir = root / "03_polish" / "racon" / sc / sample
                    final_fa = wdir / f"{tag}_racon_final.fasta"
                    jobs.append(
                        Job(
                            job_id=jid,
                            step_name="polish_racon",
                            script=STEPS_DIR / "step_polish_racon.sh",
                            env={
                                "THREADS": str(resources_for(cfg, "polish").get("cpus", 40)),
                                "ONT_READS": str(Path(inputs["ont_reads"]).resolve()),
                                "ASM_IN": asm_scaff,
                                "WORK_DIR": str(wdir),
                                "TAG": tag,
                                "FINAL_OUT": str(final_fa),
                            },
                            resources_key="polish",
                            outputs=[final_fa],
                            depends_on=scaff_dep.copy(),
                        )
                    )
                    polish_outputs.append((sc, "racon", tag, final_fa, jid))
                    append_qc_bundle(
                        jobs,
                        cfg,
                        root,
                        job_id=f"qc_polish_racon_{sc}_{sample}{ps}",
                        label=f"{tag}_racon",
                        fasta=final_fa,
                        depends_on=[jid],
                        qc_subtag=f"post_polish_racon_{phase}_{sc}",
                    )
                elif pol == "medaka":
                    jid = f"polish_medaka_{sc}_{sample}{ps}"
                    wdir = root / "03_polish" / "medaka" / sc / sample
                    final_fa = wdir / f"{tag}_medaka_consensus.fasta"
                    jobs.append(
                        Job(
                            job_id=jid,
                            step_name="polish_medaka",
                            script=STEPS_DIR / "step_polish_medaka.sh",
                            env={
                                "THREADS": str(resources_for(cfg, "polish").get("cpus", 40)),
                                "ONT_READS": str(Path(inputs["ont_reads"]).resolve()),
                                "ASM_IN": asm_scaff,
                                "WORK_DIR": str(wdir),
                                "TAG": tag,
                                "FINAL_OUT": str(final_fa),
                                "MEDAKA_MODEL": inputs.get("medaka_model", "r941_min_sup_g507"),
                            },
                            resources_key="medaka",
                            outputs=[final_fa],
                            depends_on=scaff_dep.copy(),
                        )
                    )
                    polish_outputs.append((sc, "medaka", tag, final_fa, jid))
                    append_qc_bundle(
                        jobs,
                        cfg,
                        root,
                        job_id=f"qc_polish_medaka_{sc}_{sample}{ps}",
                        label=f"{tag}_medaka",
                        fasta=final_fa,
                        depends_on=[jid],
                        qc_subtag=f"post_polish_medaka_{phase}_{sc}",
                    )

        pol_job_ids: Dict[Tuple[str, str], str] = {}
        for sc, pol, tag, final_fa, pjid in polish_outputs:
            pol_job_ids[(sc, pol)] = pjid

        gapfillers = stages.get("gapfillers") or []
        ont = str(Path(inputs["ont_reads"]).resolve())

        for sc, pol, tag, pol_fa, _pj in polish_outputs:
            pol_key = (sc, pol)
            if pol_key not in pol_job_ids:
                continue
            pdeps = [pol_job_ids[pol_key]]
            racon_wdir = root / "03_polish" / "racon" / sc / sample
            asm_r1 = racon_wdir / f"{tag}_racon_r1.fasta"
            for gf in gapfillers:
                jid = f"gap_{gf}_{pol}_{sc}_{sample}{ps}"
                if gf == "ntlink":
                    wdir = root / "04_gapfill" / "ntlink" / pol / sc / sample
                    final_out = wdir / f"{tag}_{pol}_ntlink.fa"
                    asm_in = str(pol_fa) if pol == "medaka" else str(asm_r1)
                    jobs.append(
                        Job(
                            job_id=jid,
                            step_name="gap_ntlink",
                            script=STEPS_DIR / "step_gap_ntlink.sh",
                            env={
                                "THREADS": str(resources_for(cfg, "gapfill").get("cpus", 40)),
                                "NTLINK_BIN": str(inputs.get("ntlink_bin") or "ntLink"),
                                "ASM_IN": asm_in,
                                "ONT_FASTQ": ont,
                                "WORK_DIR": str(wdir),
                                "TAG": tag,
                                "POLISHER": pol,
                                "FINAL_OUT": str(final_out),
                            },
                            resources_key="gapfill",
                            outputs=[final_out],
                            depends_on=pdeps.copy(),
                        )
                    )
                    append_qc_bundle(
                        jobs, cfg, root, job_id=f"qc_{jid}", label=jid,
                        fasta=final_out, depends_on=[jid], qc_subtag=f"post_{jid}",
                    )
                elif gf == "tgsgapcloser":
                    wdir = root / "04_gapfill" / "tgsgapcloser" / pol / sc / sample
                    final_out = wdir / f"{tag}_{pol}_tgsgapcloser.fa"
                    jobs.append(
                        Job(
                            job_id=jid,
                            step_name="gap_tgsgapcloser",
                            script=STEPS_DIR / "step_gap_tgsgapcloser.sh",
                            env={
                                "THREADS": str(resources_for(cfg, "gapfill").get("cpus", 40)),
                                "ASM_IN": str(pol_fa),
                                "ONT_FASTQ": ont,
                                "WORK_DIR": str(wdir),
                                "TAG": tag,
                                "FINAL_OUT": str(final_out),
                            },
                            resources_key="gapfill",
                            outputs=[final_out],
                            depends_on=pdeps.copy(),
                        )
                    )
                    append_qc_bundle(
                        jobs, cfg, root, job_id=f"qc_{jid}", label=jid,
                        fasta=final_out, depends_on=[jid], qc_subtag=f"post_{jid}",
                    )

        if stages.get("funannotate_clean_full_tree", False):
            clean_nodes: List[Tuple[str, Path, List[str]]] = []
            clean_nodes.append(
                (f"{sample}_{phase}_contigs", Path(contig_fa_p), primary_graph_deps.copy())
            )
            for sc in scaffolders:
                if sc in scaff_paths:
                    clean_nodes.append(
                        (
                            f"{sample}_{phase}_hic_{sc}",
                            Path(scaff_paths[sc]),
                            [scaffold_job_id(sc, sample, phase, phases)],
                        )
                    )
            for j in jobs[phase_jobs_mark:]:
                if j.step_name.startswith("polish_") or j.step_name.startswith("gap_"):
                    if j.outputs:
                        clean_nodes.append((j.job_id, j.outputs[0], [j.job_id]))

            seen: Set[str] = set()
            for aid, src, deps in clean_nodes:
                key = str(src.resolve())
                if key in seen:
                    continue
                seen.add(key)
                outc = root / "06_clean" / f"{aid}.cleaned.fasta"
                cjid = f"clean_{aid}"
                jobs.append(
                    Job(
                        job_id=cjid,
                        step_name="funannotate_clean",
                        script=STEPS_DIR / "step_funannotate_clean.sh",
                        env={
                            "GENOME": str(src.resolve()),
                            "OUT_FASTA": str(outc),
                            "FUNANNOTATE_RUNNER": str(
                                inputs.get("funannotate_runner", "micromamba run -n funannotate_clean_env")
                            ),
                        },
                        resources_key="funannotate_clean",
                        outputs=[outc],
                        depends_on=deps,
                    )
                )
                all_clean_job_ids.append(cjid)
                append_qc_bundle(
                    jobs,
                    cfg,
                    root,
                    job_id=f"qc_{cjid}",
                    label=f"{aid}_cleaned",
                    fasta=outc,
                    depends_on=[cjid],
                    qc_subtag=f"post_clean_{aid}",
                )

    if stages.get("final_quast_multi", False):
        dep_ids = all_clean_job_ids if all_clean_job_ids else [j.job_id for j in jobs if j.outputs]
        list_file = workdir / "07_quast_all" / "genome_list.txt"
        labels_file = workdir / "07_quast_all" / "labels.txt"
        jobs.append(
            Job(
                job_id="quast_multi_final",
                step_name="quast_multi",
                script=STEPS_DIR / "step_quast_multi.sh",
                env={
                    "OUT_DIR": str(workdir / "07_quast_all" / "quast_report"),
                    "THREADS": str(resources_for(cfg, "quast_all").get("cpus", 16)),
                    "QUAST_FASTA_LIST": str(list_file),
                    "QUAST_LABELS_FILE": str(labels_file),
                },
                resources_key="quast_all",
                outputs=[workdir / "07_quast_all" / "quast_report" / "report.txt"],
                depends_on=dep_ids,
            )
        )

    return jobs


def _qc_env(
    cfg: dict,
    fa: str,
    label: str,
    busco_root: Path,
    log_root: Path,
    quast_root: Path,
) -> Dict[str, str]:
    inputs = cfg["inputs"]
    cpus = str(resources_for(cfg, "qc").get("cpus", 40))
    ref = inputs.get("reference_fasta")
    dotprep = inputs.get("dotprep_script")
    env = {
        "FA": fa,
        "LABEL": label,
        "L1_PATH": str(Path(inputs["busco_lineage_1"]).resolve()),
        "L2_PATH": str(Path(inputs["busco_lineage_2"]).resolve()),
        "BUSCO_OUT_ROOT": str(busco_root),
        "LOG_ROOT": str(log_root),
        "CPUS": cpus,
        "QUAST_OUT_ROOT": str(quast_root),
    }
    if ref and dotprep:
        outpre = log_root.parent / "dotplot" / label / f"{label}__vs_ref"
        _mkdir(outpre.parent)
        env["REFERENCE_FASTA"] = str(Path(ref).resolve())
        env["DOTPREP"] = str(Path(dotprep).resolve())
        env["DOTPLOT_OUTPREFIX"] = str(outpre)
        env["THREADS"] = cpus
    return env


def write_quast_inputs_files(workdir: Path) -> None:
    """Write genome_list.txt + labels.txt (sorted, same order) for final QUAST."""
    outd = workdir / "07_quast_all"
    _mkdir(outd)
    fps: List[Path] = []
    top_clean = workdir / "06_clean"
    if top_clean.is_dir():
        fps.extend(top_clean.glob("*.cleaned.fasta"))
    phases_root = workdir / "phases"
    if phases_root.is_dir():
        for phase_dir in sorted(phases_root.iterdir()):
            pc = phase_dir / "06_clean"
            if pc.is_dir():
                fps.extend(pc.glob("*.cleaned.fasta"))
    fps = sorted(fps, key=lambda p: str(p.resolve()))
    genomes = [str(p.resolve()) for p in fps]
    labels = [p.stem for p in fps]
    (outd / "genome_list.txt").write_text("\n".join(genomes) + ("\n" if genomes else ""))
    (outd / "labels.txt").write_text("\n".join(labels) + ("\n" if labels else ""))


def run_pipeline(cfg_path: Path, dry_run: bool) -> int:
    cfg = load_config(cfg_path)
    workdir = Path(cfg["project"]["workdir"]).expanduser().resolve()
    if not dry_run:
        _mkdir(workdir)
    log_cfg = cfg.get("logging", {})
    if dry_run:
        log_dir = Path.cwd() / ".najdi_pipeline_dryrun_logs"
        _mkdir(log_dir)
    else:
        log_dir = workdir / (log_cfg.get("dir") or "logs")
        _mkdir(log_dir)
    driver_log = log_dir / f"driver_{_ts()}.log"
    failures_tsv = log_dir / "failures.tsv"
    manifest_path = log_dir / "manifest.jsonl"
    tail_n = int(log_cfg.get("failure_tail_lines", 80))

    logging.basicConfig(
        level=getattr(logging, str(log_cfg.get("verbosity", "INFO")).upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(driver_log), logging.StreamHandler(sys.stdout)],
    )
    log = logging.getLogger("najdi")

    jobs = build_jobs(cfg)
    if not jobs:
        log.warning("No jobs built from config.")
        return 0

    state_path = workdir / (cfg.get("resume", {}).get("state_file") or ".pipeline_state.json")
    state = load_state(state_path) if cfg.get("resume", {}).get("enabled", True) else {}

    waves = topological_waves(jobs)
    log.info("Planned %d jobs in %d waves", len(jobs), len(waves))

    ex = cfg.get("execution", {})
    mode = ex.get("mode", "slurm")
    slurm_mod = (ex.get("slurm") or {}).get("modules_by_step", {})

    for wi, wave in enumerate(waves):
        log.info("--- Wave %d: %d jobs ---", wi, len(wave))
        for job in wave:
            if job.outputs and cfg.get("resume", {}).get("enabled", True):
                if outputs_ready(job.outputs) and state.get(job.job_id) == "done":
                    log.info("Resume skip (done): %s", job.job_id)
                    continue
                if outputs_ready(job.outputs):
                    log.info("Resume skip (outputs exist): %s", job.job_id)
                    state[job.job_id] = "done"
                    continue

            if job.job_id == "quast_multi_final":
                write_quast_inputs_files(workdir)

            job_log = log_dir / "jobs" / f"{job.job_id}.log"
            _mkdir(job_log.parent)
            log.info("START %s step=%s", job.job_id, job.step_name)
            append_manifest(
                manifest_path,
                {"job_id": job.job_id, "step": job.step_name, "log": str(job_log), "ts": _ts()},
            )

            modules = slurm_mod.get(job.resources_key) or slurm_mod.get(job.step_name)

            if dry_run:
                log.info("[dry-run] Would run %s with script %s", job.job_id, job.script)
                continue

            slurm_id = ""
            if mode == "slurm":
                body = slurm_wrap(job, cfg, job_log, modules)
                with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".sbatch", delete=False, dir=str(job_log.parent)
                ) as tf:
                    tf.write(body)
                    sbatch_path = tf.name
                try:
                    out = subprocess.check_output(
                        ["sbatch", "--parsable", sbatch_path], text=True
                    ).strip()
                    slurm_id = out.split(";")[0]
                    log.info("Submitted %s SLURM %s", job.job_id, slurm_id)
                    # Wait for job completion (sacct) — simple poll
                    rc = _wait_slurm(slurm_id, log)
                finally:
                    Path(sbatch_path).unlink(missing_ok=True)
            else:
                local = ex.get("local", {})
                use_docker = bool(local.get("use_docker", False))
                image = (local.get("images") or {}).get(job.resources_key) or local.get("images", {}).get("default")
                binds = [str(Path(b).resolve()) for b in local.get("bind_mounts", [])]
                res = resources_for(cfg, job.resources_key)
                rc = run_local(
                    job,
                    job_log,
                    use_docker,
                    image,
                    binds,
                    int(res.get("cpus", 4)),
                    str(res.get("memory", "8G")),
                )

            tail = tail_file(job_log, tail_n)
            if rc != 0:
                log.error("FAILED %s exit=%s log=%s", job.job_id, rc, job_log)
                write_failure_row(
                    failures_tsv,
                    {
                        "timestamp": _ts(),
                        "job_id": job.job_id,
                        "assembly_id": job.env.get("LABEL", job.env.get("TAG", "")),
                        "step": job.step_name,
                        "exit_code": rc,
                        "slurm_job_id": slurm_id,
                        "log_path": str(job_log),
                        "tail_stderr": tail,
                    },
                )
                state[job.job_id] = "failed"
                save_state(state_path, state)
                return rc
            state[job.job_id] = "done"
            save_state(state_path, state)
            log.info("DONE %s", job.job_id)

    log.info("All waves complete.")
    if not dry_run:
        summ = {
            "finished_at": _ts(),
            "total_jobs": len(jobs),
            "state": state,
        }
        sd = workdir / (log_cfg.get("dir") or "logs")
        _mkdir(sd)
        (sd / "run_summary.json").write_text(json.dumps(summ, indent=2))
    return 0


def _wait_slurm(job_id: str, log: logging.Logger, poll: int = 30) -> int:
    """Poll sacct until job completes; return exit code approximation."""
    while True:
        proc = subprocess.run(
            ["sacct", "-j", job_id, "--format=State,ExitCode", "--noheader", "--parsable2"],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            log.warning("sacct failed: %s", proc.stderr)
            time.sleep(poll)
            continue
        lines = [ln.strip() for ln in proc.stdout.strip().splitlines() if ln.strip()]
        if not lines:
            time.sleep(poll)
            continue
        # first line often batch job
        state = lines[0].split("|")[0]
        if state in ("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL"):
            raw = lines[0].split("|")[-1] if "|" in lines[0] else "1"
            try:
                ec = int(raw.split(":")[-1])
            except ValueError:
                ec = 0
            if state == "COMPLETED":
                return ec
            return ec if ec != 0 else 1
        time.sleep(poll)


def main() -> None:
    ap = argparse.ArgumentParser(description="Najdi T2T assembly pipeline driver")
    ap.add_argument("--config", required=True, type=Path, help="Path to YAML config")
    ap.add_argument("--dry-run", action="store_true", help="Print plan only")
    args = ap.parse_args()
    sys.exit(run_pipeline(args.config, args.dry_run))


if __name__ == "__main__":
    main()
