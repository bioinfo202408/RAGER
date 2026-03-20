from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple
import subprocess

import streamlit as st

try:
    import yaml
except Exception:
    yaml = None

try:
    from pdf2image import convert_from_path  # requires poppler
except Exception:
    convert_from_path = None


ROOT = Path(__file__).resolve().parent
MAX_LOG_LINES = 5000
LOG_REFRESH_EVERY = 50


@dataclass(frozen=True)
class Workflow:
    key: str
    label: str
    snakefile: Path
    configfile: Path
    output_dirs: List[Path]


def read_text(p: Path) -> str:
    return p.read_text(encoding="utf-8") if p.exists() else ""


def write_text(p: Path, text: str) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text.rstrip() + "\n", encoding="utf-8")


def validate_yaml(text: str) -> Tuple[bool, str]:
    if yaml is None:
        return True, ""
    try:
        yaml.safe_load(text)
        return True, ""
    except Exception as e:
        return False, str(e)


def workflows() -> List[Workflow]:
    base = ROOT / "scripts" / "snakemake"
    return [
        Workflow(
            "rna",
            "Preprocess RNAseq",
            base / "Preprocess_RNAseq_data" / "RNAseq_snakefile.py",
            base / "Preprocess_RNAseq_data" / "config.yaml",
            [ROOT / "datasets" / "RNAseq"],
        ),
        Workflow(
            "atac",
            "Preprocess ATACseq",
            base / "Preprocess_ATACseq_data" / "ATACseq_snakefile.py",
            base / "Preprocess_ATACseq_data" / "config.yaml",
            [ROOT / "datasets" / "ATACseq"],
        ),
        Workflow(
            "joint",
            "Joint Analysis",
            base / "Joint_analysis" / "Analysis_snakefile.py",
            base / "Joint_analysis" / "config.yaml",
            [
                ROOT / "datasets" / "Promoter_region_analysis",
                ROOT / "datasets" / "Enhancer_region_analysis",
            ],
        ),
        Workflow(
            "custom",
            "Custom genes analysis",
            base / "Custom_genes_analysis" / "Custom_gene_analysis.py",
            base / "Custom_genes_analysis" / "config.yaml",
            [ROOT / "datasets" / "Custom_genes_analysis"],
        ),
    ]


def run_snakemake(wf: Workflow, threads: int) -> bool:
    cmd = [
        "snakemake",
        "--snakefile",
        wf.snakefile.as_posix(),
        "--configfile",
        wf.configfile.as_posix(),
        "-j",
        str(threads),
    ]

    st.info("Analyzing...")
    out_box = st.empty()
    log_lines: List[str] = []
    n = 0

    try:
        proc = subprocess.Popen(
            cmd,
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
    except FileNotFoundError:
        st.error("snakemake was not found in PATH.")
        return False

    with st.spinner("Running..."):
        assert proc.stdout is not None
        for line in proc.stdout:
            n += 1
            log_lines.append(line.rstrip("\n"))
            if len(log_lines) > MAX_LOG_LINES:
                log_lines = log_lines[-MAX_LOG_LINES:]

            if n % LOG_REFRESH_EVERY == 0:
                out_box.code("\n".join(log_lines))

        rc = proc.wait()
        out_box.code("\n".join(log_lines))

    with st.expander("Logs"):
        st.code("\n".join(log_lines) if log_lines else "(no output)")

    if rc == 0:
        st.success("Finished. Please check the output directories.")
        return True

    st.error(f"Run failed (exit code {rc}).")
    return False


def collect_files(dirs: List[Path]) -> List[Path]:
    exts = {".csv", ".txt", ".pdf"}
    files: List[Path] = []
    for d in dirs:
        if d.exists():
            files.extend([p for p in d.rglob("*") if p.is_file() and p.suffix.lower() in exts])
    files.sort(key=lambda p: p.stat().st_mtime if p.exists() else 0, reverse=True)
    return files


def download_button_for_file(p: Path, key: str) -> None:
    st.download_button(
        label="Download",
        data=p.read_bytes(),
        file_name=p.name,
        mime="application/octet-stream",
        key=key,
    )


@st.cache_data(show_spinner=False)
def pdf_to_images(pdf_path: str, mtime: float):
    if convert_from_path is None:
        return None
    return convert_from_path(pdf_path)


def render_results(wf: Workflow) -> None:
    st.subheader("Results")
    for d in wf.output_dirs:
        st.caption(f"Output dir: {d}")

    files = collect_files(wf.output_dirs)
    if not files:
        st.info("No result files found (.csv/.txt/.pdf).")
        return

    tab_files = [p for p in files if p.suffix.lower() in {".csv", ".txt"}]
    pdf_files = [p for p in files if p.suffix.lower() == ".pdf"]

    if tab_files:
        st.markdown("### Tables (.csv / .txt)")
        for i, p in enumerate(tab_files, 1):
            rel = p.relative_to(ROOT) if p.is_relative_to(ROOT) else p
            with st.expander(f"{i}. {rel}", expanded=False):
                c1, c2 = st.columns([1, 5])
                with c1:
                    download_button_for_file(p, key=f"dl_tab_{wf.key}_{i}")
                with c2:
                    if p.suffix.lower() == ".csv":
                        try:
                            import pandas as pd

                            df = pd.read_csv(p)
                            st.dataframe(df, use_container_width=True)
                        except Exception:
                            st.code(p.read_text(encoding="utf-8", errors="ignore")[:20000])
                    else:
                        st.code(p.read_text(encoding="utf-8", errors="ignore")[:20000])

    if pdf_files:
        st.markdown("### PDF")
        if convert_from_path is None:
            st.warning("PDF preview is unavailable (pdf2image/poppler not installed). Downloads are still available.")
        for i, p in enumerate(pdf_files, 1):
            rel = p.relative_to(ROOT) if p.is_relative_to(ROOT) else p
            with st.expander(f"{i}. {rel}", expanded=False):
                c1, c2 = st.columns([1, 5])
                with c1:
                    download_button_for_file(p, key=f"dl_pdf_{wf.key}_{i}")
                with c2:
                    if convert_from_path is None:
                        continue
                    try:
                        imgs = pdf_to_images(p.as_posix(), p.stat().st_mtime)
                        if not imgs:
                            st.info("No pages to preview.")
                        else:
                            for page_idx, img in enumerate(imgs, 1):
                                st.image(img, caption=f"Page {page_idx}", use_container_width=True)
                    except Exception as e:
                        st.warning(f"Failed to render PDF preview: {e}")


def render(wf: Workflow) -> None:
    st.subheader(wf.label)
    st.caption(f"Config: {wf.configfile}")
    st.caption(f"Snakefile: {wf.snakefile}")

    if not wf.snakefile.exists():
        st.error("Snakefile not found.")
        return
    if not wf.configfile.exists():
        st.error("Config file not found.")
        return

    cfg_key = f"cfg_{wf.key}"
    if cfg_key not in st.session_state:
        st.session_state[cfg_key] = read_text(wf.configfile)

    c1, c2, c3 = st.columns([1, 1, 2])
    with c1:
        threads = st.number_input("Threads (-j)", 1, 256, 10, 1, key=f"t_{wf.key}")
    with c2:
        if st.button("Reload", key=f"reload_{wf.key}"):
            st.session_state[cfg_key] = read_text(wf.configfile)
            st.success("Reloaded.")
    with c3:
        if st.button("Save", key=f"save_{wf.key}"):
            ok, err = validate_yaml(st.session_state[cfg_key])
            if not ok:
                st.error(f"Invalid YAML: {err}")
            else:
                write_text(wf.configfile, st.session_state[cfg_key])
                st.success("Saved.")

    cfg_text = st.text_area("Edit YAML", height=420, key=cfg_key)

    st.divider()
    if st.button("Run", key=f"run_{wf.key}"):
        ok, err = validate_yaml(cfg_text)
        if not ok:
            st.error(f"Invalid YAML: {err}")
        else:
            write_text(wf.configfile, cfg_text)
            run_snakemake(wf, int(threads))

    st.divider()
    render_results(wf)


def main() -> None:
    st.set_page_config(page_title="RAGER Snakemake Pipeline UI", layout="wide")
    st.title("RAGER Snakemake Pipeline UI")

    wfs = workflows()
    tabs = st.tabs([w.label for w in wfs])
    for tab, wf in zip(tabs, wfs):
        with tab:
            render(wf)


if __name__ == "__main__":
    main()
