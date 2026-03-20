#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import argparse
import subprocess
import shutil
from pathlib import Path


def run_cmd(cmd, log_file=None, cwd=None):
    cmd_str = " ".join(cmd)
    print(f"[RUN] {cmd_str}")
    if log_file:
        with open(log_file, "a") as log:
            log.write(f"\n[RUN] {cmd_str}\n")
            process = subprocess.run(cmd, stdout=log, stderr=log, cwd=cwd)
    else:
        process = subprocess.run(cmd, cwd=cwd)

    if process.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd_str}")


def get_cmd_output(cmd):
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{result.stderr}")
    return result.stdout.strip()


def check_executable(name):
    if shutil.which(name) is None:
        raise FileNotFoundError(f"Required executable not found in PATH: {name}")


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def file_nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def infer_input_type(input_file):
    lower = input_file.lower()
    if lower.endswith(".sam"):
        return "sam"
    elif lower.endswith(".bam"):
        return "bam"
    else:
        raise ValueError("Input file must end with .sam or .bam")


def sam_to_bam(input_sam, output_bam, threads, log_file):
    run_cmd(
        [
            "samtools", "view",
            "-@", str(threads),
            "-bS", input_sam,
            "-o", output_bam
        ],
        log_file=log_file
    )


def count_fasta_records(fasta_file):
    if not file_nonempty(fasta_file):
        return 0

    count = 0
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def count_mapped_records_from_sam(sam_file, threads):
    if not file_nonempty(sam_file):
        return 0

    out = get_cmd_output(
        [
            "samtools", "view",
            "-@", str(threads),
            "-c",
            "-F", "4",
            sam_file
        ]
    )
    return int(out)


def make_fasta_headers_unique(input_fasta, output_fasta):
    seen = {}
    with open(input_fasta, "r") as fin, open(output_fasta, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                if not header:
                    seq_id = "seq"
                    desc = ""
                else:
                    parts = header.split(None, 1)
                    seq_id = parts[0]
                    desc = parts[1] if len(parts) > 1 else ""

                seen[seq_id] = seen.get(seq_id, 0) + 1
                new_id = f"{seq_id}__{seen[seq_id]}"

                if desc:
                    fout.write(f">{new_id} {desc}\n")
                else:
                    fout.write(f">{new_id}\n")
            else:
                fout.write(line)


def remove_bowtie2_index_files(index_base):
    parent = Path(index_base).parent
    stem = Path(index_base).name
    for p in parent.glob(f"{stem}*.bt2*"):
        if p.exists():
            p.unlink()


def build_bowtie2_index_if_needed(ref_fasta, index_base, log_file):
    bt2_small = index_base + ".1.bt2"
    bt2_large = index_base + ".1.bt2l"
    unique_ref_fasta = index_base + ".unique_headers.fa"

    index_exists = os.path.exists(bt2_small) or os.path.exists(bt2_large)

    if not os.path.exists(unique_ref_fasta):
        make_fasta_headers_unique(ref_fasta, unique_ref_fasta)
        if index_exists:
            remove_bowtie2_index_files(index_base)
            index_exists = False

    if not index_exists:
        run_cmd(
            [
                "bowtie2-build",
                unique_ref_fasta,
                index_base
            ],
            log_file=log_file
        )

    return unique_ref_fasta


def extract_unmapped_to_fasta(input_align, input_type, outdir, threads, sample, log_file):
    """
    Extract unmapped records from SAM/BAM and export them to FASTA.
    Secondary and supplementary alignments are excluded.
    """
    workdir = os.path.join(outdir, "01_extract_unmapped")
    ensure_dir(workdir)

    if input_type == "sam":
        input_bam = os.path.join(workdir, f"{sample}.input.bam")
        sam_to_bam(input_align, input_bam, threads, log_file)
    else:
        input_bam = input_align

    unmapped_bam = os.path.join(workdir, f"{sample}.unmapped.bam")
    unmapped_fa = os.path.join(workdir, f"{sample}.unmapped.fa")

    run_cmd(
        [
            "samtools", "view",
            "-@", str(threads),
            "-b",
            "-f", "4",
            "-F", "2304",
            input_bam,
            "-o", unmapped_bam
        ],
        log_file=log_file
    )

    with open(unmapped_fa, "w") as outfa:
        process = subprocess.run(
            ["samtools", "fasta", "-@", str(threads), unmapped_bam],
            stdout=outfa,
            stderr=subprocess.PIPE,
            text=True
        )
    if process.returncode != 0:
        raise RuntimeError(f"Command failed: samtools fasta -@ {threads} {unmapped_bam}\n{process.stderr}")

    return unmapped_bam, unmapped_fa


def align_and_collect_unaligned(
    input_fasta,
    ref_fasta,
    ref_name,
    outdir,
    sample,
    threads,
    log_file,
    sensitivity="--very-sensitive"
):
    """
    Align input FASTA reads to a reference using bowtie2 in single-end mode.
    Return:
      sam_out, unaligned_fasta, mapped_count
    """
    ensure_dir(outdir)
    index_dir = os.path.join(outdir, "indexes")
    ensure_dir(index_dir)

    index_base = os.path.join(index_dir, f"{ref_name}")
    build_bowtie2_index_if_needed(ref_fasta, index_base, log_file)

    sam_out = os.path.join(outdir, f"{sample}.vs_{ref_name}.sam")
    unaligned_fa = os.path.join(outdir, f"{sample}.after_{ref_name}.unaligned.fa")

    if not file_nonempty(input_fasta):
        with open(unaligned_fa, "w") as f:
            pass
        return sam_out, unaligned_fa, 0

    run_cmd(
        [
            "bowtie2",
            sensitivity,
            "-p", str(threads),
            "-f",
            "-x", index_base,
            "-U", input_fasta,
            "--un", unaligned_fa,
            "-S", sam_out
        ],
        log_file=log_file
    )

    mapped_count = count_mapped_records_from_sam(sam_out, threads)
    return sam_out, unaligned_fa, mapped_count


def write_summary(summary_file, sample, total_unmapped, cp_count, mt_count, te_count, unclassified_count):
    def safe_pct(x, total):
        return 0.0 if total == 0 else (x / total) * 100.0

    rows = [
        {
            "sample": sample,
            "category": "chloroplast",
            "read_count": cp_count,
            "percentage": f"{safe_pct(cp_count, total_unmapped):.6f}"
        },
        {
            "sample": sample,
            "category": "mitochondrial",
            "read_count": mt_count,
            "percentage": f"{safe_pct(mt_count, total_unmapped):.6f}"
        },
        {
            "sample": sample,
            "category": "repeat_TE",
            "read_count": te_count,
            "percentage": f"{safe_pct(te_count, total_unmapped):.6f}"
        },
        {
            "sample": sample,
            "category": "unclassified",
            "read_count": unclassified_count,
            "percentage": f"{safe_pct(unclassified_count, total_unmapped):.6f}"
        }
    ]

    with open(summary_file, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["sample", "category", "read_count", "percentage"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_wide_summary(summary_file, sample, total_unmapped, cp_count, mt_count, te_count, unclassified_count):
    def safe_pct(x, total):
        return 0.0 if total == 0 else (x / total) * 100.0

    with open(summary_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "sample",
            "total_unmapped_reads",
            "chloroplast_reads",
            "chloroplast_pct",
            "mitochondrial_reads",
            "mitochondrial_pct",
            "repeat_TE_reads",
            "repeat_TE_pct",
            "unclassified_reads",
            "unclassified_pct"
        ])
        writer.writerow([
            sample,
            total_unmapped,
            cp_count,
            f"{safe_pct(cp_count, total_unmapped):.6f}",
            mt_count,
            f"{safe_pct(mt_count, total_unmapped):.6f}",
            te_count,
            f"{safe_pct(te_count, total_unmapped):.6f}",
            unclassified_count,
            f"{safe_pct(unclassified_count, total_unmapped):.6f}"
        ])


def cleanup_temp_files(paths):
    for p in paths:
        if p and os.path.exists(p):
            if os.path.isdir(p):
                shutil.rmtree(p)
            else:
                os.remove(p)


def parse_args():
    parser = argparse.ArgumentParser(
        description="ATAC-seq unmapped-read module: extract unmapped reads, classify them against chloroplast, mitochondrial, and repeat/TE references, and summarize counts and percentages."
    )

    parser.add_argument("--input", required=True, help="Input alignment file (.sam or .bam) containing unmapped records.")
    parser.add_argument("--sample", required=True, help="Sample name.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads. Default: 8")

    parser.add_argument("--chloroplast-ref", required=True, help="Chloroplast genome FASTA.")
    parser.add_argument("--mitochondria-ref", required=True, help="Mitochondrial genome FASTA.")
    parser.add_argument("--te-ref", required=True, help="Repeat/TE reference FASTA.")

    parser.add_argument(
        "--bowtie2-sensitivity",
        default="--very-sensitive",
        choices=["--very-fast", "--fast", "--sensitive", "--very-sensitive", "--very-sensitive-local"],
        help="Bowtie2 sensitivity preset. Default: --very-sensitive"
    )

    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate files.")
    return parser.parse_args()


def main():
    args = parse_args()

    ensure_dir(args.outdir)
    log_file = os.path.join(args.outdir, f"{args.sample}.atac_unmapped_module.log")

    check_executable("samtools")
    check_executable("bowtie2")
    check_executable("bowtie2-build")

    input_type = infer_input_type(args.input)

    unmapped_bam, unmapped_fa = extract_unmapped_to_fasta(
        input_align=args.input,
        input_type=input_type,
        outdir=args.outdir,
        threads=args.threads,
        sample=args.sample,
        log_file=log_file
    )

    total_unmapped = count_fasta_records(unmapped_fa)

    summary_dir = os.path.join(args.outdir, "03_summary")
    ensure_dir(summary_dir)

    if total_unmapped == 0:
        long_summary = os.path.join(summary_dir, f"{args.sample}.unmapped_classification.long.tsv")
        wide_summary = os.path.join(summary_dir, f"{args.sample}.unmapped_classification.wide.tsv")
        write_summary(long_summary, args.sample, 0, 0, 0, 0, 0)
        write_wide_summary(wide_summary, args.sample, 0, 0, 0, 0, 0)

        print("[INFO] No unmapped reads found.")
        print(f"[INFO] Long summary: {long_summary}")
        print(f"[INFO] Wide summary: {wide_summary}")
        print(f"[INFO] Log file: {log_file}")
        return

    classify_dir = os.path.join(args.outdir, "02_classification")
    ensure_dir(classify_dir)

    cp_dir = os.path.join(classify_dir, "chloroplast")
    cp_sam, after_cp_fa, cp_count = align_and_collect_unaligned(
        input_fasta=unmapped_fa,
        ref_fasta=args.chloroplast_ref,
        ref_name="chloroplast",
        outdir=cp_dir,
        sample=args.sample,
        threads=args.threads,
        log_file=log_file,
        sensitivity=args.bowtie2_sensitivity
    )

    mt_dir = os.path.join(classify_dir, "mitochondrial")
    mt_sam, after_mt_fa, mt_count = align_and_collect_unaligned(
        input_fasta=after_cp_fa,
        ref_fasta=args.mitochondria_ref,
        ref_name="mitochondrial",
        outdir=mt_dir,
        sample=args.sample,
        threads=args.threads,
        log_file=log_file,
        sensitivity=args.bowtie2_sensitivity
    )

    te_dir = os.path.join(classify_dir, "repeat_TE")
    te_sam, after_te_fa, te_count = align_and_collect_unaligned(
        input_fasta=after_mt_fa,
        ref_fasta=args.te_ref,
        ref_name="repeat_TE",
        outdir=te_dir,
        sample=args.sample,
        threads=args.threads,
        log_file=log_file,
        sensitivity=args.bowtie2_sensitivity
    )

    unclassified_count = count_fasta_records(after_te_fa)

    assigned_total = cp_count + mt_count + te_count + unclassified_count
    if assigned_total != total_unmapped:
        diff = total_unmapped - assigned_total
        unclassified_count += diff

    long_summary = os.path.join(summary_dir, f"{args.sample}.unmapped_classification.long.tsv")
    wide_summary = os.path.join(summary_dir, f"{args.sample}.unmapped_classification.wide.tsv")

    write_summary(
        long_summary,
        args.sample,
        total_unmapped,
        cp_count,
        mt_count,
        te_count,
        unclassified_count
    )

    write_wide_summary(
        wide_summary,
        args.sample,
        total_unmapped,
        cp_count,
        mt_count,
        te_count,
        unclassified_count
    )

    if not args.keep_temp:
        temp_to_remove = [
            unmapped_bam,
            cp_sam,
            mt_sam,
            te_sam
        ]
        cleanup_temp_files(temp_to_remove)

    print("[DONE] ATAC-seq unmapped-read module finished successfully.")
    print(f"[INFO] Total unmapped reads: {total_unmapped}")
    print(f"[INFO] Chloroplast reads: {cp_count}")
    print(f"[INFO] Mitochondrial reads: {mt_count}")
    print(f"[INFO] Repeat/TE reads: {te_count}")
    print(f"[INFO] Unclassified reads: {unclassified_count}")
    print(f"[INFO] Long summary: {long_summary}")
    print(f"[INFO] Wide summary: {wide_summary}")
    print(f"[INFO] Log file: {log_file}")


if __name__ == "__main__":
    main()