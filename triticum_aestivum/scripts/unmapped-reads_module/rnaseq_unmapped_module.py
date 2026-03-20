#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
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


def extract_unmapped_reads(input_align, input_type, layout, outdir, threads, sample, log_file):
    workdir = os.path.join(outdir, "01_extract_unmapped")
    ensure_dir(workdir)

    if input_type == "sam":
        input_bam = os.path.join(workdir, f"{sample}.input.bam")
        sam_to_bam(input_align, input_bam, threads, log_file)
    else:
        input_bam = input_align

    unmapped_bam = os.path.join(workdir, f"{sample}.unmapped.bam")

    run_cmd(
        [
            "samtools", "view",
            "-@", str(threads),
            "-b",
            "-f", "4",
            input_bam,
            "-o", unmapped_bam
        ],
        log_file=log_file
    )

    if layout == "PE":
        name_sorted_bam = os.path.join(workdir, f"{sample}.unmapped.name_sorted.bam")
        run_cmd(
            [
                "samtools", "sort",
                "-@", str(threads),
                "-n",
                "-o", name_sorted_bam,
                unmapped_bam
            ],
            log_file=log_file
        )

        fq1 = os.path.join(workdir, f"{sample}.unmapped.R1.fastq.gz")
        fq2 = os.path.join(workdir, f"{sample}.unmapped.R2.fastq.gz")
        fqs = os.path.join(workdir, f"{sample}.unmapped.singletons.fastq.gz")

        run_cmd(
            [
                "samtools", "fastq",
                "-@", str(threads),
                "-1", fq1,
                "-2", fq2,
                "-s", fqs,
                "-0", "/dev/null",
                name_sorted_bam
            ],
            log_file=log_file
        )

        return {
            "left": fq1,
            "right": fq2,
            "single": fqs
        }

    else:
        fq = os.path.join(workdir, f"{sample}.unmapped.fastq.gz")
        run_cmd(
            [
                "samtools", "fastq",
                "-@", str(threads),
                "-0", fq,
                unmapped_bam
            ],
            log_file=log_file
        )
        return {
            "single": fq
        }


def run_trinity(layout, fq_dict, outdir, sample, threads, max_memory, min_contig_length, log_file):
    trinity_dir = os.path.join(outdir, "02_trinity", f"{sample}.trinity")
    ensure_dir(os.path.dirname(trinity_dir))

    cmd = [
        "Trinity",
        "--seqType", "fq",
        "--max_memory", max_memory,
        "--CPU", str(threads),
        "--min_contig_length", str(min_contig_length),
        "--output", trinity_dir
    ]

    if layout == "PE":
        if not file_nonempty(fq_dict["left"]) or not file_nonempty(fq_dict["right"]):
            raise RuntimeError("Paired-end mode selected, but R1/R2 unmapped FASTQ is empty or missing.")

        cmd.extend(["--left", fq_dict["left"], "--right", fq_dict["right"]])

        if "single" in fq_dict and file_nonempty(fq_dict["single"]):
            print(f"[INFO] Singleton reads detected: {fq_dict['single']}")
            print("[INFO] Trinity does not allow mixing PE and SE inputs in one run; singleton reads will be ignored in PE assembly.")
            with open(log_file, "a") as log:
                log.write(f"[INFO] Singleton reads detected: {fq_dict['single']}\n")
                log.write("[INFO] Trinity does not allow mixing PE and SE inputs in one run; singleton reads will be ignored in PE assembly.\n")

    else:
        if not file_nonempty(fq_dict["single"]):
            raise RuntimeError("Single-end unmapped FASTQ is empty or missing.")
        cmd.extend(["--single", fq_dict["single"]])

    run_cmd(cmd, log_file=log_file)

    trinity_fasta = os.path.join(trinity_dir, "Trinity.fasta")
    if not file_nonempty(trinity_fasta):
        raise RuntimeError(f"Trinity assembly output not found: {trinity_fasta}")

    final_trinity_fasta = os.path.join(outdir, "02_trinity", f"{sample}.trinity.Trinity.fasta")
    shutil.copyfile(trinity_fasta, final_trinity_fasta)

    if not file_nonempty(final_trinity_fasta):
        raise RuntimeError(f"Trinity assembly output not found: {final_trinity_fasta}")

    return final_trinity_fasta


def build_blast_db_if_needed(fasta, dbtype, log_file):
    if dbtype == "nucl":
        index_file = fasta + ".nin"
    else:
        index_file = fasta + ".pin"

    if not os.path.exists(index_file):
        run_cmd(
            [
                "makeblastdb",
                "-in", fasta,
                "-dbtype", dbtype
            ],
            log_file=log_file
        )


def run_blastn(query_fasta, db_fasta, outdir, sample, threads, evalue, max_targets, log_file):
    ensure_dir(outdir)
    build_blast_db_if_needed(db_fasta, "nucl", log_file)

    out_file = os.path.join(outdir, f"{sample}.blastn.tsv")
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    run_cmd(
        [
            "blastn",
            "-query", query_fasta,
            "-db", db_fasta,
            "-out", out_file,
            "-outfmt", outfmt,
            "-evalue", str(evalue),
            "-max_target_seqs", str(max_targets),
            "-num_threads", str(threads)
        ],
        log_file=log_file
    )
    return out_file


def run_blastx(query_fasta, db_fasta, outdir, sample, threads, evalue, max_targets, log_file):
    ensure_dir(outdir)
    build_blast_db_if_needed(db_fasta, "prot", log_file)

    out_file = os.path.join(outdir, f"{sample}.blastx.tsv")
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    run_cmd(
        [
            "blastx",
            "-query", query_fasta,
            "-db", db_fasta,
            "-out", out_file,
            "-outfmt", outfmt,
            "-evalue", str(evalue),
            "-max_target_seqs", str(max_targets),
            "-num_threads", str(threads)
        ],
        log_file=log_file
    )
    return out_file


def summarize_best_hits(tsv_file, summary_file):
    if not file_nonempty(tsv_file):
        with open(summary_file, "w") as out:
            out.write("qseqid\tsseqid\tpident\tlength\tevalue\tbitscore\n")
        return

    best = {}
    with open(tsv_file) as f:
        for line in f:
            arr = line.rstrip("\n").split("\t")
            if len(arr) < 12:
                continue
            qseqid = arr[0]
            sseqid = arr[1]
            pident = arr[2]
            length = arr[3]
            evalue = float(arr[10])
            bitscore = float(arr[11])

            if qseqid not in best:
                best[qseqid] = (sseqid, pident, length, evalue, bitscore)
            else:
                old = best[qseqid]
                if bitscore > old[4] or (bitscore == old[4] and evalue < old[3]):
                    best[qseqid] = (sseqid, pident, length, evalue, bitscore)

    with open(summary_file, "w") as out:
        out.write("qseqid\tsseqid\tpident\tlength\tevalue\tbitscore\n")
        for qseqid, val in best.items():
            out.write(
                f"{qseqid}\t{val[0]}\t{val[1]}\t{val[2]}\t{val[3]}\t{val[4]}\n"
            )


def parse_args():
    parser = argparse.ArgumentParser(
        description="RNA-seq unmapped-read module: extract unmapped reads, Trinity assembly, and simple annotation."
    )

    parser.add_argument("--input", required=True, help="Input alignment file (.sam or .bam) containing unmapped records.")
    parser.add_argument("--sample", required=True, help="Sample name.")
    parser.add_argument("--layout", required=True, choices=["SE", "PE"], help="Library layout: SE or PE.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads. Default: 8")
    parser.add_argument("--max-memory", default="50G", help="Trinity max memory, e.g. 50G")
    parser.add_argument("--min-contig-length", type=int, default=200, help="Minimum Trinity contig length. Default: 200")
    parser.add_argument("--evalue", type=float, default=1e-5, help="BLAST e-value cutoff. Default: 1e-5")
    parser.add_argument("--max-target-seqs", type=int, default=5, help="BLAST max target seqs. Default: 5")

    parser.add_argument("--cdna-db", default=None, help="Reference transcript FASTA for blastn annotation.")
    parser.add_argument("--protein-db", default=None, help="Reference protein FASTA for blastx annotation.")
    parser.add_argument("--skip-annotation", action="store_true", help="Skip annotation step.")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary intermediate BAM files.")

    return parser.parse_args()


def main():
    args = parse_args()

    ensure_dir(args.outdir)
    log_file = os.path.join(args.outdir, f"{args.sample}.unmapped_module.log")

    check_executable("samtools")
    check_executable("Trinity")

    if not args.skip_annotation:
        if args.cdna_db is not None:
            check_executable("makeblastdb")
            check_executable("blastn")
        if args.protein_db is not None:
            check_executable("makeblastdb")
            check_executable("blastx")

    input_type = infer_input_type(args.input)

    fq_dict = extract_unmapped_reads(
        input_align=args.input,
        input_type=input_type,
        layout=args.layout,
        outdir=args.outdir,
        threads=args.threads,
        sample=args.sample,
        log_file=log_file
    )

    trinity_fasta = run_trinity(
        layout=args.layout,
        fq_dict=fq_dict,
        outdir=args.outdir,
        sample=args.sample,
        threads=args.threads,
        max_memory=args.max_memory,
        min_contig_length=args.min_contig_length,
        log_file=log_file
    )

    ann_dir = os.path.join(args.outdir, "03_annotation", args.sample)
    ensure_dir(ann_dir)

    if args.skip_annotation:
        print("[INFO] Annotation step skipped.")
    else:
        if args.cdna_db:
            blastn_out = run_blastn(
                query_fasta=trinity_fasta,
                db_fasta=args.cdna_db,
                outdir=ann_dir,
                sample=args.sample,
                threads=args.threads,
                evalue=args.evalue,
                max_targets=args.max_target_seqs,
                log_file=log_file
            )
            summarize_best_hits(
                blastn_out,
                os.path.join(ann_dir, f"{args.sample}.blastn.besthit.tsv")
            )

        if args.protein_db:
            blastx_out = run_blastx(
                query_fasta=trinity_fasta,
                db_fasta=args.protein_db,
                outdir=ann_dir,
                sample=args.sample,
                threads=args.threads,
                evalue=args.evalue,
                max_targets=args.max_target_seqs,
                log_file=log_file
            )
            summarize_best_hits(
                blastx_out,
                os.path.join(ann_dir, f"{args.sample}.blastx.besthit.tsv")
            )

    print("[DONE] RNA-seq unmapped-read module finished successfully.")
    print(f"[INFO] Trinity assembly: {trinity_fasta}")
    print(f"[INFO] Log file: {log_file}")


if __name__ == "__main__":
    main()