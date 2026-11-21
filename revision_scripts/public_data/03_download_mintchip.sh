#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Downloads the fold-change bigWig and the narrowPeak file for the six
# Mint-ChIP-seq experiments (POLR2A-degron HCT116):
#   H3K27ac, H3K4me3, H3K9me3 x (1 uM 5-Ph-IAA 6 h, DMSO)
# For each experiment the released GRCh38 files are read from ENCODE;
# the download is checked against ENCODE's md5sum before it is renamed.
# Experiments whose target or condition does not match are skipped. The DMSO
# experiments have no treatment recorded on ENCODE; they are accepted only when
# they share an ENCODE series with the matching IAA experiment.
#
# Output (data/):
#   <mark>_IAA.bw   <mark>_DMSO.bw
#   <mark>_IAA.narrowPeak   <mark>_DMSO.narrowPeak
#   download_log.tsv
#
# Run live in a terminal (progress is printed as each file downloads):
#   bash /projects/b1042/BackmanLab/Lucas/loop_chromatin_regression/scripts/03_download_mintchip.sh
##################################################################

module load python

BASE=/projects/b1042/BackmanLab/Lucas/loop_chromatin_regression
DATA=${BASE}/data
mkdir -p "${DATA}"

python3 -u - "${DATA}" <<'PY'
import json, sys, os, gzip, hashlib, subprocess, collections

data_dir = sys.argv[1]
BASE_URL = "https://www.encodeproject.org"

EXPERIMENTS = [
    ("H3K27ac", "IAA",  "ENCSR728KQX"),
    ("H3K27ac", "DMSO", "ENCSR930EMO"),
    ("H3K4me3", "IAA",  "ENCSR523IAZ"),
    ("H3K4me3", "DMSO", "ENCSR280DQH"),
    ("H3K9me3", "IAA",  "ENCSR906ZFR"),
    ("H3K9me3", "DMSO", "ENCSR649DZS"),
]
# Files seen on the six ENCODE experiment pages (used only as a cross-check)
SEEN_BW = {("H3K27ac", "IAA"): "ENCFF888OJS", ("H3K27ac", "DMSO"): "ENCFF206GDZ",
           ("H3K4me3", "IAA"): "ENCFF921VVL", ("H3K4me3", "DMSO"): "ENCFF064VUP",
           ("H3K9me3", "IAA"): "ENCFF070SYQ", ("H3K9me3", "DMSO"): "ENCFF392QHG"}
SEEN_NP = {("H3K27ac", "IAA"): "ENCFF282FCC", ("H3K27ac", "DMSO"): "ENCFF805FRN",
           ("H3K4me3", "IAA"): "ENCFF752JSF", ("H3K4me3", "DMSO"): "ENCFF636MAK",
           ("H3K9me3", "IAA"): "ENCFF555MYR", ("H3K9me3", "DMSO"): "ENCFF261TUC"}
iaa_series = {}

problems = []
log_rows = []
step = 0
total_steps = len(EXPERIMENTS) * 2


def run_curl(url, dest, max_s, show=False):
    if show:
        r = subprocess.run(["curl", "-L", "--progress-bar", "-m", str(max_s), "-o", dest, url])
        return r.returncode == 0, "curl exit code %d" % r.returncode
    r = subprocess.run(["curl", "-sSL", "-m", str(max_s), "-o", dest, url],
                       capture_output=True, text=True)
    return r.returncode == 0, r.stderr.strip()


def get_experiment(acc):
    dest = os.path.join(data_dir, "." + acc + ".json")
    ok, err = run_curl("%s/experiments/%s/?format=json&frame=embedded" % (BASE_URL, acc), dest, 180)
    if not ok:
        return None, "curl failed: " + err
    try:
        with open(dest) as fh:
            d = json.load(fh)
    except Exception as e:
        return None, "not valid JSON (%s)" % e
    finally:
        if os.path.exists(dest):
            os.remove(dest)
    return d, None


def treat_text(rep):
    bio = (rep.get("library") or {}).get("biosample") or {}
    parts = []
    for t in bio.get("treatments") or []:
        if isinstance(t, dict):
            bits = [str(t.get("treatment_term_name") or "")]
            if t.get("amount") is not None:
                bits.append("%s %s" % (t.get("amount"), t.get("amount_units") or ""))
            if t.get("duration") is not None:
                bits.append("%s %s" % (t.get("duration"), t.get("duration_units") or ""))
            parts.append(" ".join(b.strip() for b in bits if b.strip()))
    return "; ".join(parts)


def cond_of(text):
    t = text.lower()
    has_dmso = "dmso" in t or "dimethyl" in t
    has_iaa = "indole" in t or "iaa" in t or "auxin" in t
    if has_dmso and has_iaa:
        return "IAA+DMSO"
    if has_dmso:
        return "DMSO"
    if has_iaa:
        return "IAA"
    return "OTHER"


def candidates(files, kind):
    c = []
    for f in files:
        if not isinstance(f, dict):
            continue
        if f.get("status") != "released" or f.get("assembly") != "GRCh38":
            continue
        if kind == "bw":
            if f.get("file_format") == "bigWig" and "fold change" in (f.get("output_type") or "").lower():
                c.append(f)
        else:
            if f.get("file_format") == "bed" and f.get("file_format_type") == "narrowPeak":
                c.append(f)
    if kind == "np" and len(c) > 1:
        pref = [f for f in c if f.get("preferred_default")]
        if len(pref) == 1:
            c = pref
    return c


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_verified(f, final_path, kind):
    """Download, check md5 against ENCODE, check the file type, then move into place."""
    if os.path.exists(final_path):
        return "exists"
    part = final_path + ".part"
    ok, err = run_curl(BASE_URL + f["href"], part, 3600, show=True)
    if not ok:
        return "FAILED: curl: " + err
    want = f.get("md5sum")
    got = md5_of(part)
    if want and got != want:
        os.remove(part)
        return "FAILED: md5 %s != ENCODE %s" % (got, want)
    if not want:
        problems.append("%s: ENCODE metadata has no md5sum; md5 not verified" % f.get("accession"))
    if kind == "bw":
        with open(part, "rb") as fh:
            magic = fh.read(4)
        if magic != bytes([0x26, 0xFC, 0x8F, 0x88]):
            os.remove(part)
            return "FAILED: not a bigWig (magic bytes %r)" % magic
        os.rename(part, final_path)
    else:
        try:
            with gzip.open(part, "rt") as fin:
                lines = fin.read().splitlines()
        except Exception as e:
            os.remove(part)
            return "FAILED: not a gzipped bed (%s)" % e
        os.remove(part)
        lines = [ln for ln in lines if ln.strip()]
        if not lines or len(lines[0].split("\t")) < 10:
            return "FAILED: narrowPeak has %d lines, first line has %d columns (need 10)" % (
                len(lines), len(lines[0].split("\t")) if lines else 0)
        with open(final_path, "w") as fout:
            fout.write("\n".join(lines) + "\n")
    return "downloaded"


for mark, cond, acc in EXPERIMENTS:
    print("\n== %s %s %s" % (mark, cond, acc))
    d, err = get_experiment(acc)
    if err:
        problems.append("%s %s %s: %s" % (mark, cond, acc, err))
        continue
    tgt = d.get("target")
    tlabel = (tgt.get("label") or tgt.get("name") or "") if isinstance(tgt, dict) else str(tgt or "")
    texts = sorted({treat_text(r) for r in d.get("replicates") or [] if isinstance(r, dict)})
    conds = sorted({cond_of(t) for t in texts})
    print("   target: %s | treatment: %s" % (tlabel, " | ".join(texts)))
    if tlabel.lower() != mark.lower():
        problems.append("%s %s %s: target is '%s'; skipped" % (mark, cond, acc, tlabel))
        continue
    series_accs = set()
    for s_ in d.get("related_series") or []:
        series_accs.add(s_.get("accession") if isinstance(s_, dict) else str(s_).strip("/").split("/")[-1])
    if cond == "IAA":
        iaa_series[mark] = series_accs
        cond_ok = (conds == ["IAA"])
        why = "treatment on ENCODE is %s (%s)" % (conds, " | ".join(texts))
    else:
        # ENCODE records no treatment for the DMSO vehicle control, so the control is
        # identified by: no treatment listed, and sharing a series with the IAA experiment
        shared = series_accs & iaa_series.get(mark, set())
        cond_ok = (texts == [""] and bool(shared))
        why = "treatments %s, shared series with the IAA experiment: %s" % (texts, sorted(shared) or "none")
        if cond_ok:
            print("   vehicle control: no treatment recorded on ENCODE; same series as the IAA experiment: %s"
                  % ", ".join(sorted(shared)))
    if not cond_ok:
        problems.append("%s %s %s: %s, expected %s; skipped" % (mark, cond, acc, why, cond))
        continue
    files = d.get("files") or []
    for kind, ext in (("bw", "bw"), ("np", "narrowPeak")):
        c = candidates(files, kind)
        label = "bigWig" if kind == "bw" else "narrowPeak"
        if len(c) != 1:
            problems.append("%s %s %s: %d candidate %s files (need exactly 1): %s"
                            % (mark, cond, acc, len(c), label,
                               ", ".join("%s[%s]" % (f.get("accession"), f.get("output_type")) for f in c)))
            continue
        f = c[0]
        seen = (SEEN_BW if kind == "bw" else SEEN_NP).get((mark, cond))
        if seen and f.get("accession") != seen:
            problems.append("%s %s %s: picked %s but the ENCODE page showed %s" % (mark, cond, label, f.get("accession"), seen))
        name = "%s_%s.%s" % (mark, cond, ext)
        final_path = os.path.join(data_dir, name)
        step += 1
        size_mb = ("%.0f MB" % (f["file_size"] / 1e6)) if f.get("file_size") else "size unknown"
        print("   [%d/%d] %s  %s [%s]  %s" % (step, total_steps, name, f.get("accession"),
                                           f.get("output_type"), size_mb))
        status = fetch_verified(f, final_path, kind)
        print("         -> %s" % status)
        if status.startswith("FAILED"):
            problems.append("%s: %s" % (name, status))
        log_rows.append((name, mark, cond, acc, f.get("accession"), f.get("output_type"),
                         f.get("md5sum", ""), status))

with open(os.path.join(data_dir, "download_log.tsv"), "w") as out:
    out.write("\t".join(["file", "mark", "condition", "experiment", "encode_file",
                         "output_type", "encode_md5sum", "status"]) + "\n")
    for r in log_rows:
        out.write("\t".join(str(x) for x in r) + "\n")

print("\nSelf-check: the 12 expected files in %s" % data_dir)
for mark, _, _ in EXPERIMENTS[::2]:
    for cond in ("DMSO", "IAA"):
        for ext in ("bw", "narrowPeak"):
            p = os.path.join(data_dir, "%s_%s.%s" % (mark, cond, ext))
            if os.path.exists(p) and os.path.getsize(p) > 0:
                extra = ""
                if ext == "narrowPeak":
                    with open(p) as fh:
                        extra = "  %d peaks" % sum(1 for _ in fh)
                print("  OK       %-24s %10.1f MB%s" % (os.path.basename(p), os.path.getsize(p) / 1e6, extra))
            else:
                print("  MISSING  %s" % os.path.basename(p))
if problems:
    print("\nPROBLEMS:")
    for p in problems:
        print("  " + p)
PY
