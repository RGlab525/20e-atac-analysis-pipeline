#!/usr/bin/env python3
"""Top de novo motifs (lowest P) from a HOMER homerResults.html, as a TSV for the bar plots of Fig 3E / 4E.

Usage: python3 Fig3/homer_top_motifs.py <HOMER output dir> <group name> <out.tsv> [top_n=10]
Columns: group, rank, p_value, log_p_value, best_match_gene, best_match_full, target_percent, background_percent
The motif logos in the figures were taken from the HOMER HTML report (SVG) and placed in Illustrator.
Source: 260612_Enh_motif/extract_top_motifs.py (parse_de_novo), reduced to the table output.
"""
import csv
import html
import re
import sys
from pathlib import Path


def clean_text(value):
    value = re.sub(r"<br\s*/?>", " ", value, flags=re.I)
    value = re.sub(r"<[^>]+>", "", value)
    value = html.unescape(value)
    return re.sub(r"\s+", " ", value).strip()


def html_rows(path):
    text = path.read_text(errors="replace")
    starts = [m.start() for m in re.finditer(r"<TR><TD>\d+", text, flags=re.I)]
    rows = []
    for i, start in enumerate(starts):
        end = starts[i + 1] if i + 1 < len(starts) else text.find("</TABLE>", start)
        rows.append(text[start:end if end != -1 else len(text)].strip())
    return rows


def split_cells(row_html):
    row_html = re.sub(r"^<TR><TD>", "", row_html.strip(), flags=re.I)
    row_html = re.sub(r"</TD></TR>$", "", row_html.strip(), flags=re.I)
    return re.split(r"</TD><TD>", row_html, flags=re.I)


def short_name(full_name):
    """Gene part of a HOMER best-match name, e.g. 'Trl(Zf)/S2-GAGA.../Homer(0.72)' -> 'Trl'."""
    split_at, depth = len(full_name), 0
    for i, ch in enumerate(full_name):
        if ch == "(":
            depth += 1
        elif ch == ")" and depth:
            depth -= 1
        elif ch == "/" and depth == 0:
            split_at = i
            break
    name = re.sub(r"\([^)]*\)", "", full_name[:split_at]).strip()
    return name or full_name


def pvalue_key(value):
    try:
        return float(value.replace("E", "e"))
    except ValueError:
        return float("inf")


def parse_de_novo(homer_dir, group):
    records = []
    for row in html_rows(Path(homer_dir) / "homerResults.html"):
        cells = split_cells(row)
        if len(cells) < 9:
            continue
        rank = re.search(r"\d+", clean_text(cells[0]))
        if not rank:
            continue
        best_full = re.sub(r"\s*More Information.*$", "", clean_text(cells[7])).strip()
        records.append({
            "group": group,
            "rank": int(rank.group(0)),
            "p_value": clean_text(cells[2]),
            "log_p_value": clean_text(cells[3]),
            "best_match_gene": short_name(best_full),
            "best_match_full": best_full,
            "target_percent": clean_text(cells[4]).replace("%", ""),
            "background_percent": clean_text(cells[5]).replace("%", ""),
        })
    return sorted(records, key=lambda r: pvalue_key(r["p_value"]))


if __name__ == "__main__":
    homer_dir, group, out = sys.argv[1:4]
    top_n = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    recs = parse_de_novo(homer_dir, group)[:top_n]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(recs[0].keys()))
        w.writeheader()
        w.writerows(recs)
