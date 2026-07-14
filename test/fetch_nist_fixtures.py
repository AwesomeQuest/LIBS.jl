#!/usr/bin/env python3
"""Fetch NIST ASD LIBS reference data for a set of elements."""

import re
import json
import urllib.request
import time
import sys
import os

BASE_URL = (
    "https://physics.nist.gov/cgi-bin/ASD/lines1.pl"
    "?composition={elem}%3A100"
    "&mytext%5B%5D={elem}"
    "&myperc%5B%5D=100"
    "&spectra={elem}{lo}-{hi}"
    "&low_w=200"
    "&limits_type=0"
    "&upp_w=600"
    "&show_av=2"
    "&unit=1"
    "&resolution=1000"
    "&temp=1"
    "&eden=1e17"
    "&maxcharge={maxcharge}"
    "&min_rel_int=0.001"
    "&int_scale=1"
    "&libs=1"
)

ELEMENTS = [
    ("H", 1, 0, 0),
    ("He", 2, 0, 1),
    ("Li", 3, 0, 1),
    ("C", 6, 0, 2),
    ("N", 7, 0, 2),
    ("O", 8, 0, 1),
    ("Al", 13, 0, 2),
    ("Si", 14, 0, 2),
    ("Ar", 18, 0, 1),
    ("Fe", 26, 0, 2),
    ("Ni", 28, 0, 2),
    ("Cu", 29, 0, 1),
    ("U", 92, 0, 2),
    ("Kr", 36, 0, 1),
    ("Zn", 30, 0, 1),
]


def parse_data_sticks_array(html):
    m = re.search(r"var dataSticksArray\s*=\s*(\[.*?\])\s*;", html, re.DOTALL)
    if not m:
        return None
    arr_str = m.group(1)
    arr_str = arr_str[1:]  # remove leading [
    rows = []
    depth = 0
    current = ""
    for ch in arr_str:
        if ch == "[":
            if depth > 0:
                current += ch
            depth += 1
        elif ch == "]":
            depth -= 1
            if depth == 0:
                rows.append(current)
                current = ""
            else:
                current += ch
        else:
            if depth >= 1:
                current += ch
    if not rows:
        return None
    header = rows[0]
    charge_labels = re.findall(r"label:'([^']+)'", header)
    data = []
    for row_str in rows[1:]:
        parts = row_str.split(",")
        if not parts:
            continue
        wl = float(parts[0].strip())
        ints = [
            None if p.strip() == "null" else float(p.strip())
            for p in parts[1:]
        ]
        data.append({"wl": wl, "intensities": ints})
    stages = []
    pops = []
    for label in charge_labels:
        m_stage = re.match(r"(.+?)\s*\(([\d.eE+-]+)\)", label)
        if m_stage and label != "Ritz wavelength (nm)":
            stages.append(m_stage.group(1).strip())
            pops.append(float(m_stage.group(2)))
        else:
            stages.append(label)
            pops.append(None)
    return {"charge_stages": stages, "saha_populations": pops, "data": data}


def main():
    output_dir = os.path.dirname(os.path.abspath(__file__))
    fixtures_dir = os.path.join(output_dir, "fixtures")
    os.makedirs(fixtures_dir, exist_ok=True)

    success = 0
    fail = 0

    for symbol, z, lo, hi in ELEMENTS:
        url = BASE_URL.format(
            elem=symbol, lo=lo, hi=hi, maxcharge=hi
        )
        fixture_path = os.path.join(fixtures_dir, f"{symbol.lower()}.json")
        html_path = os.path.join(fixtures_dir, f"{symbol.lower()}.html")

        print(f"Fetching {symbol} (Z={z})... ", end="", flush=True)

        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": (
                        "Mozilla/5.0 (X11; Linux x86_64) "
                        "AppleWebKit/537.36"
                    )
                },
            )
            with urllib.request.urlopen(req, timeout=30) as resp:
                html = resp.read().decode("utf-8", errors="replace")

            # Save raw HTML as backup
            with open(html_path, "w") as f:
                f.write(html)

            parsed = parse_data_sticks_array(html)
            if parsed is None:
                print("FAIL (no dataSticksArray)")
                fail += 1
                continue

            if not parsed["data"]:
                print("WARN (empty data)")
                # Still save the empty fixture
            else:
                print(
                    f"OK ({len(parsed['data'])} lines, "
                    f"{len([s for s in parsed['charge_stages'] if 'nm' not in s])} stages)"
                )

            parsed["element"] = symbol
            parsed["url_params"] = {
                "temp": 1,
                "eden": 1e17,
                "resolution": 1000,
                "low_w": 200,
                "upp_w": 600,
                "min_rel_int": 0.001,
            }

            with open(fixture_path, "w") as f:
                json.dump(parsed, f, indent=2)

            success += 1
            time.sleep(1.5)  # Be polite to NIST

        except Exception as e:
            print(f"FAIL ({e})")
            fail += 1

    print(f"\nDone: {success} OK, {fail} FAIL")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
