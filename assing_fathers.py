import os
from typing import List, Tuple, Dict, Any, Optional

import pandas as pd


# Configuration
CHILD_DB = r"C:\Users\user\Desktop\genetic\zrya_processed\genotypes_unified.csv"
BULLS_DB = r"C:\Users\user\Desktop\genetic\zrya_processed\fathers_registry.csv"
OUTPUT_DB = r"C:\Users\user\Desktop\genetic\zrya_processed\lokus_database_with_fathers.csv"

# Matching thresholds
MIN_MATCHED_LOCI = 11
MAX_MUTATIONS = 1


def normalize_allele(value: Any) -> str:
    if pd.isna(value):
        return ""
    s = str(value).strip()
    if s == "-" or s == ".":
        return ""
    # unify comma/dot separators, spaces
    s = s.replace(",", ".").replace(" ", "")
    return s


def get_child_loci_pairs(columns: List[str]) -> List[Tuple[str, str, str]]:
    """Return (locus, col1, col2) for child columns without parent suffixes, in canonical order."""
    pairs: List[Tuple[str, str, str]] = []
    s = set(columns)
    for col in columns:
        if col.startswith("1_") and "_otca" not in col and "_materi" not in col:
            locus = col[2:]
            c2 = f"2_{locus}"
            if c2 in s:
                pairs.append((locus, col, c2))
    return pairs


def get_father_id_column(df_bulls: pd.DataFrame) -> str:
    candidates = [
        "reganimal",
        "regotca",
        "bull_id",
        "id",
        "ID",
        "Идентификационный номер",
        "Номер",
        "Number",
    ]
    for c in candidates:
        if c in df_bulls.columns:
            return c
    # fallback: first column
    return df_bulls.columns[0]


def build_signature_counts_for_bulls(df_bulls: pd.DataFrame, child_pairs: List[Tuple[str, str, str]]) -> Dict[int, Dict[str, Tuple[str, str]]]:
    """Map bulls index -> locus -> (a1, a2) using child locus names (1_/2_)."""
    bulls: Dict[int, Dict[str, Tuple[str, str]]] = {}
    for idx, row in df_bulls.iterrows():
        per_locus: Dict[str, Tuple[str, str]] = {}
        for locus, c1, c2 in child_pairs:
            b1 = normalize_allele(row.get(c1, ""))
            b2 = normalize_allele(row.get(c2, ""))
            per_locus[locus] = (b1, b2)
        bulls[idx] = per_locus
    return bulls


def evaluate_match(child_vals: Dict[str, Tuple[str, str]], father_vals: Dict[str, Tuple[str, str]]) -> Tuple[int, int, int]:
    """Return (matches, mismatches, compared) given locus -> (c1,c2) and (f1,f2)."""
    matches = 0
    mismatches = 0
    compared = 0
    for locus, (c1, c2) in child_vals.items():
        f1, f2 = father_vals.get(locus, ("", ""))
        if not (c1 or c2) or not (f1 or f2):
            continue
        compared += 1
        child_set = {x for x in [c1, c2] if x}
        father_set = {x for x in [f1, f2] if x}
        if child_set.intersection(father_set):
            matches += 1
        else:
            mismatches += 1
    return matches, mismatches, compared


def main():
    df_children = pd.read_csv(CHILD_DB, sep=";", dtype=str).fillna("")
    df_bulls = pd.read_csv(BULLS_DB, sep=";", dtype=str).fillna("")

    child_pairs = get_child_loci_pairs(list(df_children.columns))
    if not child_pairs:
        raise RuntimeError("Не удалось определить список локусов у детей (1_/2_ столбцы)")

    father_id_col = get_father_id_column(df_bulls)

    # Build per-row dicts of locus->(a1,a2)
    bulls_loci = build_signature_counts_for_bulls(df_bulls, child_pairs)

    # Prepare child loci map
    children_loci: Dict[int, Dict[str, Tuple[str, str]]] = {}
    for idx, row in df_children.iterrows():
        per_locus: Dict[str, Tuple[str, str]] = {}
        for locus, c1, c2 in child_pairs:
            a1 = normalize_allele(row.get(c1, ""))
            a2 = normalize_allele(row.get(c2, ""))
            per_locus[locus] = (a1, a2)
        children_loci[idx] = per_locus

    # Select children without father assigned
    mask_no_father = df_children["regotca"].astype(str).fillna("").str.strip() == ""
    candidate_children_idx = [i for i in df_children.index.tolist() if mask_no_father.iloc[i]]
    pre_assigned_children = set(df_children.index.tolist()) - set(candidate_children_idx)

    # First pass: find best bull per child under thresholds and collect ALL candidates per child
    best_candidate_for_child: Dict[int, Optional[int]] = {}
    score_for_pair: Dict[Tuple[int, int], Tuple[int, int, int]] = {}
    candidates_for_child: Dict[int, List[Tuple[int, Tuple[int, int, int]]]] = {}

    # Keep diagnostics of best overall candidate per child (even if under thresholds)
    best_overall_for_child: Dict[int, Tuple[Optional[int], Tuple[int, int, int]]] = {}

    for ci in candidate_children_idx:
        cvals = children_loci[ci]

        # Only proceed if child has sufficient filled loci
        child_compared_possible = sum(1 for locus, (a1, a2) in cvals.items() if a1 or a2)
        if child_compared_possible < MIN_MATCHED_LOCI:
            best_candidate_for_child[ci] = None
            continue

        best_idx: Optional[int] = None
        best_tuple: Tuple[int, int, int] = (-1, 999, -1)  # matches desc, mismatches asc, compared desc
        overall_best_idx: Optional[int] = None
        overall_best_tuple: Tuple[int, int, int] = (-1, 999, -1)
        child_candidates: List[Tuple[int, Tuple[int, int, int]]] = []

        for bi in df_bulls.index:
            bvals = bulls_loci[bi]
            matches, mismatches, compared = evaluate_match(cvals, bvals)
            score_for_pair[(ci, bi)] = (matches, mismatches, compared)

            # Track overall best regardless of thresholds
            candidate_tuple_any = (matches, -mismatches, compared)
            overall_best_cmp = (overall_best_tuple[0], -overall_best_tuple[1], overall_best_tuple[2])
            if candidate_tuple_any > overall_best_cmp:
                overall_best_tuple = (matches, mismatches, compared)
                overall_best_idx = bi

            if matches < MIN_MATCHED_LOCI or mismatches > MAX_MUTATIONS:
                continue

            # tie-break: higher matches, then lower mismatches, then higher compared
            candidate_tuple = (matches, -mismatches, compared)
            best_tuple_cmp = (best_tuple[0], -best_tuple[1], best_tuple[2])
            if candidate_tuple > best_tuple_cmp:
                best_tuple = (matches, mismatches, compared)
                best_idx = bi

            # collect as valid candidate
            child_candidates.append((bi, (matches, mismatches, compared)))

        # sort candidates: matches desc, mismatches asc, compared desc
        child_candidates.sort(key=lambda x: (x[1][0], -x[1][1], x[1][2]), reverse=True)
        best_candidate_for_child[ci] = best_idx
        best_overall_for_child[ci] = (overall_best_idx, overall_best_tuple)
        candidates_for_child[ci] = child_candidates

    # Diagnostics
    total_candidates = len(candidate_children_idx)
    found_any = sum(1 for v in best_candidate_for_child.values() if v is not None)
    print(f"Кандидатов-детей без отца: {total_candidates}; найдено сопоставлений: {found_any}")

    # Build mapping bull -> children without enforcing minimum children rule
    bull_to_children: Dict[int, List[int]] = {}
    for ci, bi in best_candidate_for_child.items():
        if bi is None:
            continue
        bull_to_children.setdefault(bi, []).append(ci)

    # Ensure father allele columns exist
    for locus, _, _ in child_pairs:
        for suf in ["1", "2"]:
            col = f"{suf}_{locus}_otca"
            if col not in df_children.columns:
                df_children[col] = ""

    # Update regotca and father loci
    for bi, kids in bull_to_children.items():
        father_id = str(df_bulls.at[bi, father_id_col]).strip()
        bvals = bulls_loci[bi]
        for ci in kids:
            df_children.at[ci, "regotca"] = father_id
            for locus, (f1, f2) in bvals.items():
                df_children.at[ci, f"1_{locus}_otca"] = f1
                df_children.at[ci, f"2_{locus}_otca"] = f2

    # Save updated children CSV
    os.makedirs(os.path.dirname(OUTPUT_DB), exist_ok=True)
    df_children.to_csv(OUTPUT_DB, sep=";", index=False, encoding="utf-8-sig")

    # Print fathers and number of children
    print("Подтвержденные отцы и число потомков:")
    counts = (
        df_children[df_children["regotca"].astype(str).str.strip() != ""]
        .groupby("regotca").size().sort_values(ascending=False)
    )
    for reg, cnt in counts.items():
        print(f"{reg};{cnt}")

    # Build human-readable Excel report with loci horizontally and exactly two rows per pair (child, father)
    # Report will include ALL candidates per child (children without pre-assigned father only),
    # so user can choose among multiple suitable fathers.
    loci_order = [locus for locus, _, _ in child_pairs]

    # Build columns: meta + per-locus two columns
    meta_cols = ["reganimal", "father", "role"]
    locus_cols: List[str] = []
    for locus in loci_order:
        locus_cols.append(f"{locus}_1")
        locus_cols.append(f"{locus}_2")
    all_cols = meta_cols + locus_cols

    report_rows: List[Dict[str, Any]] = []

    # Iterate children without pre-assigned father; for each - include ALL candidate fathers
    total_pairs = 0
    for ci in candidate_children_idx:
        if ci in pre_assigned_children:
            continue
        reganimal = str(df_children.at[ci, "reganimal"]).strip() if "reganimal" in df_children.columns else str(ci)
        child_candidates = candidates_for_child.get(ci, [])
        if not child_candidates:
            continue
        for bi, _score in child_candidates:
            father_id = str(df_bulls.at[bi, father_id_col]).strip()
            bvals = bulls_loci[bi]
            total_pairs += 1

            # Child row
            row_child: Dict[str, Any] = {c: "" for c in all_cols}
            row_child["reganimal"] = reganimal
            row_child["father"] = father_id
            row_child["role"] = "child"
            for locus in loci_order:
                c1, c2 = children_loci.get(ci, {}).get(locus, ("", ""))
                row_child[f"{locus}_1"] = c1
                row_child[f"{locus}_2"] = c2
            report_rows.append(row_child)

            # Father row
            row_father: Dict[str, Any] = {c: "" for c in all_cols}
            row_father["reganimal"] = reganimal
            row_father["father"] = father_id
            row_father["role"] = "father"
            for locus in loci_order:
                f1, f2 = bvals.get(locus, ("", ""))
                row_father[f"{locus}_1"] = f1
                row_father[f"{locus}_2"] = f2
            report_rows.append(row_father)

            # Blank separator
            report_rows.append({c: "" for c in all_cols})

    print(f"Паров ребенок-отец для отчета (только новые): {total_pairs}")

    if total_pairs == 0:
        # Print top-1 candidate per child for diagnostics
        print("Не назначено ни одной пары. Диагностика по детям:")
        for ci in candidate_children_idx:
            reganimal = str(df_children.at[ci, "reganimal"]).strip() if "reganimal" in df_children.columns else str(ci)
            overall_idx, (m, mm, cmpd) = best_overall_for_child.get(ci, (None, (-1, -1, -1)))
            if overall_idx is None:
                print(f"  {reganimal}: нет подходящих быков с пересечением локусов")
            else:
                father_id = str(df_bulls.at[overall_idx, father_id_col]).strip()
                print(f"  {reganimal}: лучший {father_id} — совпадений={m}, несовпадений={mm}, сравнивали={cmpd} (пороги: MIN_MATCHED_LOCI={MIN_MATCHED_LOCI}, MAX_MUTATIONS={MAX_MUTATIONS})")

    # ------------------------------
    # Build ALL-children report and stats (ignoring pre-existing fathers)
    # ------------------------------
    def compute_candidates_for_child(ci: int) -> List[Tuple[int, Tuple[int, int, int]]]:
        cvals = children_loci.get(ci, {})
        # If child has too few filled loci, return empty
        child_compared_possible = sum(1 for _l, (a1, a2) in cvals.items() if a1 or a2)
        if child_compared_possible < MIN_MATCHED_LOCI:
            return []
        found: List[Tuple[int, Tuple[int, int, int]]] = []
        for bi in df_bulls.index:
            bvals = bulls_loci[bi]
            matches, mismatches, compared = evaluate_match(cvals, bvals)
            if matches >= MIN_MATCHED_LOCI and mismatches <= MAX_MUTATIONS:
                found.append((bi, (matches, mismatches, compared)))
        # sort by matches desc, mismatches asc, compared desc
        found.sort(key=lambda x: (x[1][0], -x[1][1], x[1][2]), reverse=True)
        return found

    loci_order = [locus for locus, _, _ in child_pairs]
    meta_cols = ["reganimal", "father", "role"]
    locus_cols: List[str] = []
    for locus in loci_order:
        locus_cols.append(f"{locus}_1")
        locus_cols.append(f"{locus}_2")
    all_cols = meta_cols + locus_cols

    report_all_rows: List[Dict[str, Any]] = []
    stats_diff_rows: List[Dict[str, Any]] = []
    stats_original_not_in_candidates: List[Dict[str, Any]] = []

    for ci in df_children.index.tolist():
        reganimal = str(df_children.at[ci, "reganimal"]).strip() if "reganimal" in df_children.columns else str(ci)
        original_father = str(df_children.at[ci, "regotca"]).strip() if "regotca" in df_children.columns else ""

        child_candidates_all = compute_candidates_for_child(ci)
        candidate_father_ids = [str(df_bulls.at[bi, father_id_col]).strip() for bi, _ in child_candidates_all]

        # Build stacked rows for each candidate
        for bi, _score in child_candidates_all:
            father_id = str(df_bulls.at[bi, father_id_col]).strip()
            bvals = bulls_loci[bi]

            row_child: Dict[str, Any] = {c: "" for c in all_cols}
            row_child["reganimal"] = reganimal
            row_child["father"] = father_id
            row_child["role"] = "child"
            for locus in loci_order:
                c1, c2 = children_loci.get(ci, {}).get(locus, ("", ""))
                row_child[f"{locus}_1"] = c1
                row_child[f"{locus}_2"] = c2
            report_all_rows.append(row_child)

            row_father: Dict[str, Any] = {c: "" for c in all_cols}
            row_father["reganimal"] = reganimal
            row_father["father"] = father_id
            row_father["role"] = "father"
            for locus in loci_order:
                f1, f2 = bvals.get(locus, ("", ""))
                row_father[f"{locus}_1"] = f1
                row_father[f"{locus}_2"] = f2
            report_all_rows.append(row_father)

            report_all_rows.append({c: "" for c in all_cols})

        # Stats
        best_found = candidate_father_ids[0] if candidate_father_ids else ""
        if original_father:
            if best_found and best_found != original_father:
                stats_diff_rows.append({
                    "reganimal": reganimal,
                    "original_father": original_father,
                    "best_found_father": best_found,
                })
            if original_father not in set(candidate_father_ids):
                stats_original_not_in_candidates.append({
                    "reganimal": reganimal,
                    "original_father": original_father,
                })

    # Write ALL-children report with highlighting
    report_all_df = pd.DataFrame(report_all_rows, columns=all_cols)
    report_all_path = os.path.join(os.path.dirname(OUTPUT_DB), "assigned_fathers_all_report.xlsx")
    with pd.ExcelWriter(report_all_path, engine="xlsxwriter") as writer:
        report_all_df.to_excel(writer, index=False, sheet_name="report")

        # Highlight mutations
        wb = writer.book
        ws = writer.sheets["report"]
        red_fmt = wb.add_format({"bg_color": "#FFC7CE"})

        col_to_idx = {col_name: idx for idx, col_name in enumerate(all_cols)}
        row_idx_excel = 1
        total_rows = len(report_all_rows)
        while row_idx_excel - 1 < total_rows:
            rr_index = row_idx_excel - 1
            role = str(report_all_df.at[rr_index, "role"]).strip() if rr_index < len(report_all_df) else ""
            if role == "child":
                rr_father_index = rr_index + 1
                if rr_father_index < len(report_all_df) and str(report_all_df.at[rr_father_index, "role"]).strip() == "father":
                    for locus in loci_order:
                        c1 = str(report_all_df.at[rr_index, f"{locus}_1"]).strip() if f"{locus}_1" in report_all_df.columns else ""
                        c2 = str(report_all_df.at[rr_index, f"{locus}_2"]).strip() if f"{locus}_2" in report_all_df.columns else ""
                        f1 = str(report_all_df.at[rr_father_index, f"{locus}_1"]).strip() if f"{locus}_1" in report_all_df.columns else ""
                        f2 = str(report_all_df.at[rr_father_index, f"{locus}_2"]).strip() if f"{locus}_2" in report_all_df.columns else ""
                        child_set = {x for x in [c1, c2] if x}
                        father_set = {x for x in [f1, f2] if x}
                        is_mismatch = bool(child_set) and bool(father_set) and child_set.isdisjoint(father_set)
                        if is_mismatch:
                            for suffix in ("_1", "_2"):
                                col_name = f"{locus}{suffix}"
                                col_index = col_to_idx.get(col_name)
                                if col_index is not None:
                                    val_child = report_all_df.at[rr_index, col_name]
                                    val_father = report_all_df.at[rr_father_index, col_name]
                                    ws.write(row_idx_excel, col_index, val_child, red_fmt)
                                    ws.write(row_idx_excel + 1, col_index, val_father, red_fmt)
                row_idx_excel += 3
                continue
            else:
                row_idx_excel += 1

        # Stats sheet
        stats_start_row = 0
        stats_sheet = wb.add_worksheet("stats")
        # Part 1: best found father differs from original
        stats_sheet.write(stats_start_row, 0, "reganimal")
        stats_sheet.write(stats_start_row, 1, "original_father")
        stats_sheet.write(stats_start_row, 2, "best_found_father")
        r = stats_start_row + 1
        for it in stats_diff_rows:
            stats_sheet.write(r, 0, it.get("reganimal", ""))
            stats_sheet.write(r, 1, it.get("original_father", ""))
            stats_sheet.write(r, 2, it.get("best_found_father", ""))
            r += 1
        # Blank row
        r += 1
        # Part 2: original father not among found candidates
        stats_sheet.write(r, 0, "reganimal")
        stats_sheet.write(r, 1, "original_father")
        r += 1
        for it in stats_original_not_in_candidates:
            stats_sheet.write(r, 0, it.get("reganimal", ""))
            stats_sheet.write(r, 1, it.get("original_father", ""))
            r += 1

    print(f"Полный отчет по всем детям: {report_all_path}")

    report_df = pd.DataFrame(report_rows, columns=all_cols)
    report_path = os.path.join(os.path.dirname(OUTPUT_DB), "assigned_fathers_report.xlsx")
    with pd.ExcelWriter(report_path, engine="xlsxwriter") as writer:
        report_df.to_excel(writer, index=False, sheet_name="report")

        # Highlight mutations (loci where child and father share no alleles) in red fill
        wb = writer.book
        ws = writer.sheets["report"]
        red_fmt = wb.add_format({"bg_color": "#FFC7CE"})

        # Build column index mapping
        col_to_idx = {col_name: idx for idx, col_name in enumerate(all_cols)}

        # Iterate over rows: data starts at row 1 (row 0 is header)
        row_idx_excel = 1
        total_rows = len(report_rows)
        while row_idx_excel - 1 < total_rows:
            # Map excel row to our report_rows index
            rr_index = row_idx_excel - 1
            role = str(report_df.at[rr_index, "role"]).strip() if rr_index < len(report_df) else ""

            if role == "child":
                # father row should be next (child then father then blank)
                rr_father_index = rr_index + 1
                if rr_father_index < len(report_df) and str(report_df.at[rr_father_index, "role"]).strip() == "father":
                    # for each locus, check mismatch
                    for locus in loci_order:
                        c1 = str(report_df.at[rr_index, f"{locus}_1"]).strip() if f"{locus}_1" in report_df.columns else ""
                        c2 = str(report_df.at[rr_index, f"{locus}_2"]).strip() if f"{locus}_2" in report_df.columns else ""
                        f1 = str(report_df.at[rr_father_index, f"{locus}_1"]).strip() if f"{locus}_1" in report_df.columns else ""
                        f2 = str(report_df.at[rr_father_index, f"{locus}_2"]).strip() if f"{locus}_2" in report_df.columns else ""

                        child_set = {x for x in [c1, c2] if x}
                        father_set = {x for x in [f1, f2] if x}
                        is_mismatch = bool(child_set) and bool(father_set) and child_set.isdisjoint(father_set)
                        if is_mismatch:
                            # Color both child's and father's allele cells for this locus
                            for suffix in ("_1", "_2"):
                                col_name = f"{locus}{suffix}"
                                col_index = col_to_idx.get(col_name)
                                if col_index is not None:
                                    # write back same value but with format
                                    val_child = report_df.at[rr_index, col_name]
                                    val_father = report_df.at[rr_father_index, col_name]
                                    ws.write(row_idx_excel, col_index, val_child, red_fmt)
                                    ws.write(row_idx_excel + 1, col_index, val_father, red_fmt)

                # advance: child row, father row, blank row → +3
                row_idx_excel += 3
                continue
            else:
                # not a child row, advance by 1
                row_idx_excel += 1

    print(f"\nГотово. Обновленный файл: {OUTPUT_DB}")
    print(f"Отчет: {report_path}")


if __name__ == "__main__":
    main()
