import os
import re
import pandas as pd

# -------------------------
# Настройки
# -------------------------
raw_folder = r"C:\Users\user\Desktop\genetic\zrya_raw"
output_folder = r"C:\Users\user\Desktop\genetic\zrya_processed"
os.makedirs(output_folder, exist_ok=True)

# -------------------------
# Вспомогательные функции
# -------------------------
def clean_hoz_name(filename: str) -> str:
    """Очистить имя файла -> читаемое имя хозяйства (убираем 'партия' и т.п.)"""
    name = os.path.splitext(filename)[0]
    name = re.sub(r"\s*\(.*партия.*\)\s*", "", name, flags=re.IGNORECASE)
    return name.strip()

def cell_text(v):
    """Преобразовать значение в строку без лишних пробелов"""
    if pd.isna(v):
        return ""
    return str(v).strip()

def split_ids(raw_id: str) -> list:
    """Разбить строку идентификатора на отдельные ID только по явным разделителям '/', ',', ';' или '\\'.
    Пробелы считаем частью идентификатора (например: '9061 Маяк Рр' — один ID)."""
    s = (raw_id or "").strip()
    if not s:
        return []
    # split ONLY on explicit separators, not on whitespace
    parts = re.split(r"[/,;\\]+", s)
    return [p.strip() for p in parts if p.strip()]

# -------------------------
# Фиксированный список локусов (ровно 16 для пар аллелей) в заданном порядке.
# AMEL исключаем из пар, чтобы получить 16*2=32 столбца. При необходимости AMEL можно добавить отдельно.
# -------------------------
CANONICAL_LOCI_ORDER = [
    "TGLA227","BM2113","TGLA53","ETH10","SPS115","TGLA122",
    "INRA23","TGLA126","BM1818","ETH225","BM1824","CSRM60",
    "CSSM43","ETH3","ILST006","HAUT27","AMEL"
]
FIXED_LOCI = CANONICAL_LOCI_ORDER[:16]

# -------------------------
# Основной проход по файлам
# -------------------------
all_data = []
hoz_mapping = {}
hoz_counter = 1
errors = []
father_registry = {}
registry_loci_set = set()

for fname in sorted(os.listdir(raw_folder)):
    if not (fname.lower().endswith(".xlsx") or fname.lower().endswith(".xls")):
        continue

    path = os.path.join(raw_folder, fname)
    hoz_name = clean_hoz_name(fname)
    if hoz_name not in hoz_mapping:
        hoz_mapping[hoz_name] = hoz_counter
        hoz_counter += 1
    nomhoz = hoz_mapping[hoz_name]

    print(f"Обрабатываю файл: {fname} (хоз: {hoz_name} → {nomhoz})")

    try:
        # Читаем все листы для устойчивости к разметке
        sheets = pd.read_excel(path, header=None, dtype=object, sheet_name=None)
    except Exception as e:
        errors.append(f"Ошибка чтения {fname}: {e}")
        print(errors[-1])
        continue

    file_animals = 0

    for sheet_name, df in sheets.items():
        if df is None or df.empty:
            print(f"  Лист '{sheet_name}': пустой")
            continue

        # Используем фиксированный порядок локусов (без AMEL): ровно 16
        loci = FIXED_LOCI
        print(f"  Лист '{sheet_name}': используем фиксированные локусы = {len(loci)}")

        nrows = df.shape[0]
        i = 0
        while i < nrows:
            # Ищем слово 'потомок' в строке (не только в столбце B)
            row_text_joined = " ".join([cell_text(x).lower() for x in df.iloc[i, :].tolist()])
            if "потомок" in row_text_joined:
                # Проверка на границы для блока из 6 строк
                if i + 5 >= nrows:
                    errors.append(f"{fname} / лист '{sheet_name}': неполный блок потомка начиная со строки {i+1}")
                    print(f"  Лист '{sheet_name}': неполный блок потомка (i={i})")
                    i += 1
                    continue

                reganimal = cell_text(df.iloc[i, 2])

                # потомок
                values_child_1 = df.iloc[i, 3:3+len(loci)].tolist()
                values_child_2 = df.iloc[i+1, 3:3+len(loci)].tolist()

                # мать (идёт перед отцом)
                regmateri = cell_text(df.iloc[i+2, 2])
                values_mother_1 = df.iloc[i+2, 3:3+len(loci)].tolist()
                values_mother_2 = df.iloc[i+3, 3:3+len(loci)].tolist()

                # отец
                regotca = cell_text(df.iloc[i+4, 2])
                values_father_1 = df.iloc[i+4, 3:3+len(loci)].tolist()
                values_father_2 = df.iloc[i+5, 3:3+len(loci)].tolist()

                # статус: ищем по всей строке потомка
                status_row_text = row_text_joined
                if "по отцу и по матери" in status_row_text:
                    status = 1
                elif "по отцу" in status_row_text:
                    status = 2
                elif "по матери" in status_row_text:
                    status = 3
                else:
                    status = 0  # не указан

                rec = {
                    "nomanimal": len(all_data) + 1,
                    "reganimal": reganimal,
                    "nomhoz": nomhoz,
                    "regotca": regotca,
                    "regmateri": regmateri,
                    "status": status
                }

                # потомок
                for j, locus in enumerate(loci):
                    rec[f"1_{locus}"] = cell_text(values_child_1[j]) if j < len(values_child_1) else ""
                    rec[f"2_{locus}"] = cell_text(values_child_2[j]) if j < len(values_child_2) else ""

                # отец
                for j, locus in enumerate(loci):
                    rec[f"1_{locus}_otca"] = cell_text(values_father_1[j]) if j < len(values_father_1) else ""
                    rec[f"2_{locus}_otca"] = cell_text(values_father_2[j]) if j < len(values_father_2) else ""

                # мать
                for j, locus in enumerate(loci):
                    rec[f"1_{locus}_materi"] = cell_text(values_mother_1[j]) if j < len(values_mother_1) else ""
                    rec[f"2_{locus}_materi"] = cell_text(values_mother_2[j]) if j < len(values_mother_2) else ""

                all_data.append(rec)
                # accumulate father registry (возможны множественные ID через '/', ',', пробел)
                for fid in split_ids(regotca):
                    entry = father_registry.setdefault(fid, {})
                    for j, locus in enumerate(loci):
                        f1 = cell_text(values_father_1[j]) if j < len(values_father_1) else ""
                        f2 = cell_text(values_father_2[j]) if j < len(values_father_2) else ""
                        if f1 or f2:
                            registry_loci_set.add(locus)
                            if f"1_{locus}" not in entry or not entry[f"1_{locus}"]:
                                entry[f"1_{locus}"] = f1
                            if f"2_{locus}" not in entry or not entry[f"2_{locus}"]:
                                entry[f"2_{locus}"] = f2
                file_animals += 1
                i += 6
                continue

            i += 1

    print(f"  Найдено животных в файле: {file_animals}")

# -------------------------
# Сохранение
# -------------------------
df_all = pd.DataFrame(all_data)
out_csv = os.path.join(output_folder, "genotypes_unified.csv")
df_all.to_csv(out_csv, index=False, sep=";", encoding="utf-8-sig")

hoz_list = pd.DataFrame([{"nomhoz": v, "name_hoz": k} for k, v in hoz_mapping.items()])
hoz_csv = os.path.join(output_folder, "hoz_list.csv")
hoz_list.to_csv(hoz_csv, index=False, sep=";", encoding="utf-8-sig")

# fathers registry CSV (строго по фиксированным локусам)
if father_registry:
    loci_order = FIXED_LOCI
    cols = ["Идентификационный номер"]
    for l in loci_order:
        cols.append(f"1_{l}")
        cols.append(f"2_{l}")
    rows = []
    for fid, data in father_registry.items():
        row = {c: "" for c in cols}
        row["Идентификационный номер"] = fid
        for l in loci_order:
            row[f"1_{l}"] = data.get(f"1_{l}", "")
            row[f"2_{l}"] = data.get(f"2_{l}", "")
        rows.append(row)
    fathers_csv = os.path.join(output_folder, "fathers_registry.csv")
    pd.DataFrame(rows, columns=cols).to_csv(fathers_csv, index=False, sep=";", encoding="utf-8-sig")

err_log = os.path.join(output_folder, "processing_errors.txt")
with open(err_log, "w", encoding="utf-8") as f:
    f.write("\n".join(errors))

print(f"\nГотово! Записано {len(df_all)} животных.")
print(f"Файлы: {out_csv}, {hoz_csv}")
if father_registry:
    print("Сформирован реестр отцов: fathers_registry.csv")
