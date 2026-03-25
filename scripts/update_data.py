#!/usr/bin/env python3
"""
GEO CNS Tumor (中枢神经系统肿瘤) 数据更新脚本

默认执行增量更新；当检测到本地数据只有单一年份或只有 bulk RNA-seq 时，
自动回退到全量历史抓取并重建数据文件，避免首次错误初始化后一直缺失历史数据。
"""

import argparse
import json
import os
import re
import time
from collections import Counter
from datetime import datetime, timedelta
from pathlib import Path

import requests
from Bio import Entrez

def resolve_project_root():
    script_path = Path(__file__).resolve()
    candidates = [
        script_path.parent,
        script_path.parent.parent,
        Path.cwd(),
    ]

    for candidate in candidates:
        if (candidate / "data").is_dir() and (candidate / "requirements.txt").is_file():
            return candidate

    return script_path.parent.parent


PROJECT_ROOT = resolve_project_root()
DATA_FILE = PROJECT_ROOT / "data" / "geo_data.json"

NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
MINIMAX_API_KEY = os.environ.get("MINIMAX_API_KEY", "")
HTTP_SESSION = requests.Session()

SEARCH_CONFIG = {
    "keywords": [
        "glioma", "glioblastoma", "GBM", "astrocytoma",
        "oligodendroglioma", "diffuse glioma", "low grade glioma",
        "high grade glioma", "IDH mutant glioma",
        "brain metastasis", "brain metastases",
        "cerebral metastasis", "cerebral metastases",
        "leptomeningeal metastasis", "leptomeningeal metastases",
        "intracranial metastasis", "intracranial metastases",
        "secondary brain tumor", "secondary brain tumour",
        "metastatic brain tumor", "metastatic brain tumour",
        "brain colonization",
        "lung cancer brain metastas", "NSCLC brain metastas",
        "breast cancer brain metastas", "melanoma brain metastas",
        "renal cell carcinoma brain", "colorectal cancer brain",
        "HER2 brain metastas",
        "meningioma",
        "pituitary adenoma", "pituitary tumor", "pituitary tumour",
        "pituitary neuroendocrine tumor",
        "vestibular schwannoma", "acoustic neuroma",
        "medulloblastoma", "ependymoma", "craniopharyngioma",
        "CNS lymphoma", "brain tumor", "brain tumour",
        "central nervous system tumor", "central nervous system tumour",
        "intracranial tumor", "intracranial tumour",
        "DIPG", "diffuse intrinsic pontine glioma",
        "choroid plexus tumor", "neurocytoma",
    ],
    "organisms": ["Homo sapiens", "Mus musculus"],
    "data_types": [
        "Expression profiling by high throughput sequencing",
        "Expression profiling by array",
        "Methylation profiling by array",
        "Methylation profiling by high throughput sequencing",
        "Genome binding/occupancy profiling by high throughput sequencing",
        "Non-coding RNA profiling by high throughput sequencing",
        "Genome variation profiling by high throughput sequencing",
    ],
}

CNS_TERMS = [
    "glioma", "glioblastoma", "gbm", "astrocytoma", "oligodendroglioma",
    "glial tumor", "glial tumour", "gliogenesis",
    "brain metastas", "cerebral metastas", "intracranial metastas",
    "leptomeningeal metastas", "secondary brain tumor", "secondary brain tumour",
    "metastatic brain", "brain coloniz", "brm ",
    "brain metastatic niche", "metastatic niche brain",
    "meningioma",
    "pituitary",
    "schwannoma", "acoustic neuroma",
    "medulloblastoma", "ependymoma", "craniopharyngioma",
    "cns lymphoma", "cns tumor", "cns tumour",
    "intracranial tumor", "intracranial tumour", "intracranial neoplasm",
    "brain cancer", "brain tumor", "brain tumour",
    "neuro-oncol", "neurooncol",
    "dipg", "pontine glioma", "choroid plexus", "neurocytoma",
]

EXCLUDE_TERMS = ["skin fibrosis", "dermal fibrosis", "cardiac fibrosis"]

DATA_TYPE_RULES = [
    ("scRNA-seq", [
        "single-cell rna", "single cell rna", "single-cell transcriptom",
        "single cell transcriptom", "scrna-seq", "scrna seq", "scrnaseq",
        "10x genomics", "10x chromium", "drop-seq", "dropseq",
        "smart-seq", "smartseq", "cel-seq", "celseq", "indrops", "sci-rna",
    ]),
    ("snRNA-seq", [
        "single-nucleus rna", "single nucleus rna", "single-nuclei rna",
        "single nuclei rna", "snrna-seq", "snrna seq", "snrnaseq",
    ]),
    ("spatial transcriptomics", [
        "spatial transcriptom", "spatially resolved transcriptom", "spatial rna",
        "visium", "10x visium", "slide-seq", "slideseq", "merfish",
        "seqfish", "stereo-seq", "dbit-seq", "hdst",
    ]),
    ("bulk RNA-seq", [
        "rna-seq", "rna seq", "rnaseq", "mrna-seq", "mrna seq",
        "transcriptome sequencing", "expression profiling by high throughput sequencing",
    ]),
]

DATA_TYPE_ORDER = [
    "scRNA-seq",
    "snRNA-seq",
    "spatial transcriptomics",
    "bulk RNA-seq",
]

# 只关注这四种数据类型
TARGET_DATA_TYPES = {
    "scRNA-seq",
    "snRNA-seq",
    "spatial transcriptomics",
    "bulk RNA-seq",
}

TRANSCRIPTOME_SPECIFIC_TYPES = {
    "scRNA-seq",
    "snRNA-seq",
    "spatial transcriptomics",
}

TUMOR_TYPE_RULES = [
    ("Glioma - GBM", ["glioblastoma", "gbm "]),
    ("Glioma - Astrocytoma", ["astrocytoma"]),
    ("Glioma - Oligodendroglioma", ["oligodendroglioma"]),
    ("Glioma - DIPG", ["dipg", "diffuse intrinsic pontine glioma"]),
    ("Glioma", ["glioma", "glial tumor", "glial tumour"]),
    ("Brain Metastasis - Lung", [
        "lung cancer brain", "nsclc brain", "lung adenocarcinoma brain",
        "small cell lung cancer brain", "sclc brain",
    ]),
    ("Brain Metastasis - Breast", ["breast cancer brain", "her2 brain", "triple negative brain"]),
    ("Brain Metastasis - Melanoma", [
        "melanoma brain metastas", "melanoma brain coloniz",
        "melanoma cerebral", "melanoma leptomeningeal",
    ]),
    ("Brain Metastasis - Renal", [
        "renal cell carcinoma brain", "renal brain metastas", "kidney cancer brain",
    ]),
    ("Brain Metastasis - Colorectal", ["colorectal cancer brain", "colon cancer brain"]),
    ("Brain Metastasis", [
        "brain metastas", "cerebral metastas", "intracranial metastas",
        "leptomeningeal metastas", "secondary brain tumor", "secondary brain tumour",
        "metastatic brain", "brain coloniz", "brm ",
    ]),
    ("Meningioma", ["meningioma"]),
    ("Pituitary Tumor", [
        "pituitary adenoma", "pituitary tumor", "pituitary tumour",
        "pituitary neuroendocrine", "pituitary carcinoma", "pituitary blastoma",
    ]),
    ("Vestibular Schwannoma", [
        "vestibular schwannoma", "acoustic neuroma", "acoustic schwannoma",
    ]),
    ("Medulloblastoma", ["medulloblastoma"]),
    ("Ependymoma", ["ependymoma"]),
    ("Craniopharyngioma", ["craniopharyngioma"]),
    ("CNS Lymphoma", [
        "cns lymphoma", "primary central nervous system lymphoma", "pcnsl",
    ]),
    ("Choroid Plexus Tumor", ["choroid plexus"]),
    ("Neurocytoma", ["neurocytoma"]),
    ("Other CNS Tumor", [
        "brain tumor", "brain tumour", "brain cancer",
        "central nervous system tumor", "central nervous system tumour",
        "intracranial tumor", "intracranial tumour", "intracranial neoplasm",
    ]),
]


def log(message):
    print(message, flush=True)


def parse_args():
    parser = argparse.ArgumentParser(description="Update GEO CNS tumor dataset cache.")
    parser.add_argument(
        "--full",
        action="store_true",
        help="Ignore incremental mode and rebuild the dataset from all historical GEO results.",
    )
    parser.add_argument(
        "--recent-days",
        type=int,
        default=30,
        help="Incremental mode lookback window in days. Default: 30.",
    )
    parser.add_argument(
        "--skip-ai",
        action="store_true",
        help="Skip AI summary generation even if MINIMAX_API_KEY is set.",
    )
    return parser.parse_args()


def setup_entrez():
    Entrez.email = NCBI_EMAIL
    if NCBI_API_KEY:
        Entrez.api_key = NCBI_API_KEY


def build_query():
    keyword_query = " OR ".join(f'"{kw}"' for kw in SEARCH_CONFIG["keywords"])
    org_query = " OR ".join(f'"{org}"[Organism]' for org in SEARCH_CONFIG["organisms"])
    type_query = " OR ".join(f'"{item}"[DataSet Type]' for item in SEARCH_CONFIG["data_types"])
    return f"({keyword_query}) AND ({org_query}) AND ({type_query})"


def extract_year(date_value):
    match = re.search(r"(19|20)\d{2}", str(date_value or ""))
    return int(match.group()) if match else None


def load_existing_data():
    if DATA_FILE.exists():
        return json.loads(DATA_FILE.read_text(encoding="utf-8"))
    return []


def save_data(data):
    DATA_FILE.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def sort_data_records(data):
    def sort_key(item):
        date_text = str(item.get("Submission_Date", ""))
        try:
            return datetime.strptime(date_text, "%Y/%m/%d")
        except ValueError:
            year = extract_year(date_text) or 0
            return datetime(year=max(year, 1), month=1, day=1)

    return sorted(data, key=sort_key, reverse=True)


def summarize_existing_data(data):
    years = sorted({year for year in (extract_year(row.get("Submission_Date")) for row in data) if year})
    data_types = sorted({(row.get("Data_Type") or "").strip() for row in data if (row.get("Data_Type") or "").strip()})
    return years, data_types


def should_force_full_refresh(data):
    if not data:
        return True, "本地数据为空，需要全量初始化"

    years, data_types = summarize_existing_data(data)
    if len(years) <= 1:
        return True, f"本地数据仅覆盖 {years[0] if years else '未知'} 年，判定为历史数据不完整"

    if set(data_types) == {"bulk RNA-seq"}:
        return True, "本地数据类型全部为 bulk RNA-seq，判定为旧版错误分类结果"

    return False, "保留增量更新"


def is_cns_relevant(record):
    title = record.get("title", "").lower()
    summary = record.get("summary", "").lower()
    combined = f"{title} {summary}"

    if any(term in combined for term in CNS_TERMS):
        return True

    if any(term in combined for term in EXCLUDE_TERMS):
        return False

    return False


def classify_tumor_type(title, summary):
    combined = f"{title} {summary}".lower()
    for standard_name, keywords in TUMOR_TYPE_RULES:
        if any(keyword in combined for keyword in keywords):
            main_type = standard_name.split(" - ")[0] if " - " in standard_name else standard_name
            return main_type, standard_name
    return "Other CNS Tumor", "Other CNS Tumor"


def classify_data_type(title, summary, overall_design, geo_data_type=""):
    combined = " ".join([
        title or "",
        summary or "",
        overall_design or "",
        geo_data_type or "",
    ]).lower()

    matched = []

    has_single_cell_rna_context = (
        ("single-cell" in combined or "single cell" in combined)
        and any(token in combined for token in ["transcriptom", "gene expression", "rna"])
    )
    has_single_nucleus_rna_context = (
        any(token in combined for token in ["single-nucleus", "single nucleus", "single-nuclei", "single nuclei"])
        and any(token in combined for token in ["transcriptom", "gene expression", "rna"])
    )

    if has_single_cell_rna_context:
        matched.append("scRNA-seq")
    if has_single_nucleus_rna_context:
        matched.append("snRNA-seq")

    for label, keywords in DATA_TYPE_RULES:
        if any(keyword in combined for keyword in keywords):
            matched.append(label)

    unique_matches = []
    for label in matched:
        if label not in unique_matches:
            unique_matches.append(label)

    # 只保留目标数据类型
    unique_matches = [label for label in unique_matches if label in TARGET_DATA_TYPES]

    # 如果匹配到更具体的转录组类型（sc/sn/spatial），去掉 bulk RNA-seq
    if "bulk RNA-seq" in unique_matches and any(label in unique_matches for label in TRANSCRIPTOME_SPECIFIC_TYPES):
        unique_matches.remove("bulk RNA-seq")

    ordered_matches = [label for label in DATA_TYPE_ORDER if label in unique_matches]
    if ordered_matches:
        return "; ".join(ordered_matches)
    return ""


def search_geo(full_refresh=False, recent_days=30, max_retries=3):
    query = build_query()
    end_date = datetime.now().strftime("%Y/%m/%d")
    search_params = {
        "db": "gds",
        "term": query,
        "retmax": 10000,
        "usehistory": "y",
    }

    if full_refresh:
        log("执行全量历史抓取：不限制日期范围")
    else:
        start_date = (datetime.now() - timedelta(days=recent_days)).strftime("%Y/%m/%d")
        search_params.update({
            "mindate": start_date,
            "maxdate": end_date,
            "datetype": "pdat",
        })
        log(f"执行增量更新，日期范围: {start_date} - {end_date}")

    log(f"搜索查询: {query[:120]}...")

    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(**search_params)
            results = Entrez.read(handle)
            handle.close()
            return results.get("IdList", [])
        except Exception as exc:
            log(f"搜索失败 (尝试 {attempt + 1}/{max_retries}): {exc}")
            if attempt < max_retries - 1:
                time.sleep(10)
            else:
                return []


def batched(items, size):
    for start in range(0, len(items), size):
        yield items[start:start + size]


def fetch_summaries(id_list, max_retries=3, batch_size=200):
    if not id_list:
        return []

    all_records = []
    for batch in batched(id_list, batch_size):
        for attempt in range(max_retries):
            try:
                handle = Entrez.esummary(db="gds", id=",".join(batch))
                records = Entrez.read(handle)
                handle.close()
                all_records.extend(records)
                log(f"已获取摘要批次: {len(all_records)}/{len(id_list)}")
                break
            except Exception as exc:
                log(f"获取摘要失败 (尝试 {attempt + 1}/{max_retries}): {exc}")
                if attempt < max_retries - 1:
                    time.sleep(10)
                else:
                    log(f"跳过批次: {batch[:3]}... 共 {len(batch)} 条")
        time.sleep(0.1)
    return all_records


def clean_pubmed_ids(pubmed_str):
    if not pubmed_str:
        return ""

    numbers = re.findall(r"IntegerElement\((\d+)", str(pubmed_str))
    if not numbers:
        numbers = re.findall(r"\d+", str(pubmed_str))
    return "; ".join(numbers) if numbers else str(pubmed_str)


def generate_ai_summary(title, summary, data_type):
    if not MINIMAX_API_KEY:
        return ""

    prompt = f"""请用中文为以下GEO数据集生成一个精炼的科研摘要（80-120字）：

标题: {title}
数据类型: {data_type}
研究摘要: {summary[:800]}

请直接输出中文摘要："""

    try:
        response = HTTP_SESSION.post(
            "https://api.minimaxi.com/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {MINIMAX_API_KEY}",
                "Content-Type": "application/json",
            },
            json={
                "model": "MiniMax-M2.1",
                "messages": [{"role": "user", "content": prompt}],
                "max_tokens": 1500,
                "temperature": 0.7,
            },
            timeout=60,
        )
        if response.status_code == 200:
            content = response.json()["choices"][0]["message"]["content"]
            return re.sub(r"<think>.*?</think>", "", content, flags=re.DOTALL).strip()
    except Exception as exc:
        log(f"AI 摘要生成失败: {exc}")
    return ""


def fetch_geo_soft(accession):
    url = (
        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
        "&targ=self&form=text&view=full"
    )
    try:
        response = HTTP_SESSION.get(url, timeout=15)
        if response.status_code != 200:
            return {}

        info = {
            "overall_design": "",
            "contributors": [],
            "lab": "",
            "institute": "",
            "country": "",
            "series_type": [],
        }

        for raw_line in response.text.split("\n"):
            line = raw_line.strip()
            if line.startswith("!Series_overall_design"):
                info["overall_design"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Series_type"):
                value = line.split("=", 1)[1].strip()
                if value and value not in info["series_type"]:
                    info["series_type"].append(value)
            elif line.startswith("!Series_contributor"):
                contributor = line.split("=", 1)[1].strip()
                parts = [part.strip() for part in contributor.split(",") if part.strip()]
                if len(parts) >= 2:
                    formatted = f"{parts[-1]} {parts[0]}".strip()
                    if formatted and formatted not in info["contributors"]:
                        info["contributors"].append(formatted)
            elif line.startswith("!Series_contact_laboratory"):
                info["lab"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Series_contact_institute"):
                info["institute"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Series_contact_country"):
                info["country"] = line.split("=", 1)[1].strip()

        return info
    except Exception as exc:
        log(f"    获取 SOFT 信息失败: {exc}")
        return {}


def parse_record(record, existing_entry=None, skip_ai=False):
    accession = record.get("Accession", "")
    if not accession.startswith("GSE"):
        return None

    title = record.get("title", "")
    summary = record.get("summary", "")
    pubmed_ids = record.get("PubMedIds", [])
    pubmed_str = clean_pubmed_ids("; ".join(str(item) for item in pubmed_ids) if pubmed_ids else "")

    soft_info = fetch_geo_soft(accession)

    geo_data_type = "; ".join(soft_info.get("series_type", []))
    overall_design = soft_info.get("overall_design", "")
    data_type = classify_data_type(title, summary, overall_design, geo_data_type)

    # 跳过不属于目标数据类型的记录
    if not data_type:
        return None

    tumor_main, tumor_sub = classify_tumor_type(title, summary)

    existing_ai_summary = ""
    if existing_entry:
        existing_ai_summary = existing_entry.get("AI_Summary_CN") or existing_entry.get("AI_Summary") or ""

    ai_summary = existing_ai_summary
    if not ai_summary and not skip_ai:
        ai_summary = generate_ai_summary(title, summary, data_type)
        if ai_summary:
            time.sleep(1)

    return {
        "Accession": accession,
        "Title": title,
        "Organism": record.get("taxon", ""),
        "Tumor_Type": tumor_main,
        "Tumor_Subtype": tumor_sub,
        "Data_Type": data_type,
        "Sample_Count": record.get("n_samples", 0),
        "Platform": record.get("GPL", ""),
        "Country": soft_info.get("country", ""),
        "Lab": soft_info.get("lab", ""),
        "Institute": soft_info.get("institute", ""),
        "Contributors": "; ".join(soft_info.get("contributors", [])),
        "PubMed_IDs": pubmed_str,
        "Supplementary_Size": existing_entry.get("Supplementary_Size", "N/A") if existing_entry else "N/A",
        "Summary": summary,
        "Overall_Design": overall_design,
        "AI_Summary_CN": ai_summary,
        "AI_Summary": ai_summary,
        "GEO_Link": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
        "Submission_Date": record.get("PDAT", ""),
    }


def process_full_refresh(summaries, existing_by_accession, skip_ai=False):
    rebuilt_data = []
    skipped_irrelevant = 0

    total_records = len(summaries)
    for index, record in enumerate(summaries, start=1):
        accession = record.get("Accession", "")
        if index == 1 or index % 25 == 0:
            log(f"全量重建进度: {index}/{total_records}")
        if not accession.startswith("GSE"):
            continue

        if not is_cns_relevant(record):
            skipped_irrelevant += 1
            log(f"  跳过 (非 CNS 肿瘤相关): {accession}")
            continue

        parsed = parse_record(record, existing_entry=existing_by_accession.get(accession), skip_ai=skip_ai)
        if parsed:
            rebuilt_data.append(parsed)
            log(f"  刷新: {accession} -> {parsed['Data_Type']}")

    if skipped_irrelevant:
        log(f"跳过 {skipped_irrelevant} 条非 CNS 肿瘤相关数据集")

    return sort_data_records(rebuilt_data)


def process_incremental_update(existing_data, summaries, skip_ai=False):
    existing_accessions = {row["Accession"] for row in existing_data}
    new_count = 0
    skipped_irrelevant = 0

    for record in summaries:
        accession = record.get("Accession", "")
        if accession in existing_accessions or not accession.startswith("GSE"):
            continue

        if not is_cns_relevant(record):
            skipped_irrelevant += 1
            log(f"  跳过 (非 CNS 肿瘤相关): {accession}")
            continue

        parsed = parse_record(record, skip_ai=skip_ai)
        if parsed:
            existing_data.append(parsed)
            existing_accessions.add(accession)
            new_count += 1
            log(f"  新增: {accession} -> {parsed['Data_Type']}")

    if skipped_irrelevant:
        log(f"跳过 {skipped_irrelevant} 条非 CNS 肿瘤相关数据集")

    updated = sort_data_records(existing_data)
    return updated, new_count


def print_dataset_summary(data):
    years = Counter()
    data_types = Counter()
    for row in data:
        year = extract_year(row.get("Submission_Date"))
        if year:
            years[year] += 1
        data_types[row.get("Data_Type", "Other")] += 1

    log(f"总数据集: {len(data)}")
    log(f"年份覆盖: {dict(sorted(years.items()))}")
    log(f"数据类型 Top 10: {data_types.most_common(10)}")


def main():
    args = parse_args()
    log(f"开始更新 CNS Tumor 数据 - {datetime.now()}")

    if not NCBI_EMAIL:
        log("错误: 未设置 NCBI_EMAIL")
        return

    setup_entrez()

    existing_data = load_existing_data()
    existing_by_accession = {row.get("Accession", ""): row for row in existing_data}
    log(f"现有数据集: {len(existing_data)}")

    auto_full_refresh, refresh_reason = should_force_full_refresh(existing_data)
    full_refresh = args.full or auto_full_refresh
    log(f"更新模式: {'全量重建' if full_refresh else '增量更新'} ({refresh_reason if auto_full_refresh and not args.full else '手动指定全量重建' if args.full else '本地数据完整'})")

    id_list = search_geo(full_refresh=full_refresh, recent_days=args.recent_days)
    log(f"搜索到: {len(id_list)} 条记录")
    if not id_list:
        log("没有可处理的数据")
        return

    summaries = fetch_summaries(id_list)
    log(f"获取摘要完成: {len(summaries)} 条")

    if full_refresh:
        updated_data = process_full_refresh(summaries, existing_by_accession, skip_ai=args.skip_ai)
        save_data(updated_data)
        log("全量重建完成")
    else:
        updated_data, new_count = process_incremental_update(existing_data, summaries, skip_ai=args.skip_ai)
        if new_count > 0:
            save_data(updated_data)
            log(f"增量更新完成，新增 {new_count} 条")
        else:
            log("没有新数据需要添加")
            updated_data = sort_data_records(existing_data)

    print_dataset_summary(updated_data)


if __name__ == "__main__":
    main()
