# -*- coding: utf-8 -*-
"""全局配置：阈值常量 + 字段名映射。

表名与列名已于 2026-08-02 按 IDA 实际下载文件核对（All_Subjects_*_02Aug2026）。
基线定义：每张表按各 RID 的日期最早访问（各表 VISCODE 体系不统一，不依赖 'bl'）。
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"
SEED = 42

# ---------- 阳性阈值 ----------
# p_tau217_pg_ml: 由 scripts/calibrate.py 输出 Youden 阈值，人工确认后填入
# abeta_ratio: 由 calibrate.py 输出（备用主轴）
THRESHOLDS = {
    "p_tau217_pg_ml": 0.183,  # calibrate.py Youden 阈值 (AUC 0.885) 2026-08-02 自动标定，待复核
    "abeta_ratio": 0.085,     # calibrate.py Youden 阈值 (AUC 0.790, 逆向标志物) 备用主轴
    "pet_centiloid": 19.0,    # 备用 PET 阈值（CENTILOIDS），主用官方 AMYLOID_STATUS
}

# ---------- 文件与字段映射 ----------
# 2026-09-02 切换数据源：IDA 单表 CSV → ADNIMERGE2 总和数据
# （scripts/export_adnimerge2.R 从 ADNI_MERGE/ADNIMERGE2/data/*.rda 导出到 data/raw/merged/）
# 列名已与旧 IDA 表核对一致，FIELDS 无需改动。
FILES = {
    "ptdemog": "merged/PTDEMOG.csv",
    "dxsum": "merged/DXSUM.csv",
    "cdr": "merged/CDR.csv",
    "mmse": "merged/MMSE.csv",
    "adas": "merged/ADAS.csv",
    "medhist": "merged/MEDHIST.csv",
    "inithealth": "merged/INITHEALTH.csv",
    "plasma": "merged/UPENN_PLASMA_FUJIREBIO_QUANTERIX.csv",
    "pet": "merged/UCBERKELEY_AMY_6MM.csv",
    "tau": "merged/UCBERKELEY_TAU_6MM.csv",
}

FIELDS = {
    "rid": "RID",
    "date": {  # 各表日期列
        "ptdemog": "VISDATE", "dxsum": "EXAMDATE", "cdr": "VISDATE",
        "mmse": "VISDATE", "adas": "VISDATE", "plasma": "EXAMDATE",
        "pet": "SCANDATE", "tau": "SCANDATE",
    },
    "pt_dob_yy": "PTDOBYY",   # 出生日期（算年龄用）
    "gender": "PTGENDER", "educat": "PTEDUCAT",
    "dx": "DIAGNOSIS",        # 1=CN 2=MCI 3=AD (ADNI 标准编码，含 10=ADNI4 新码)
    "cdrsb": "CDRSB", "mmse": "MMSCORE", "adas13": "TOTAL13",
    "p_tau217": "pT217_F", "abeta42_ab40": "AB42_AB40_F",
    "pet_status": "AMYLOID_STATUS", "pet_suvr": "SUMMARY_SUVR",
    "centiloids": "CENTILOIDS", "tracer": "TRACER",
    # tau PET（UC Berkeley 6mm，META_TEMPORAL 为 AD 标准区）
    "tau_suvr": "META_TEMPORAL_SUVR", "qc": "qc_flag",
}

# 合并症：MEDHIST 仅覆盖 ADNI1/GO/2 且为系统标记（MH4CARD 心血管、MH9ENDO 内分泌）
# 作为混杂检查的代理变量（第二轮启用；第一轮跳过并在报告中注明）
COMORBIDITY_PROXY = {"MH4CARD": "CARDIO", "MH9ENDO": "ENDO"}


def normalize_dx(s):
    """DXSUM DIAGNOSIS 统一为数字编码 1=CN, 2=MCI, 3=AD。
    ADNIMERGE2 用 'CN'/'MCI'/'Dementia' 字符串，老 IDA 是数字 1/2/3；两者兼容。"""
    import pandas as pd
    mapped = s.map({"CN": 1, "MCI": 2, "Dementia": 3, "AD": 3})
    numeric = pd.to_numeric(s, errors="coerce")
    return mapped.where(mapped.notna(), numeric)


def normalize_yes_no(s):
    """'Yes'/'No' 字符串 → 1/0。ADNIMERGE2 把 MEDHIST 系统标记（MH4CARD 等）
    重编码为 'Yes'/'No'，老 IDA CSV 是数字 0/1；两者兼容。"""
    import pandas as pd
    mapped = s.map({"Yes": 1, "No": 0})
    numeric = pd.to_numeric(s, errors="coerce")
    return mapped.where(mapped.notna(), numeric)


def normalize_gender(s):
    """PTGENDER 统一为数字编码 1=Male, 2=Female。
    ADNIMERGE2 打包时把 PTGENDER 重编码为 'Male'/'Female' 字符串，
    老 IDA CSV 是数字 1/2；此处两者兼容，下游回归按数字编码消费。"""
    import pandas as pd
    mapped = s.map({"Male": 1, "Female": 2})
    numeric = pd.to_numeric(s, errors="coerce")
    return mapped.where(mapped.notna(), numeric)
