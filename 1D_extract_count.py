import os
import polars as pl


def get_file_paths(directory: str) -> list[str]:
    file_paths = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file == "quant.sf":
                file_paths.append(os.path.join(root, file))
    return file_paths


def concat_dataframes(paths: list[str]) -> pl.DataFrame:
    df_res = None
    for path in paths:
        _, _, colname = path.replace("/quant.sf", "").rpartition("/")
        df = pl.read_csv(path, separator="\t")
        if df_res is None:
            df_res = df.select(["Name"])
        df_res = df_res.join(df.select(["Name", "TPM"]).rename({"TPM": colname}), on="Name", how="left")
    return df_res

paths = get_file_paths("/home/melanie/net/virus_melanie/maize/results/t4_all_v1/") # change me
df = concat_dataframes(paths)
df.write_csv("/home/melanie/net/virus_melanie/maize/results/t4_all_v1/MM20250403_TPM_basic_salmon.csv") # change me
