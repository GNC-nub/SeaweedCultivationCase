
from pathlib import Path
from tempsloading import (
    build_dataframe,
    FILE_321,
    FILE_252,
)

def export_station_csv(path: Path, min_date: str | None = None) -> Path:
    """Build dataframe from a KNMI station text file and write a CSV next to it.

    Returns the CSV path.
    """
    df = build_dataframe(path, min_date=min_date)
    csv_path = path.with_suffix(".csv")
    df.to_csv(csv_path, index=False)
    print(f"Saved {csv_path} with shape {df.shape}")
    return csv_path


if __name__ == "__main__":
    min_date = "20190101"
    export_station_csv(FILE_321, min_date=min_date)
    export_station_csv(FILE_252, min_date=min_date)
