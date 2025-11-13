from tempsloading import (
	extract_headers,
	build_dataframe,
	FILE_321,
	FILE_252, 
	FILE_201,
	FILE_214
)
import pandas as pd

headers = extract_headers(FILE_321)
print("Headers from FILE_321:")
print(headers)

df_321 = build_dataframe(FILE_321, min_date="20190101")
print("\nFiltered (dates > 20190101) DataFrame columns:", list(df_321.columns))
print("Filtered shape:", df_321.shape)
print("\nPreview:")
print(df_321.head(5))
mean_321 = pd.to_numeric(df_321['TZ'], errors='coerce').mean()
print(f'The mean TZ is {float(mean_321)}')


df_252 = build_dataframe(FILE_252, min_date="20190101")
print("DataFrame columns:", list(df_252.columns))
print("DataFrame shape:", df_252.shape)
print("Preview:")
print(df_252.head(5))
mean_252 = pd.to_numeric(df_252['TZ'], errors='coerce').mean()
print(f'The mean TZ is {float(mean_252)}')


df_214 = build_dataframe(FILE_214, min_date="20190101")
print("\nFiltered (dates > 20190101) DataFrame columns:", list(df_214.columns))
print("Filtered shape:", df_214.shape)
print("\nPreview:")
print(df_214.head(5))
mean_214 = pd.to_numeric(df_214['TZ'], errors='coerce').mean()
print(f'The mean TZ is {float(mean_214)}')

df_201 = build_dataframe(FILE_201, min_date="20190101")
print("\nFiltered (dates > 20190101) DataFrame columns:", list(df_201.columns))
print("Filtered shape:", df_201.shape)
print("\nPreview:")
print(df_201.head(5))
mean_201 = pd.to_numeric(df_201['TZ'], errors='coerce').mean()
print(f'The mean TZ is {float(mean_201)}')
