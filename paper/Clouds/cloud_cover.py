import xarray as xr
import os

# -----------------------------
# User Inputs
# -----------------------------
data_folder = "cm_saf_cfc_new"
cloud_threshold = 40   # percent cloud fraction considered "clear"

# -----------------------------
# Ground station coordinates
# -----------------------------
stations = {
    "Namea": (37.84, 22.62),
    "Heraklion": (35.33, 25.13),
    "Teide Observatory": (28.30, -16.51),
    "FOGATA": (37.09, -2.36),
    "LaBoT": (52.85, 10.13),
    "OGS NBB": (48.08, 11.64),
    "OGS NSG": (53.33, 13.07),
    "OGSOP-NG": (48.08, 11.28),
}

# -----------------------------
# Step 1: List all .nc files
# -----------------------------
nc_files = [
    os.path.join(data_folder, f)
    for f in os.listdir(data_folder)
    if f.endswith(".nc")
]

if not nc_files:
    raise FileNotFoundError(f"No .nc files found in {data_folder}")

# -----------------------------
# Step 2: Open dataset
# -----------------------------
ds = xr.open_mfdataset(nc_files, combine="by_coords")
print(ds["cfc"])
print(ds["cfc"].attrs)
results = []

# -----------------------------
# Step 3: Loop through stations
# -----------------------------
for name, (lat, lon) in stations.items():

    station = ds.sel(lat=lat, lon=lon, method="nearest")

    cfc = station["cfc"]  # Percent Cloud Fraction

    clear = cfc < cloud_threshold

    availability = clear.mean().compute().item()

    results.append((name, lat, lon, availability))

# -----------------------------
# Step 4: Print results
# -----------------------------
print("\nMean Optical Availability (CFC < {}%)".format(cloud_threshold))
print("-----------------------------------------------------")

for name, lat, lon, availability in results:
    print(f"{name:15s} ({lat:6.2f}, {lon:7.2f})  ->  {availability*100:6.2f}%")
