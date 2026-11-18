"""Acquire the six EPA NRSA 2018–2019 files used by the nutrient audit."""
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
from urllib.request import urlopen

EPA = "https://www.epa.gov"
FILES = {
    "site.csv": ("system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv", "1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3"),
    "site-metadata.txt": ("system/files/other-files/2023-01/NRSA18_19_Site_Information_Metadata.txt", "f9ad8f8b9ee2fbbbd42c1657a123f8d79b7682d599d8806742e0df3298795dd4"),
    "chemistry.csv": ("sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_data.csv", "cb9a3349ae3708834eab533e71423d74c1778bf76aba0737f0ae1bcc8f50baa1"),
    "chemistry-metadata.txt": ("sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_metadata.txt", "c54c4bf4d6b43ca22d4ae4e86cc2d4bec6a14f1610dabee8390fa236112662fd"),
    "population.csv": ("system/files/other-files/2023-11/nrsa1819_data_for_populationestimates.csv", "fca1ea7c5ae8475c93ea2b915cd5829b36ce9da0a68c0bc2753df38214cc1bea"),
    "population-metadata.txt": ("system/files/other-files/2023-11/nrsa1819_metadata_for_populationestimates.txt", "96fd36673a121f4187057f5b1f18a0fcb329b3473d93e78d9a68dd9c865e6601"),
}
OUT = Path(__file__).resolve().parent / "data/raw"


def fetch(name, endpoint, expected):
    url = f"{EPA}/{endpoint}"
    path = OUT / name
    if path.exists() and sha256(path.read_bytes()).hexdigest() == expected:
        return {"file": name, "url": url, "bytes": path.stat().st_size, "sha256": expected}
    with urlopen(url, timeout=60) as response:
        data = response.read(8_000_001)
    if sha256(data).hexdigest() != expected:
        raise ValueError(f"EPA changed {name}; review the new release before updating its fingerprint")
    intermediate = path.with_suffix(".part")
    intermediate.write_bytes(data)
    intermediate.replace(path)
    return {"file": name, "url": url, "bytes": len(data), "sha256": expected}


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(fetch, name, *contents) for name, contents in FILES.items()]
        manifest = [future.result() for future in futures]
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")
