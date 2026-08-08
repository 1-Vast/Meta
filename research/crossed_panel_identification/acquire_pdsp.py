import urllib.request, ssl, hashlib, os, sys

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

out_dir = r"D:\MetaSieve\dataset\raw\crossed_panels\pdsp_kidb"
os.makedirs(out_dir, exist_ok=True)

url = "https://pdsp.unc.edu/databases/kiDownload/download.php"
req = urllib.request.Request(url, headers=UA)
r = urllib.request.urlopen(req, timeout=300, context=ctx)
print("status", r.status, r.headers.get("Content-Type"), r.headers.get("Content-Disposition"))
data = r.read()
print("bytes", len(data))
h = hashlib.sha256(data).hexdigest()
print("sha256", h)
path = os.path.join(out_dir, "KiDatabase.csv")
with open(path, "wb") as f:
    f.write(data)
print("wrote", path)
print("--- head ---")
sys.stdout.write(data[:1500].decode("utf-8", "ignore"))
