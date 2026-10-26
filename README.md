# Techie Lab — Fixed Deployment Package

An educational Streamlit cheminformatics app for exploring molecular structures and descriptors.

## Files
- `app.py` — Streamlit application
- `requirements.txt` — Python dependencies
- `packages.txt` — Linux system libraries needed by RDKit drawing components
- `README.md` — setup notes

## Deploy to Streamlit Community Cloud
1. Upload these files to the root of your GitHub repository.
2. In Streamlit Community Cloud, set the main file path to `app.py`.
3. Let the app rebuild after committing the changes.

The `packages.txt` file addresses the missing `libXrender.so.1` system-library error seen in some deployments. If an import error persists, check the full Streamlit Cloud logs.

## Run locally
Use Python 3.10 or 3.11 if possible:
```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

This is an updated project package with a compatible simplified app implementation; review it before replacing files in an existing repository.
