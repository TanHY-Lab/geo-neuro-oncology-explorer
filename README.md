# GEO CNS Tumor Explorer

Static GEO dataset explorer for CNS tumors, designed to run on GitHub Pages with scheduled data updates via GitHub Actions.

## Required GitHub setup

Add these repository secrets in `Settings -> Secrets and variables -> Actions`:

- `NCBI_EMAIL` (required)
- `NCBI_API_KEY` (optional)
- `MINIMAX_API_KEY` (optional)

## Workflows

- `Daily Data Update`
  - runs every day
  - executes `python update_data.py --skip-ai`
  - commits updated `data/geo_data.json` back to `main`

- `Deploy GitHub Pages`
  - deploys the static site from the repository contents
  - runs automatically on pushes to `main`

## GitHub Pages

In `Settings -> Pages`:

1. Set `Source` to `GitHub Actions`.
2. Save the setting.

After that, every successful data-update commit will trigger a Pages deployment automatically.

## Local run

```bash
pip install -r requirements.txt
export NCBI_EMAIL="your-email@example.com"
python update_data.py --skip-ai
```

## Notes

- The first successful run may perform a full historical rebuild if the current `data/geo_data.json` is incomplete.
- The frontend fetches `data/geo_data.json` with cache-busting, so GitHub Pages visitors should see fresh data after deployment.
