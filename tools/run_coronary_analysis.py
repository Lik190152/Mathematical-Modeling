import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from coronary_analysis.pipeline import run_full_analysis


def main():
    root = ROOT
    out_dir = root / "analysis_results"
    dataset = run_full_analysis(root, out_dir)
    print(f"analyzed_masks={len(dataset['results'])}")
    print(f"output_dir={out_dir}")


if __name__ == "__main__":
    main()
